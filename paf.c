#define _DEFAULT_SOURCE
#define _POSIX_C_SOURCE 200809L

#include "paf.h"

#include <ctype.h>
#include <inttypes.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

char *PAF_PATH = NULL, *PAF_REF_FAI = NULL, *PAF_READ_FAI = NULL;
int   PAF_HISTORY = 0, PAF_KEEP_LAS = 0, PAF_STATS = 0;

static FILE     *output;
static Paf_Name *query_names, *ref_names;

struct Paf_Source
  { void      *idx, *hdr;
    size_t     idx_size, hdr_size;
    size_t     idx_released, hdr_released;
    DAZZ_READ *raw;
    int        n, cursor, dam, file;
    DAZZ_STUB *stub;
    FILE      *fai;
    char      *fai_line, *last_name;
    size_t     fai_cap;
    int64_t    last_len;
  };

static void die(const char *why)
{ fprintf(stderr,"himapper PAF: %s\n",why);
  exit(1);
}

static void *checked(void *p)
{ if (!p)
    die("allocation failed");
  return p;
}

static void *map_file(const char *path, size_t *size)
{ FILE       *f = fopen(path,"rb");
  struct stat st;

  if (!f || fstat(fileno(f),&st))
    die(path);
  *size = st.st_size;

  void *p = mmap(NULL,*size,PROT_READ,MAP_PRIVATE,fileno(f),0);

  fclose(f);
  if (p == MAP_FAILED)
    die(path);
  return p;
}

Paf_Source *Paf_Open_Source(DAZZ_DB *db, const char *stub, int dam, const char *fai)
{ Paf_Source *s = checked(calloc(1,sizeof(*s)));
  char       *p = checked(malloc(strlen(db->path)+8));

  sprintf(p,"%s.idx",db->path);
  s->idx = map_file(p,&s->idx_size);
  if (s->idx_size < sizeof(DAZZ_DB))
    die("truncated DB index");
  s->n = ((DAZZ_DB *) s->idx)->ureads;
  if (s->idx_size < sizeof(DAZZ_DB)+(size_t) s->n*sizeof(DAZZ_READ))
    die("truncated DB records");
  s->raw = (DAZZ_READ *) ((char *) s->idx+sizeof(DAZZ_DB));
  s->dam = dam;
  if (dam)
    { if (!fai)
        die("DAM PAF requires --paf-ref-fai / --paf-read-fai to preserve original FASTA lengths (including terminal Ns)");
      sprintf(p,"%s.hdr",db->path);
      s->hdr = map_file(p,&s->hdr_size);
      s->fai = fopen(fai,"r");
      if (!s->fai)
        die(fai);
    }
  else
    { s->stub = Read_DB_Stub((char *) stub,DB_STUB_NREADS | DB_STUB_PROLOGS);
      if (!s->stub)
        die("cannot read DB naming metadata");
    }
  free(p);
  return s;
}

static int64_t original_length(Paf_Source *s, const char *name)
{ if (s->last_name && !strcmp(name,s->last_name))
    return s->last_len;
  while (getline(&s->fai_line,&s->fai_cap,s->fai) >= 0)
    { char *tab = strchr(s->fai_line,'\t'), *end;

      if (!tab)
        die("invalid FASTA index");
      *tab = 0;

      int64_t len = strtoll(tab+1,&end,10);

      if (*end != '\t' || len < 1)
        die("invalid FASTA index length");
      if (strcmp(name,s->fai_line))
        continue;
      free(s->last_name);
      s->last_name = checked(strdup(name));
      s->last_len  = len;
      return len;
    }
  fprintf(stderr,"himapper PAF: FASTA index entry missing/out of DB input order: %s\n",name);
  exit(1);
}

static void release_prefix(void *base, size_t consumed, size_t *released)
{ size_t page = (size_t) sysconf(_SC_PAGESIZE), end = consumed/page*page;

  if (end > *released)
    { /* Metadata is visited once in DB order. Keep resident memory bounded by
         the current block even when the original read DB spans a whole genome. */
      madvise((char *) base+*released,end-*released,MADV_DONTNEED);
      *released = end;
    }
}

Paf_Name *Paf_Read_Names(Paf_Source *s, DAZZ_DB *db)
{ Paf_Name *out = checked(calloc((size_t) db->nreads+1,sizeof(*out)));
  size_t    header_end = s->hdr_released;

  for (int i = 0; i < db->nreads; ++i)
    { DAZZ_READ *r = db->reads+i, *raw;

      for (;;)
        { if (s->cursor >= s->n)
            die("DB metadata and mapped reads disagree");
          raw = s->raw+s->cursor++;
          if (raw->origin == r->origin && raw->rlen == r->rlen && raw->fpulse == r->fpulse &&
              (!s->dam || raw->coff == r->coff))
            break;
        }
      if (s->dam)
        { if (raw->coff < 0 || (uint64_t) raw->coff >= s->hdr_size)
            die("invalid DAM header offset");

          const char *start = (char *) s->hdr+raw->coff, *end = (char *) s->hdr+s->hdr_size;

          if (start < end && *start == '>')
            ++start;

          const char *p = start;

          while (p < end && !isspace((unsigned char) *p) && *p)
            ++p;
          if (p == start)
            die("empty sequence name");
          out[i].name = checked(strndup(start,p-start));
          header_end = p-(char *) s->hdr;
          out[i].length = original_length(s,out[i].name);
          out[i].offset = r->fpulse;
          if (out[i].length < (int64_t) r->fpulse+r->rlen)
            die("FASTA index length disagrees with DAM");
        }
      else
        { while (s->file+1 < s->stub->nfiles && s->cursor-1 >= s->stub->nreads[s->file])
            ++s->file;

          const char *prefix = s->stub->prolog[s->file];

          out[i].name = checked(malloc(strlen(prefix)+80));
          sprintf(out[i].name,"%s/%d/%d_%d",prefix,r->origin,r->fpulse,r->fpulse+r->rlen);
          out[i].length = r->rlen;
          out[i].offset = 0;
        }
    }
  release_prefix(s->idx,sizeof(DAZZ_DB)+(size_t) s->cursor*sizeof(DAZZ_READ),&s->idx_released);
  if (s->hdr)
    release_prefix(s->hdr,header_end,&s->hdr_released);
  return out;
}

void Paf_Free_Names(Paf_Name *names, int n)
{ if (names)
    { for (int i = 0; i < n; ++i)
        free(names[i].name);
      free(names);
    }
}

void Paf_Close_Source(Paf_Source *s)
{ if (!s)
    return;
  munmap(s->idx,s->idx_size);
  if (s->hdr)
    munmap(s->hdr,s->hdr_size);
  if (s->fai)
    fclose(s->fai);
  if (s->stub)
    Free_DB_Stub(s->stub);
  free(s->fai_line);
  free(s->last_name);
  free(s);
}

void Paf_Set_Names(Paf_Name *reads, Paf_Name *reference)
{ query_names = reads;
  ref_names   = reference;
}

void Paf_Open_Output(void)
{ output = fopen(PAF_PATH,"w");
  if (!output)
    die(PAF_PATH);
  setvbuf(output,NULL,_IOFBF,1 << 20);
}

void Paf_Close_Output(void)
{ if (output && fclose(output))
    die("failed to close PAF output");
  output = NULL;
}

static Cigar_Entry *lookup(Cigar_Cache *s, uint64_t key)
{ size_t lo = 0, hi = s->n;

  while (lo < hi)
    { size_t m = lo+(hi-lo)/2;

      if (s->entry[m].key < key)
        lo = m+1;
      else
        hi = m;
    }
  return lo < s->n && s->entry[lo].key == key ? s->entry+lo : NULL;
}

const Cigar *Cigar_Cache_Get(Cigar_Cache *s, uint64_t key)
{ Cigar_Entry *e = lookup(s,key);

  return e ? &e->cigar : NULL;
}

void Cigar_Cache_Put(Cigar_Cache *s, uint64_t key, const Path *p, const Cigar *c)
{ if (s->n && key <= s->entry[s->n-1].key)
    die("unordered traceback cache");
  if (s->n == s->cap)
    { s->cap = s->cap+s->cap/2+32;
      s->entry = checked(realloc(s->entry,s->cap*sizeof(*s->entry)));
    }

  Cigar_Entry *e = s->entry+s->n++;

  *e = (Cigar_Entry) {.key = key, .ab = p->abpos, .bb = p->bbpos};
  Cigar_Append(&e->cigar,c);
  s->bytes += e->cigar.cap*sizeof(uint32_t);

  size_t bytes = s->bytes+s->cap*sizeof(*e);

  if (bytes > s->peak)
    s->peak = bytes;
}

void Cigar_Cache_Clear(Cigar_Cache *s)
{ for (size_t i = 0; i < s->n; ++i)
    Cigar_Free(&s->entry[i].cigar);
  s->n = 0;
  s->bytes = 0;
}

void Cigar_Cache_Free(Cigar_Cache *s)
{ Cigar_Cache_Clear(s);
  free(s->entry);
  *s = (Cigar_Cache) {0};
}

void Cigar_Cache_Fuse(Cigar_Cache *s, uint64_t key, const Path *a,
                       const Path *b, int q, int r)
{ Cigar_Entry *x = lookup(s,(uint64_t) a->trace), *y = lookup(s,(uint64_t) b->trace);

  if (!x || !y)
    return;

  Cigar left = {0}, right = {0};

  if (Cigar_Slice(&left,&x->cigar,q-x->ab,r-x->bb,0) &&
      Cigar_Slice(&right,&y->cigar,q-y->ab,r-y->bb,1))
    { Cigar_Append(&left,&right);
      Cigar_Cache_Put(s,key,a,&left);
    }
  /* Sparse checkpoints can differ from the discovery path at an indel.
     Such fused paths are reconstructed on demand, never silently spliced. */
  Cigar_Free(&left);
  Cigar_Free(&right);
}

double Paf_Thread_Seconds(void)
{ struct timespec t;

  clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);
  return t.tv_sec+t.tv_nsec*1e-9;
}

double Paf_Wall_Seconds(void)
{ struct timespec t;

  clock_gettime(CLOCK_MONOTONIC,&t);
  return t.tv_sec+t.tv_nsec*1e-9;
}

static const Cigar *reconstruct(Paf_Work *w, const Overlap *ovl,
                                const char *query, const char *reverse_query, int qlen,
                                const char *reference, int rlen, int spacing)
{ Path p = ovl->path;
  int  reverse = COMP(ovl->flags);

  /* Reporter already has the reverse-complemented query for reverse hits.
     Traverse sparse tiles backwards and reuse it, avoiding a reference-slice
     copy and any additional sequence buffer even for ultra-long reads. */
  const char *a = reverse ? reverse_query : query;

  if (!a)
    die("reverse query unavailable");
  /* A zero-difference checkpoint interval is already known to be an exact
     match. Reconstruct only the other (normally ~100 bp) intervals. This
     bounds DP scratch space by checkpoint spacing, not read/reference length. */
  uint16 *points = p.trace;
  int     nt = p.tlen, aq = reverse ? qlen-p.aepos : p.abpos;
  int     br = reverse ? rlen-p.bepos : p.bbpos;
  int     qend = reverse ? qlen-p.abpos : p.aepos;
  int     rend = reverse ? rlen-p.bbpos : p.bepos;

  w->cigar.n = 0;
  for (int order = 0; order < nt; order += 2)
    { int tile = reverse ? nt-2-order : order;
      int lo = tile == 0 ? p.abpos : (p.abpos/spacing+tile/2)*spacing;
      int hi = tile+2 == nt ? p.aepos : (p.abpos/spacing+tile/2+1)*spacing;
      int ae = aq+hi-lo;
      int be = br+points[tile+1];

      if (ae < aq || ae > qend || be < br || be > rend)
        die("invalid sparse traceback");
      if (points[tile] == 0)
        { if (ae-aq != be-br)
            die("zero-difference checkpoint has unequal lengths");
          Cigar_Push(&w->cigar,CG_EQ,ae-aq);
          ++w->exact_tiles;
          w->exact_bases += ae-aq;
        }
      else
        { ++w->traced_tiles;
          w->traced_bases += ae-aq;

          Path part = {.trace = points+tile, .tlen = 2, .diffs = points[tile],
                       .abpos = 0, .aepos = ae-aq, .bbpos = 0, .bepos = be-br};
          const char *qa = a+aq, *rb = reference+br;
          Alignment aln = {&part,0,(char *) qa,(char *) rb,ae-aq,be-br};

          if (!w->trace_work)
            w->trace_work = New_Work_Data();
          if (Compute_Trace_PTS(&aln,w->trace_work,spacing,GREEDIEST))
            die("trace reconstruction failed");

          int q = 0, r = 0, *trace = part.trace;

          for (int i = 0; i <= part.tlen; ++i)
            { int qe, re;

              if (i == part.tlen)
                { qe = ae-aq;
                  re = be-br;
                }
              else if (trace[i] < 0)
                { qe = -trace[i]-1;
                  re = r+qe-q;
                }
              else
                { re = trace[i]-1;
                  qe = q+re-r;
                }
              if (qe < q || re < r || qe > ae-aq || re > be-br)
                die("invalid detailed traceback");
              while (q < qe)
                { int start = q, equal = qa[q] == rb[r];

                  do
                    { ++q;
                      ++r;
                    }
                  while (q < qe && (qa[q] == rb[r]) == equal);
                  Cigar_Push(&w->cigar,equal ? CG_EQ : CG_X,q-start);
                }
              if (i < part.tlen)
                { if (trace[i] < 0)
                    { Cigar_Push(&w->cigar,CG_D,1);
                      ++r;
                    }
                  else
                    { Cigar_Push(&w->cigar,CG_I,1);
                      ++q;
                    }
                }
            }
          if (q != ae-aq || r != be-br)
            die("traceback does not consume checkpoint bounds");
        }
      aq = ae;
      br = be;
    }
  if (aq != qend || br != rend)
    die("sparse traceback does not consume alignment bounds");
  return &w->cigar;
}

/* CIGAR has many tiny runs on noisy long reads. Avoid a printf parser for
   every run; the row remains one buffered, thread-safe fwrite operation. */

static char *print_run(char *out, unsigned n, char op)
{ char digits[10];
  int  count = 0;

  do
    { digits[count++] = '0'+n%10;
      n /= 10;
    }
  while (n);
  while (count)
    *out++ = digits[--count];
  *out++ = op;
  return out;
}

void Paf_Emit(Paf_Work *w, const Overlap *ovl, int read_index,
              const char *query, const char *reverse_query, int qlen,
              const char *reference, int rlen,
              int spacing, const Cigar *saved, int best, int chain, int segment, int mapq)
{ double       start = PAF_STATS ? Paf_Thread_Seconds() : 0;
  const Cigar *cg = saved;

  if ((mapq < 0 || mapq > 60) && mapq != 255)
    die("MAPQ must be 0-60 or 255");

  if (!cg)
    { cg = reconstruct(w,ovl,query,reverse_query,qlen,reference,rlen,spacing);
      if (PAF_HISTORY)
        ++w->fallback;
    }

  int64_t q, r, matches, columns, edits;

  Cigar_Count(cg,&q,&r,&matches,&columns,&edits);

  const Path *p = &ovl->path;

  if (q != p->aepos-p->abpos || r != p->bepos-p->bbpos)
    die("CIGAR length mismatch");

  const Paf_Name *qn = query_names+read_index, *rn = ref_names+ovl->bread;
  int            rb = COMP(ovl->flags) ? rlen-p->bepos : p->bbpos;
  size_t         need = strlen(qn->name)+strlen(rn->name)+cg->n*12+512;

  if (need > w->text_cap)
    { w->text_cap = need+4096;
      w->text = checked(realloc(w->text,w->text_cap));
    }

  int n = snprintf(w->text,w->text_cap,
                   "%s\t%" PRId64 "\t%d\t%d\t%c\t%s\t%" PRId64 "\t%d\t%d\t%" PRId64
                   "\t%" PRId64 "\t%d\ttp:A:%c\tNM:i:%" PRId64 "\thi:i:%d\tsi:i:%d\tcg:Z:",
                   qn->name,qn->length,qn->offset+p->abpos,qn->offset+p->aepos,
                   COMP(ovl->flags) ? '-' : '+',rn->name,rn->length,
                   rn->offset+rb,rn->offset+rb+(int) r,
                   matches,columns,mapq,best ? 'P' : 'S',edits,chain,segment);
  static const char ops[] = "MIDNSHP=X";
  char *end = w->text+n;

  for (size_t i = 0; i < cg->n; ++i)
    { uint32_t v = cg->ops[saved && COMP(ovl->flags) ? cg->n-1-i : i];

      end = print_run(end,v >> 4,ops[v & 15]);
    }
  *end++ = '\n';

  size_t length = end-w->text;

  if (fwrite(w->text,1,length,output) != length)
    die("PAF write failed");
  ++w->selected;
  w->output_ops += cg->n;
  if (PAF_STATS)
    w->finish_cpu += Paf_Thread_Seconds()-start;

  size_t bytes = w->text_cap+w->cigar.cap*sizeof(uint32_t)+
                 (w->trace_work ? Alignment_Work_Bytes(w->trace_work) : 0);

  if (bytes > w->peak)
    w->peak = bytes;
}

void Paf_Free_Work(Paf_Work *w)
{ if (w->trace_work)
    Free_Work_Data(w->trace_work);
  Cigar_Free(&w->cigar);
  free(w->text);
}
