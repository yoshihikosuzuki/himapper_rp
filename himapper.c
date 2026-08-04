/*********************************************************************************************\
 *
 *  Find best map of long, HiFi DNA reads to a set of reference contigs:
 *
 *  Author:  Gene Myers
 *  Date  :  Feb 9, 2016
 *
 *********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#include <DB.h>
#include "map.h"
#include "minimaxmer.h"

static char *Usage[] =
  { "[-1vpzCNG] [-k<int(20)>] [-t<int>] [-M<int>] [-T<int(4)>] [-P<dir(/tmp)>]",
    "         [-e<double(.85)] [-s<int(100)>] [-n<double(1.00)>]",
    "         [seed generation long options]",
    "         [-m<track>]+  <reference:dam> <ref.kpos|all> <reads:db> <reads.kpos|all|generate>",
  };

enum
  { SEED_MODE_NONE = 0,
    SEED_MODE_SERIAL = 1,
    SEED_MODE_PIPELINE = 2
  };

int     VERBOSE;   //   Globally visible to map.c
int     PROFILE;
int     SPACING;
double  BEST_TIE;
int     PRIMARY_ONLY;
int     JUNCTION_PEN;
char   *SORT_PATH;

uint64  MEM_LIMIT;
uint64  MEM_PHYSICAL;

typedef struct
  { uint64 reference_db_load_ns;
    uint64 reference_kpos_cache_ns;
    uint64 serial_seed_batch_ns;
    uint64 pipeline_seed_wait_ns;
    uint64 read_db_load_ns;
    uint64 read_kpos_write_ns;
    uint64 read_seed_prepare_ns;
    uint64 filter_workers_ns;
    uint64 match_filter_ns;
    uint64 reporter_ns;
    uint64 lasort_ns;
    uint64 lacat_ns;
    uint64 lamerge_ns;
    uint64 read_blocks;
  } Himapper_Phase_Stats;

static Himapper_Phase_Stats HIMAPPER_PHASE_STATS;

static uint64 himapper_now_ns(void)
{
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC,&ts) != 0)
    return 0;
  return (uint64)ts.tv_sec * 1000000000llu + (uint64)ts.tv_nsec;
}

static void himapper_add_elapsed(uint64 *total, uint64 start)
{
  uint64 end;

  if (start == 0)
    return;
  end = himapper_now_ns();
  if (end >= start)
    *total += end-start;
}

static void print_himapper_phase_stats(FILE *out, int enabled)
{
  if (!enabled || out == NULL)
    return;
#define PRINT_OUTER_PHASE(field)                                        \
  fprintf(out,"PHASE_V1\thimapper\touter.%s\tserial_sum\t%.9f\t%llu\n", \
          #field,(double)HIMAPPER_PHASE_STATS.field##_ns/1e9,           \
          (unsigned long long)HIMAPPER_PHASE_STATS.read_blocks)
  PRINT_OUTER_PHASE(reference_db_load);
  PRINT_OUTER_PHASE(reference_kpos_cache);
  PRINT_OUTER_PHASE(serial_seed_batch);
  PRINT_OUTER_PHASE(pipeline_seed_wait);
  PRINT_OUTER_PHASE(read_db_load);
  PRINT_OUTER_PHASE(read_kpos_write);
  PRINT_OUTER_PHASE(read_seed_prepare);
  PRINT_OUTER_PHASE(filter_workers);
  PRINT_OUTER_PHASE(match_filter);
  PRINT_OUTER_PHASE(reporter);
  PRINT_OUTER_PHASE(lasort);
  PRINT_OUTER_PHASE(lacat);
  PRINT_OUTER_PHASE(lamerge);
#undef PRINT_OUTER_PHASE
}

typedef struct
  { DAZZ_DB view;
    void   *forward;
    void   *reverse;
    int     forward_len;
    int     reverse_len;
  } Ref_Block_Cache;

/*  Adapted from code by David Robert Nadeau (http://NadeauSoftware.com) licensed under
 *     "Creative Commons Attribution 3.0 Unported License"
 *          (http://creativecommons.org/licenses/by/3.0/deed.en_US)
 *
 *   I removed Windows options, reformated, and return int64 instead of size_t
 */

static int64 getMemorySize( )
{
#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))

  // OSX, NetBSD, OpenBSD

  int     mib[2];
  size_t  size = 0;
  size_t  len = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
  mib[1] = HW_MEMSIZE;            // OSX
#elif defined(HW_PHYSMEM64)
  mib[1] = HW_PHYSMEM64;          // NetBSD, OpenBSD
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return ((size_t) size);
  return (0);

#elif defined(_SC_AIX_REALMEM)

  // AIX

  return ((size_t) sysconf( _SC_AIX_REALMEM ) * ((size_t) 1024L));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  // FreeBSD, Linux, OpenBSD, & Solaris

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGESIZE)));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)

  // ? Legacy ?

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGE_SIZE)));

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))

  // DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX

  int          mib[2];
  unsigned int size = 0;
  size_t       len  = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_REALMEM)
  mib[1] = HW_REALMEM;		// FreeBSD
#elif defined(HW_PYSMEM)
  mib[1] = HW_PHYSMEM;		// Others
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return (size_t)size;
  return (0);

#else

  return (0);

#endif
}

typedef struct
  { int *ano;
    int *end;
    int  idx;
    int  out;
  } Event;

// TODO: hl,dlにはない
static void reheap(int s, Event **heap, int hsize)
{ int      c, l, r;
  Event   *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      hr = heap[r];
      if (hr->idx > hl->idx)
        { if (hs->idx > hl->idx)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (hs->idx > hr->idx)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

// TODO: hlにはない
static int64 merge_size(DAZZ_DB *block, int mtop)
{ Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64       nsize;

  { DAZZ_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      DAZZ_TRACK *track;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                nsize += 1;
            }
          else
            { if (level == 0)
                nsize += 1;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
    }

  return (nsize);
}

// TODO: hlにはない
static DAZZ_TRACK *merge_tracks(DAZZ_DB *block, int mtop, int64 nsize)
{ DAZZ_TRACK *ntrack;
  Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64      *anno;
  int        *alen;
  int        *data;

  ntrack = (DAZZ_TRACK *) Malloc(sizeof(DAZZ_TRACK),"Allocating merged track");
  if (ntrack == NULL)
    Clean_Exit(1);
  ntrack->name = Strdup("merge","Allocating merged track");
  ntrack->anno = anno = (int64 *) Malloc(sizeof(int64)*(block->nreads+1),"Allocating merged track");
  ntrack->alen = alen = (int *) Malloc(sizeof(int)*block->nreads,"Allocating merged track");
  ntrack->data = data = (int *) Malloc(sizeof(int)*nsize,"Allocating merged track");
  ntrack->size = sizeof(int);
  ntrack->next = NULL;
  ntrack->loaded = 1;
  if (anno == NULL || alen == NULL ||  data == NULL || ntrack->name == NULL)
    Clean_Exit(1);

  { DAZZ_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      DAZZ_TRACK *track;

      anno[r] = nsize;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                data[nsize++] = p->idx;
            }
          else
            { if (level == 0)
                data[nsize++] = p->idx;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
      alen[r] = (int) (nsize - anno[r]);
    }
  anno[r] = nsize;

  return (ntrack);
}

static int read_DB(DAZZ_DB *block, char *name, char **mask, int *mstat, int mtop, int kmer)
{ int i, isdam, status, kind, stop;

  isdam = Minimaxmer_Open_DB_Serialized(name,block);
  if (isdam < 0)
    Clean_Exit(1);

  for (i = 0; i < mtop; i++)
    { status = Check_Track(block,mask[i],&kind);
      if (status >= 0)
        { if (kind != MASK_TRACK)
            { fprintf(stderr,"%s: %s track is not a mask track.\n",Prog_Name,mask[i]);
              exit (1);
            }
          if (status == 0)
            Open_Track(block,mask[i]);
          mstat[i] = 1;
        }
      else if (status == -1)
        { printf("%s: Warning: %s track not sync'd with db %s, ignored.\n",
                      Prog_Name,mask[i],name);
        }
    }

  Trim_DB(block);

  stop = 0;
  for (i = 0; i < mtop; i++)
    { DAZZ_TRACK *track;
      int64      *anno;
      int         j;

      status = Check_Track(block,mask[i],&kind);
      if (status < 0)
        continue;

      stop += 1;
      track = Open_Track(block,mask[i]);
      Load_All_Track_Data(track);

      anno = (int64 *) (track->anno); 
      for (j = 0; j <= block->nreads; j++)
        anno[j] /= sizeof(int);
    }

  if (stop > 1)
    { int64       nsize;
      DAZZ_TRACK *track;

      nsize = merge_size(block,stop);
      track = merge_tracks(block,stop,nsize);

      while (block->tracks != NULL)
        Close_Track(block,block->tracks);

      block->tracks = track;
    }

  if (block->cutoff < kmer)
    { for (i = 0; i < block->nreads; i++)
        if (block->reads[i].rlen < kmer)
          { fprintf(stderr,"%s: Block %s contains reads < %dbp long !  Run DBsplit -x%d\n",
                           Prog_Name,name,kmer,kmer);
            Clean_Exit(1);
          }
    }

  Load_All_Reads(block,0);

  return (isdam);
}

static void reference_block_view(DAZZ_DB *view, DAZZ_DB *whole,
                                 int part, int first, int last)
{ int i;

  *view = *whole;
  view->part   = part;
  view->ufirst = first;
  view->tfirst = first;
  view->nreads = last-first;
  view->reads  = whole->reads+first;
  view->tracks = NULL;
  view->totlen = 0;
  view->maxlen = 0;

  for (i = first; i < last; i++)
    { int len = whole->reads[i].rlen;
      view->totlen += len;
      if (len > view->maxlen)
        view->maxlen = len;
    }
}

// TODO: hl,dlにはない
static void complement(char *s, int len)
{ char *t;
  int   c;

  t = s + (len-1);
  while (s < t)
    { c = *s;
      *s = (char) (3-*t);
      *t = (char) (3-c);
      s += 1;
      t -= 1;
    }
  if (s == t)
    *s = (char) (3-*s);
}

// TODO: hl,dlにはない
static DAZZ_DB *complement_DB(DAZZ_DB *block, int inplace)
{ static DAZZ_DB _cblock, *cblock = &_cblock;
  int            nreads;
  DAZZ_READ     *reads;
  char          *seq;
  
  nreads = block->nreads;
  reads  = block->reads;
  if (inplace)
    { seq = (char *) block->bases;
      cblock = block;
    }
  else
    { seq  = (char *) Malloc(block->reads[nreads].boff+1,"Allocating dazzler sequence block");
      if (seq == NULL)
        Clean_Exit(1);
      *seq++ = 4;
      memmove(seq,block->bases,block->reads[nreads].boff);
      *cblock = *block;
      cblock->bases  = (void *) seq;
      cblock->tracks = NULL;
    }

  { int   i;
    float x;

    x = cblock->freq[0];
    cblock->freq[0] = cblock->freq[3];
    cblock->freq[3] = x;

    x = cblock->freq[1];
    cblock->freq[1] = cblock->freq[2];
    cblock->freq[2] = x;

    for (i = 0; i < nreads; i++)
      complement(seq+reads[i].boff,reads[i].rlen);
  }

  { DAZZ_TRACK *src, *trg;
    int        *data, *tata;
    int         i, x, rlen;
    int64      *tano, *anno;
    int64       j, k;

    for (src = block->tracks; src != NULL; src = src->next)
      { tano = (int64 *) src->anno;
        tata = (int   *) src->data;

        if (inplace)
          { data = tata;
            anno = tano;
            trg  = src;
          }
        else
          { data = (int *) Malloc(sizeof(int)*tano[nreads],
                                  "Allocating dazzler interval track data");
            anno = (int64 *) Malloc(sizeof(int64)*(nreads+1),
                                    "Allocating dazzler interval track index");
            trg  = (DAZZ_TRACK *) Malloc(sizeof(DAZZ_TRACK),
                                         "Allocating dazzler interval track header");
            if (data == NULL || trg == NULL || anno == NULL)
              Clean_Exit(1);

            trg->name = Strdup(src->name,"Copying track name");
            if (trg->name == NULL)
              Clean_Exit(1);

            trg->size = 4;
            trg->anno = (void *) anno;
            trg->data = (void *) data;
            trg->next = cblock->tracks;
            cblock->tracks = trg;
          }

        for (i = 0; i < nreads; i++)
          { rlen = reads[i].rlen;
            anno[i] = tano[i];
            j = tano[i+1]-1;
            k = tano[i];
            while (k < j)
              { x = tata[j];
                data[j--] = rlen - tata[k];
                data[k++] = rlen - x;
              }
            if (k == j)
              data[k] = rlen - tata[k];
          }
        anno[nreads] = tano[nreads];
      }
  }

  return (cblock);
}

static char *CommandBuffer(char *aname, char *bname, char *spath)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = 2*(strlen(aname) + strlen(bname) + strlen(spath)) + 200;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          Clean_Exit(1);
        }
    }
  return (cat);
}

void Clean_Exit(int val)
{ char *command;

  command = CommandBuffer("","",SORT_PATH);
  sprintf(command,"rm -r %s",SORT_PATH);
  if (system(command) != 0)
    { fprintf(stderr,"%s: Command Failed:\n%*s      %s\n",
                     Prog_Name,(int) strlen(Prog_Name),"",command);
      exit (1);
    }
  exit (val);
}

int main(int argc, char *argv[])
{ DAZZ_DB    _refdb, _ablock, _bblock;
  DAZZ_DB    *refdb = &_refdb, *ablock = &_ablock, *bblock = &_bblock;
  // TODO: ここに[a|b]class
  char       *afile,  *bfile;
  char       *aroot,  *broot, *apwd;
  char       *kpos_path;
  char       *bkpos_path;
  Kpos_Index *ref_kpos = NULL, *read_kpos = NULL;
  Minimaxmer_Context *seed_context = NULL;
  Minimaxmer_Kpos_Writer *seed_writer = NULL;
  Ref_Block_Cache *ref_cache = NULL;
  DAZZ_STUB  *ref_stub = NULL, *read_stub = NULL;
  void       *aindex, *bindex;
  int         alen,    blen;
  int         isdam, nblocks;   //  of reference (a-argument)
  int         mflag;
  char       *command;
  char       *seed_genome_profile = NULL;
  char       *seed_reads_profile = NULL;
  char       *seed_save_path = NULL;
  int         seed_mode = SEED_MODE_NONE;
  int         seed_threads = 0;
  int         seed_window = 100;
  int         seed_spacing = 100;
  int         seed_window_explicit = 0;
  int         seed_spacing_explicit = 0;
  int         phase_timing = 0;
  int         integrated_seed = 0;
  int         memory_limit_explicit = 0;
  int         map_threads;

  Align_Spec *settings;

  int    COVER, NOMAP;
  int    KMER_LEN;
  int    MAX_REPS;
  double AVE_ERROR;
  int    NTHREADS;
  int    MAP_ORDER;
  int    MTOP, *MSTAT;
  char **MASK;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    DIR   *dirp;

    ARG_INIT("himapper")

    // Consume the generated-seed long options before the historical short
    // option parser compacts argv.  Legacy invocations remain unchanged.
    j = 1;
    for (i = 1; i < argc; i++)
      { if (strcmp(argv[i],"--read-seed-profiles") == 0)
          { if (i+2 >= argc)
              { fprintf(stderr,"%s: --read-seed-profiles needs two paths\n",
                               Prog_Name);
                exit (1);
              }
            seed_genome_profile = argv[++i];
            seed_reads_profile = argv[++i];
          }
        else if (strcmp(argv[i],"--seed-mode") == 0)
          { if (++i >= argc)
              { fprintf(stderr,"%s: --seed-mode needs serial or pipeline\n",
                               Prog_Name);
                exit (1);
              }
            if (strcmp(argv[i],"serial") == 0)
              seed_mode = SEED_MODE_SERIAL;
            else if (strcmp(argv[i],"pipeline") == 0)
              seed_mode = SEED_MODE_PIPELINE;
            else
              { fprintf(stderr,"%s: invalid --seed-mode '%s'\n",
                               Prog_Name,argv[i]);
                exit (1);
              }
          }
        else if (strcmp(argv[i],"--seed-threads") == 0)
          { char *end;
            long value;
            if (++i >= argc)
              { fprintf(stderr,"%s: --seed-threads needs an integer\n",
                               Prog_Name);
                exit (1);
              }
            value = strtol(argv[i],&end,10);
            if (*argv[i] == '\0' || *end != '\0' || value < 1 ||
                value > INT32_MAX)
              { fprintf(stderr,"%s: invalid --seed-threads '%s'\n",
                               Prog_Name,argv[i]);
                exit (1);
              }
            seed_threads = (int)value;
          }
        else if (strcmp(argv[i],"--seed-window") == 0 ||
                 strcmp(argv[i],"--seed-spacing") == 0)
          { int is_window = argv[i][7] == 'w';
            char *end;
            long value;
            const char *option = argv[i];
            if (++i >= argc)
              { fprintf(stderr,"%s: %s needs an integer\n",Prog_Name,option);
                exit (1);
              }
            value = strtol(argv[i],&end,10);
            if (*argv[i] == '\0' || *end != '\0' ||
                value < (is_window ? 1 : 0) || value > INT32_MAX)
              { fprintf(stderr,"%s: invalid %s value '%s'\n",
                               Prog_Name,option,argv[i]);
                exit (1);
              }
            if (is_window)
              { seed_window = (int)value;
                seed_window_explicit = 1;
              }
            else
              { seed_spacing = (int)value;
                seed_spacing_explicit = 1;
              }
          }
        else if (strcmp(argv[i],"--save-read-kpos") == 0)
          { if (++i >= argc)
              { fprintf(stderr,"%s: --save-read-kpos needs a path\n",
                               Prog_Name);
                exit (1);
              }
            seed_save_path = argv[i];
          }
        else if (strcmp(argv[i],"--phase-timing") == 0)
          phase_timing = 1;
        else
          argv[j++] = argv[i];
      }
    argc = j;

    KMER_LEN  = 32;
    MAX_REPS  = 0;
    AVE_ERROR = .97;
    SPACING   = 100;       //   Globally visible to map.c
    BEST_TIE  = 1.0;   // -n: near-optimal threshold; report a same-region
                       //     secondary only if its score >= BEST_TIE*best.
                       //     1.0 => only co-optimal (exact-tie) secondaries.
    PRIMARY_ONLY = 0;  // -1: report only the single best alignment per A-read
                       //     region (suppress all same-region secondaries).
    JUNCTION_PEN = 1;  // collinear-gap penalties (chain drift + linker junction) ON by
                       //     default; -G restores the legacy (no-penalty) damapper
                       //     selection.
    NTHREADS  = 4;
    SORT_PATH = "/tmp";

    MEM_PHYSICAL = getMemorySize();
    MEM_LIMIT    = MEM_PHYSICAL;
    if (MEM_PHYSICAL == 0)
      { fprintf(stderr,"\nWarning: Could not get physical memory size\n");
        fflush(stderr);
      }

    MTOP  = 0;        // mask-track feature disabled, but MTOP/MSTAT/MASK are still
    MASK  = NULL;     // passed to read_DB: initialize them so the mask loop is a no-op.
    MSTAT = NULL;
    // MMAX  = 10;
    // MASK  = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    // MSTAT = (int *) Malloc(MMAX*sizeof(int),"Allocating mask status array");
    // if (MASK == NULL || MSTAT == NULL)
    //   exit (1);

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vpzCN1G")
            break;
          case 'e':
            ARG_REAL(AVE_ERROR)
            if (AVE_ERROR < .7 || AVE_ERROR >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",
                               Prog_Name,AVE_ERROR);
                exit (1);
              }
            break;
          case 'k':
            ARG_POSITIVE(KMER_LEN,"K-mer length")
            // if (KMER_LEN > 32)
            //   { fprintf(stderr,"%s: K-mer length must be 32 or less\n",Prog_Name);
            //     exit (1);
            //   }
            if (KMER_LEN <= 16 || KMER_LEN > 48)
              { fprintf(stderr,"%s: K-mer length must be between 17 and 48\n",Prog_Name);
                exit (1);
              }   // TODO: このKMER_LENに対応しているか確認
            break;
            // TODO: hlだと -m は MASK_DIFF_THRES。こちらもそうする？
          // case 'm':
          //   if (MTOP >= MMAX)
          //     { MMAX  = 1.2*MTOP + 10;
          //       MASK  = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
          //       MSTAT = (int *) Realloc(MSTAT,MMAX*sizeof(int),"Reallocating mask status array");
          //       if (MASK == NULL || MSTAT == NULL)
          //         exit (1);
          //     }
          //   MASK[MTOP++] = argv[i]+2;
          //   break;
          case 'n':
            ARG_REAL(BEST_TIE)
            if (BEST_TIE < .7 || BEST_TIE > 1.)
              { fprintf(stderr,"%s: Near optimal threshold must be in [.7,1.] (%g)\n",
                               Prog_Name,BEST_TIE);
                exit (1);
              }
            break;
          case 's':
            ARG_POSITIVE(SPACING,"Trace spacing")
            break;
          case 't':
            ARG_POSITIVE(MAX_REPS,"Tuple supression frequency")
            break;
          case 'M':
            { int limit;

              ARG_NON_NEGATIVE(limit,"Memory allocation (in Gb)")
              MEM_LIMIT = limit * 0x40000000ll;
              memory_limit_explicit = 1;
              break;
            }
          case 'P':
            SORT_PATH = argv[i]+2;
            if ((dirp = opendir(SORT_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];   //  Globally declared in map.h
    // PROFILE   = flags['p'];
    COVER     = flags['C'];
    NOMAP     = flags['N'];
    MAP_ORDER = 1-flags['z'];
    PRIMARY_ONLY = flags['1'];   //  -1: suppress same-region secondaries (cf. minimap2 --secondary=no)
    JUNCTION_PEN = 1 - flags['G'];  //  -G: disable collinear-gap penalties (legacy damapper)

    if (argc != 5)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: k-mer size (must be <= 32).\n");
        // fprintf(stderr,"      -t: Ignore k-mers that occur >= -t times in a block.\n");
        fprintf(stderr,"      -M: Use only -M GB of memory by ignoring most frequent k-mers.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: Look for alignments with -e percent similarity.\n");
        fprintf(stderr,"      -s: Use -s as the trace point spacing for encoding alignments.\n");
        fprintf(stderr,"      -n: Output all matches within this %% of the best\n");
        fprintf(stderr,"      -1: Report only the single best alignment per read region\n");
        fprintf(stderr,"          (suppress same-region secondaries; cf. minimap2 --secondary=no)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -P: Do sorts and merges in directory -P.\n");
        // fprintf(stderr,"      -m: Soft mask the blocks with the specified mask.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -z: sort .las by A,B-read pairs (overlap piles)\n");
        fprintf(stderr,"          off => sort .las by A-read,A-position pairs");
        fprintf(stderr," (default for mapping)\n");
        // fprintf(stderr,"      -p: Output repeat profile track\n");
        fprintf(stderr,"      -C: Output reference vs reads .las.\n");
        fprintf(stderr,"      -N: Do not output reads vs reference .las.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      For <ref.kpos|all> and <reads.kpos|all>:\n");
        fprintf(stderr,"        Specify block-indexed KPOS containers,\n");
        fprintf(stderr,"        or 'all' to use all k-mers as seeds.\n");
        fprintf(stderr,"        Use 'generate' for reads with:\n");
        fprintf(stderr,"          --read-seed-profiles <genome.prof> <reads.prof>\n");
        fprintf(stderr,"          --seed-mode <serial|pipeline>\n");
        fprintf(stderr,"          --seed-threads <N> --seed-window <W>\n");
        fprintf(stderr,"          --seed-spacing <S> [--save-read-kpos <path>]\n");
        fprintf(stderr,"          [--phase-timing]\n");
        exit (1);
      }

    // argv[1]=ref.dam  argv[2]=ref.kpos|all
    // argv[3]=reads.db argv[4]=reads.kpos|all
      kpos_path  = argv[2];
      bkpos_path = argv[4];

    integrated_seed = strcmp(bkpos_path,"generate") == 0;
    if (getenv("HIMAPPER_PHASE_TIMING") != NULL &&
        strcmp(getenv("HIMAPPER_PHASE_TIMING"),"0") != 0)
      phase_timing = 1;
    if (integrated_seed)
      { if (seed_genome_profile == NULL || seed_reads_profile == NULL ||
            seed_mode == SEED_MODE_NONE)
          { fprintf(stderr,
                    "%s: generated reads seeds require --read-seed-profiles "
                    "and --seed-mode\n",Prog_Name);
            exit (1);
          }
        if (KMER_LEN > 32)
          { fprintf(stderr,
                    "%s: generated KPOS supports k-mer lengths up to 32\n",
                    Prog_Name);
            exit (1);
          }
        if (NTHREADS > HIMAPPER_MAX_THREADS)
          { fprintf(stderr,
                    "%s: integrated seed/map supports a total -T budget "
                    "from 1 through %d\n",Prog_Name,HIMAPPER_MAX_THREADS);
            exit (1);
          }
        if (seed_threads == 0)
          { if (seed_mode == SEED_MODE_PIPELINE)
              { fprintf(stderr,
                        "%s: pipeline mode requires explicit --seed-threads\n",
                        Prog_Name);
                exit (1);
              }
            seed_threads = NTHREADS;
          }
        if (seed_mode == SEED_MODE_SERIAL && seed_threads != NTHREADS)
          { fprintf(stderr,
                    "%s: serial mode gives seed and map the full -T budget; "
                    "--seed-threads must equal -T (%d)\n",
                    Prog_Name,NTHREADS);
            exit (1);
          }
        if (seed_threads > NTHREADS)
          { fprintf(stderr,
                    "%s: seed threads (%d) exceed the total -T budget (%d)\n",
                    Prog_Name,seed_threads,NTHREADS);
            exit (1);
          }
        if (seed_mode == SEED_MODE_PIPELINE && seed_threads >= NTHREADS)
          { fprintf(stderr,
                    "%s: pipeline seed threads (%d) must be below the total "
                    "-T budget (%d)\n",Prog_Name,seed_threads,NTHREADS);
            exit (1);
          }
        if (memory_limit_explicit)
          fprintf(stderr,
                  "%s: warning: integrated seed DB/KPOS memory is included "
                  "in hit-cap planning, but -M does not strictly bound "
                  "profile indexes, worker scratch, or concurrent growth\n",
                  Prog_Name);
      }
    else if (seed_genome_profile != NULL || seed_reads_profile != NULL ||
             seed_mode != SEED_MODE_NONE || seed_threads != 0 ||
             seed_window_explicit || seed_spacing_explicit ||
             seed_save_path != NULL)
      { fprintf(stderr,
                "%s: seed generation options require reads argument 'generate'\n",
                Prog_Name);
        exit (1);
      }

    map_threads = NTHREADS;
    if (seed_mode == SEED_MODE_PIPELINE)
      map_threads -= seed_threads;

    if (COVER)
      if (NOMAP)
        mflag = FLAG_DOB;
      else
        mflag = FLAG_DOA | FLAG_DOB;
    else
      if (NOMAP)
        { fprintf(stderr,"%s: Cannot specify N flag without C also\n",Prog_Name);
          exit (1);
        }
      else
        mflag = FLAG_DOA;

    // if (NOMAP && PROFILE)
    //   { fprintf(stderr,"%s: Cannot specify both N and p flags together\n",Prog_Name);
    //     exit (1);
    //   }

    // for (j = 0; j < MTOP; j++)
    //   MSTAT[j] = -2;
  }

  //  For reference argv[1]: get file path parts (apwd/aroot), whether it is a dam or db,
  //    number of blocks, and set alignment specification with its frequencies.

  // TODO: ここからmapper-specific
  { int     p, part, nfiles, tlast, plen;
    char   *cat, *bptr, *fptr;
    char    fname[MAX_NAME], prolog[MAX_NAME];
    FILE   *dbvis, *index;
    DAZZ_DB dbase;

    plen = strlen(argv[1]);
    if (strcmp(argv[1]+(plen-4),".dam") == 0)
      aroot = Root(argv[1],".dam");
    else
      aroot = Root(argv[1],".db");
    apwd = PathTo(argv[1]);   // TODO: hlだとapath

    bptr = rindex(aroot,'.');
    if (bptr != NULL && bptr[1] != '\0' && bptr[1] != '-')
      { part = strtol(bptr+1,&fptr,10);
        if (*fptr == '\0' && part > 0)
          { fprintf(stderr,"%s: first argument '%s' cannot be a block\n",Prog_Name,argv[1]);
            exit (1);
          }
      }

    isdam = 0;
    cat   = Catenate(apwd,"/",aroot,".db");
    if (cat == NULL)
      return (-1);
    if ((dbvis = fopen(cat,"r")) == NULL)
      { cat = Catenate(apwd,"/",aroot,".dam");
        if (cat == NULL)
          return (-1);
        if ((dbvis = fopen(cat,"r")) == NULL)
          { fprintf(stderr,"%s: Could not open database %s\n",Prog_Name,argv[1]);
            exit (1);
          }
        isdam = 1;
      }

    if (fscanf(dbvis,DB_NFILE,&nfiles) != 1)
      { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,aroot);
        exit (1);
      }
    for (p = 0; p < nfiles; p++)
      if (fscanf(dbvis,DB_FDATA,&tlast,fname,prolog) != 3)
        { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,aroot);
          exit (1);
        }
    if (fscanf(dbvis,DB_NBLOCK,&nblocks) != 1)
      { fprintf(stderr,"%s: DB %s has not yet been partitioned, cannot request a block !\n",
                       Prog_Name,aroot);
        exit (1);
      }

    fclose(dbvis);
    ref_stub = Read_DB_Stub(cat,DB_STUB_BLOCKS);
    if (ref_stub == NULL || ref_stub->nblocks != nblocks)
      { fprintf(stderr,"%s: Cannot read reference block boundaries for %s\n",
                       Prog_Name,aroot);
        exit (1);
      }

    if ((index = Fopen(Catenate(apwd,"/.",aroot,".idx"),"r")) == NULL)
      exit (1);
    if (fread(&dbase,sizeof(DAZZ_DB),1,index) != 1)
      { fprintf(stderr,"%s: Index file (.idx) of %s is junk\n",Prog_Name,aroot);
        exit (1);
      }
    fclose(index);

    settings = New_Align_Spec(AVE_ERROR, SPACING, dbase.freq, 1);
  }
  // TODO: ここまでmapper-specific

  Set_Map_Phase_Timing(phase_timing);
  if ((integrated_seed ?
       Set_Filter_Params_Exact(KMER_LEN,MAX_REPS,map_threads) :
       Set_Filter_Params(KMER_LEN,MAX_REPS,map_threads)))
    { fprintf(stderr,"Illegal combination of filter parameters\n");
      exit (1);
    }
  if (integrated_seed)
    { fprintf(stderr,
              "Integrated read seeds: mode=%s total_threads=%d "
              "seed_threads=%d map_thread_budget=%d map_threads=%d "
              "window=%d spacing=%d "
              "save_kpos=%s phase_timing=%d\n",
              seed_mode == SEED_MODE_SERIAL ? "serial" : "pipeline",
              NTHREADS,seed_threads,map_threads,Filter_Thread_Count(),
              seed_window,seed_spacing,
              seed_save_path == NULL ? "off" : seed_save_path,phase_timing);
    }

  // Create directory in SORT_PATH for file operations

  { char *newpath;

    newpath = (char *) Malloc(strlen(SORT_PATH)+30,"Allocating sort path");
    if (newpath == NULL)
      exit (1);
    sprintf(newpath,"%s/himapper.%d",SORT_PATH,getpid());
    if (mkdir(newpath,S_IRWXU) !=  0)
      { fprintf(stderr,"%s: Could not create directory %s\n",Prog_Name,newpath);
        exit (1);
      }
    SORT_PATH = newpath;
  }

  //  Load the whole reference once.  Reporter needs all reference bases, while
  //  Match_Filter receives lightweight block views and cached block-local KPOS.

  { DAZZ_DB probe;
    char *bpwd, *base;
    char *stub_path;
    char *bsuffix;
    int bdam, b, k, started, seed_i;
    int read_blocks;
    int seed_batch_first = 0;
    int seed_batch_count = 0;
    Minimaxmer_Block *seed_batch = NULL;
    uint64 ref_bytes, cache_bytes;
    uint64 phase_start;

    phase_start = phase_timing ? himapper_now_ns() : 0;
    read_DB(refdb,argv[1],MASK,MSTAT,0,KMER_LEN);
    if (phase_timing)
      himapper_add_elapsed(&HIMAPPER_PHASE_STATS.reference_db_load_ns,
                           phase_start);
    ref_bytes = sizeof_DB(refdb);

    ref_cache = (Ref_Block_Cache *)
        calloc((size_t)nblocks,sizeof(Ref_Block_Cache));
    if (ref_cache == NULL)
      Clean_Exit(1);

    if (strcmp(kpos_path,"all") != 0)
      { phase_start = phase_timing ? himapper_now_ns() : 0;
        ref_kpos = Open_Kpos_Index(kpos_path,nblocks);
        for (k = 0; k < nblocks; k++)
          { Ref_Block_Cache *cache = ref_cache+k;
            reference_block_view(&cache->view,refdb,k+1,
                                 ref_stub->tblocks[k],
                                 ref_stub->tblocks[k+1]);
            cache->forward =
                Load_Kpos_Block(ref_kpos,k+1,0,&cache->view,
                                &cache->forward_len);
            cache->reverse =
                Load_Kpos_Block(ref_kpos,k+1,1,&cache->view,
                                &cache->reverse_len);
          }
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.reference_kpos_cache_ns,
                               phase_start);
      }

    bdam = Open_DB(argv[3],&probe);
    if (bdam < 0)
      Clean_Exit(1);
    if (probe.part != 0)
      { fprintf(stderr,"%s: reads argument '%s' must be a whole DB\n",
                       Prog_Name,argv[3]);
        Close_DB(&probe);
        Clean_Exit(1);
      }
    Close_DB(&probe);

    bsuffix = bdam ? ".dam" : ".db";
    base = Root(argv[3],bsuffix);
    bpwd = PathTo(argv[3]);
    stub_path =
        Strdup(Catenate(bpwd,"/",base,bsuffix),"Allocating reads stub path");
    read_stub = Read_DB_Stub(stub_path,DB_STUB_BLOCKS);
    free(stub_path);
    if (read_stub == NULL || read_stub->nblocks < 1)
      { fprintf(stderr,"%s: Cannot read block boundaries for %s\n",
                       Prog_Name,argv[3]);
        Clean_Exit(1);
      }
    read_blocks = read_stub->nblocks;
    if (integrated_seed)
      { seed_context = Minimaxmer_Open(argv[3],seed_genome_profile,
                                       seed_reads_profile,KMER_LEN,
                                       seed_window,seed_spacing,phase_timing);
        if (seed_context == NULL ||
            Minimaxmer_Block_Count(seed_context) != read_blocks ||
            Minimaxmer_Set_Retain_Loaded_DB(seed_context,1) != 0)
          { fprintf(stderr,
                    "%s: cannot initialize integrated read seed generation\n",
                    Prog_Name);
            Clean_Exit(1);
          }
        if (seed_save_path != NULL)
          { seed_writer =
                Minimaxmer_Open_Kpos_Writer(seed_save_path,KMER_LEN);
            if (seed_writer == NULL)
              { fprintf(stderr,"%s: cannot open KPOS tee %s\n",
                               Prog_Name,seed_save_path);
                Clean_Exit(1);
              }
          }
        if (seed_mode == SEED_MODE_PIPELINE &&
            Minimaxmer_Start_Stream(seed_context,seed_threads) != 0)
          { fprintf(stderr,"%s: cannot start read seed pipeline\n",Prog_Name);
            Clean_Exit(1);
          }
      }
    else if (strcmp(bkpos_path,"all") != 0)
      read_kpos = Open_Kpos_Index(bkpos_path,read_blocks);

    cache_bytes = 0;
    for (k = 0; k < nblocks; k++)
      { cache_bytes +=
            (uint64)(ref_cache[k].forward_len+2) * 16;
        cache_bytes +=
            (uint64)(ref_cache[k].reverse_len+2) * 16;
      }

#define SYSTEM_CHECK(command)                                           \
 if (VERBOSE)                                                           \
   printf("%s\n",command);                                              \
 if (system(command) != 0)                                              \
   { fprintf(stderr,"\n%s: Command Failed:\n%*s      %s\n",             \
                    Prog_Name,(int) strlen(Prog_Name),"",command);      \
     Clean_Exit(1);                                                     \
   }

    for (b = 1; b <= read_blocks; b++)
      { Minimaxmer_Block stream_block = {0};
        Minimaxmer_Block *generated_block = NULL;
        uint64 seed_resident_bytes = 0;

        if (seed_mode == SEED_MODE_SERIAL)
          { if (b < seed_batch_first ||
                b >= seed_batch_first+seed_batch_count)
              { int remaining = read_blocks-b+1;
                int count = remaining < seed_threads ? remaining : seed_threads;

                free(seed_batch);
                seed_batch = (Minimaxmer_Block *)
                    calloc((size_t)count,sizeof(*seed_batch));
                phase_start = phase_timing ? himapper_now_ns() : 0;
                if (seed_batch == NULL ||
                    Minimaxmer_Generate_Batch(seed_context,b,count,
                                              seed_threads,seed_batch) != count)
                  { fprintf(stderr,
                            "%s: failed to generate serial seed batch at block %d\n",
                            Prog_Name,b);
                    Clean_Exit(1);
                  }
                if (phase_timing)
                  himapper_add_elapsed(
                      &HIMAPPER_PHASE_STATS.serial_seed_batch_ns,
                      phase_start);
                seed_batch_first = b;
                seed_batch_count = count;
              }
            generated_block = seed_batch+(b-seed_batch_first);
            for (seed_i = 0; seed_i < seed_batch_count; seed_i++)
              seed_resident_bytes += seed_batch[seed_i].resident_bytes;
          }

        if (seed_mode == SEED_MODE_PIPELINE)
          { phase_start = phase_timing ? himapper_now_ns() : 0;
            if (Minimaxmer_Next_Block(seed_context,&stream_block) != 1 ||
                stream_block.block != b)
              { fprintf(stderr,
                        "%s: failed to receive generated seed block %d\n",
                        Prog_Name,b);
                Clean_Exit(1);
              }
            if (phase_timing)
              himapper_add_elapsed(
                  &HIMAPPER_PHASE_STATS.pipeline_seed_wait_ns,phase_start);
            generated_block = &stream_block;
            seed_resident_bytes = Minimaxmer_Resident_Bytes(seed_context);
          }

        if (read_blocks == 1)
          bfile = Strdup(argv[3],"Allocating reads DB name");
        else
          bfile = Strdup(Catenate(bpwd,"/",base,
                                  Numbered_Suffix(".",b,bsuffix)),
                         "Allocating reads block name");

        if (integrated_seed)
          { bblock = Minimaxmer_Block_DB(generated_block);
            if (bblock == NULL)
              { fprintf(stderr,
                        "%s: generated block %d has no loaded DB\n",
                        Prog_Name,b);
                Clean_Exit(1);
              }
            broot = Root(bfile,bsuffix);
          }
        else
          { bblock = &_bblock;
            phase_start = phase_timing ? himapper_now_ns() : 0;
            if (read_DB(bblock,bfile,MASK,MSTAT,MTOP,KMER_LEN))
              broot = Root(bfile,".dam");
            else
              broot = Root(bfile,".db");
            if (phase_timing)
              himapper_add_elapsed(&HIMAPPER_PHASE_STATS.read_db_load_ns,
                                   phase_start);
          }
        free(bfile);

        phase_start = phase_timing ? himapper_now_ns() : 0;
        Start_Filter_Workers();
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.filter_workers_ns,
                               phase_start);

        if (VERBOSE)
          printf("\nBuilding index for %s\n",broot);
        if (integrated_seed)
          { if (seed_writer != NULL)
              { phase_start = phase_timing ? himapper_now_ns() : 0;
                if (Minimaxmer_Write_Kpos_Block(seed_writer,generated_block,
                                                seed_context) != 0)
                  { fprintf(stderr,
                            "%s: failed to write generated seed block %d\n",
                            Prog_Name,b);
                    Clean_Exit(1);
                  }
                if (phase_timing)
                  himapper_add_elapsed(
                      &HIMAPPER_PHASE_STATS.read_kpos_write_ns,phase_start);
              }
            phase_start = phase_timing ? himapper_now_ns() : 0;
            bindex = Prepare_Kpos_Block(generated_block->entry,
                                        generated_block->len,0,bblock,&blen);
          }
        else
          { phase_start = phase_timing ? himapper_now_ns() : 0;
            if (strcmp(bkpos_path,"all") == 0)
              bindex = Sort_Kmers(bblock,&blen);
            else
              bindex = Load_Kpos_Block(read_kpos,b,0,bblock,&blen);
          }
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.read_seed_prepare_ns,
                               phase_start);

        Begin_Mapping_Block(bblock);
        started = 0;
        phase_start = phase_timing ? himapper_now_ns() : 0;
        for (k = 1; k <= nblocks; k++)
          { uint64 active_seed_source_bytes =
                sizeof_DB(bblock)+
                (uint64)(blen+2)*sizeof(Minimaxmer_Kpos);
            uint64 seed_extra_bytes;

            if (seed_mode == SEED_MODE_PIPELINE)
              seed_resident_bytes =
                  Minimaxmer_Resident_Bytes(seed_context);
            seed_extra_bytes =
                seed_resident_bytes > active_seed_source_bytes ?
                seed_resident_bytes-active_seed_source_bytes : 0;

            if (strcmp(kpos_path,"all") == 0)
              { if (isdam)
                  afile = Strdup(Catenate(apwd,"/",aroot,
                                         Numbered_Suffix(".",k,".dam")),
                                 "Allocating reference block name");
                else
                  afile = Strdup(Catenate(apwd,"/",aroot,
                                         Numbered_Suffix(".",k,".db")),
                                 "Allocating reference block name");
                read_DB(ablock,afile,MASK,MSTAT,MTOP,KMER_LEN);
                free(afile);

                if (VERBOSE)
                  printf("\nBuilding index for %s.%d\n",aroot,k);
                aindex = Sort_Kmers(ablock,&alen);
                Set_Filter_Extra_Memory(ref_bytes+seed_extra_bytes);
                if (alen > 0 && blen > 0)
                  { if (VERBOSE)
                      printf("\nComparing %s to %s.%d\n",broot,aroot,k);
                    Match_Filter(bblock,ablock,bindex,blen,aindex,alen,
                                 0,!started);
                    started = 1;
                  }
                free(aindex);

                ablock = complement_DB(ablock,1);
                if (VERBOSE)
                  printf("\nBuilding index for c(%s.%d)\n",aroot,k);
                aindex = Sort_Kmers(ablock,&alen);
                if (alen > 0 && blen > 0)
                  { if (VERBOSE)
                      printf("\nComparing %s to c(%s.%d)\n",broot,aroot,k);
                    Match_Filter(bblock,ablock,bindex,blen,aindex,alen,
                                 1,!started);
                    started = 1;
                  }
                free(aindex);
                Close_DB(ablock);
              }
            else
              { Ref_Block_Cache *cache = ref_cache+(k-1);
                uint64 view_bytes = sizeof_DB(&cache->view);
                uint64 active_bytes;
                uint64 extra =
                    ref_bytes > view_bytes ? ref_bytes-view_bytes : 0;

                alen = cache->forward_len;
                aindex = cache->forward;
                active_bytes = (uint64)(alen+2) * 16;
                if (cache_bytes > active_bytes)
                  extra += cache_bytes-active_bytes;
                extra += seed_extra_bytes;
                Set_Filter_Extra_Memory(extra);
                if (alen > 0 && blen > 0)
                  { if (VERBOSE)
                      printf("\nComparing %s to %s.%d\n",broot,aroot,k);
                    Match_Filter(bblock,&cache->view,bindex,blen,aindex,alen,
                                 0,!started);
                    started = 1;
                  }

                alen = cache->reverse_len;
                aindex = cache->reverse;
                active_bytes = (uint64)(alen+2) * 16;
                extra = ref_bytes > view_bytes ? ref_bytes-view_bytes : 0;
                if (cache_bytes > active_bytes)
                  extra += cache_bytes-active_bytes;
                extra += seed_extra_bytes;
                Set_Filter_Extra_Memory(extra);
                if (alen > 0 && blen > 0)
                  { if (VERBOSE)
                      printf("\nComparing %s to c(%s.%d)\n",broot,aroot,k);
                    Match_Filter(bblock,&cache->view,bindex,blen,aindex,alen,
                                 1,!started);
                    started = 1;
                  }
              }
          }
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.match_filter_ns,
                               phase_start);

        if (VERBOSE)
          printf("\nFinding best matches for block %s\n",broot);
        phase_start = phase_timing ? himapper_now_ns() : 0;
        Reporter(broot,bblock,aroot,refdb,settings,mflag);
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.reporter_ns,
                               phase_start);
        phase_start = phase_timing ? himapper_now_ns() : 0;
        Stop_Filter_Workers();
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.filter_workers_ns,
                               phase_start);
        HIMAPPER_PHASE_STATS.read_blocks += 1;

        if (integrated_seed)
          { if (seed_mode == SEED_MODE_PIPELINE)
              Minimaxmer_Release_Block(seed_context,&stream_block);
            else
              Minimaxmer_Free_Block(generated_block);
          }
        else
          { free(bindex);
            Close_DB(bblock);
          }

        command = CommandBuffer(aroot,broot,SORT_PATH);

        if ((mflag & FLAG_DOA) != 0)
          { phase_start = phase_timing ? himapper_now_ns() : 0;
            sprintf(command,"LAsort %s %s %s/%s.%s.M%c.las",VERBOSE?"-v":"",
                            MAP_ORDER?"-a":"",SORT_PATH,broot,aroot,BLOCK_SYMBOL);
            SYSTEM_CHECK(command)
            if (phase_timing)
              himapper_add_elapsed(&HIMAPPER_PHASE_STATS.lasort_ns,
                                   phase_start);

            phase_start = phase_timing ? himapper_now_ns() : 0;
            sprintf(command,"LAcat %s %s/%s.%s.M%c.S >%s.%s.las",VERBOSE?"-v":"",
                            SORT_PATH,broot,aroot,BLOCK_SYMBOL,broot,aroot);
            SYSTEM_CHECK(command)
            if (phase_timing)
              himapper_add_elapsed(&HIMAPPER_PHASE_STATS.lacat_ns,
                                   phase_start);
          }

        if ((mflag & FLAG_DOB) != 0)
          { phase_start = phase_timing ? himapper_now_ns() : 0;
            sprintf(command,"LAsort %s %s %s/%s.%s.R%c.las",VERBOSE?"-v":"",
                            MAP_ORDER?"-a":"",SORT_PATH,aroot,broot,BLOCK_SYMBOL);
            SYSTEM_CHECK(command)
            if (phase_timing)
              himapper_add_elapsed(&HIMAPPER_PHASE_STATS.lasort_ns,
                                   phase_start);

            phase_start = phase_timing ? himapper_now_ns() : 0;
            sprintf(command,"LAmerge %s %s %s.%s %s/%s.%s.R%c.S.las",VERBOSE?"-v":"",
                            MAP_ORDER?"-a":"",aroot,broot,SORT_PATH,aroot,broot,BLOCK_SYMBOL);
            SYSTEM_CHECK(command)
            if (phase_timing)
              himapper_add_elapsed(&HIMAPPER_PHASE_STATS.lamerge_ns,
                                   phase_start);
          }

        free(broot);
      }

    Set_Filter_Extra_Memory(0);
    if (seed_mode == SEED_MODE_PIPELINE &&
        Minimaxmer_Stop_Stream(seed_context) != 0)
      { fprintf(stderr,"%s: failed to stop read seed pipeline\n",Prog_Name);
        Clean_Exit(1);
      }
    if (seed_writer != NULL)
      { phase_start = phase_timing ? himapper_now_ns() : 0;
        if (Minimaxmer_Close_Kpos_Writer(seed_writer,seed_context) != 0)
          { fprintf(stderr,"%s: failed to finalize optional read KPOS\n",
                           Prog_Name);
            Clean_Exit(1);
          }
        if (phase_timing)
          himapper_add_elapsed(&HIMAPPER_PHASE_STATS.read_kpos_write_ns,
                               phase_start);
        seed_writer = NULL;
      }
    if (seed_context != NULL)
      { Minimaxmer_Print_Phase_Stats(seed_context,stderr,
                                     "integrated_read_seed");
        Minimaxmer_Close(seed_context);
        seed_context = NULL;
      }
    Print_Map_Phase_Stats(stderr);
    print_himapper_phase_stats(stderr,phase_timing);
    free(seed_batch);
    Close_Kpos_Index(read_kpos);
    Close_Kpos_Index(ref_kpos);
    for (k = 0; k < nblocks; k++)
      { free(ref_cache[k].forward);
        free(ref_cache[k].reverse);
      }
    free(ref_cache);
    Free_DB_Stub(read_stub);
    Free_DB_Stub(ref_stub);
    Close_DB(refdb);
    free(base);
    free(bpwd);
  }

  Clean_Exit(0);
}
