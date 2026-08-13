// Development-only modulo-hash seed-selector variant.
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <DB.h>
#include "nthash.h"
#include "libfastk.h"
#include "copynum.h"

int K = 32;          // k-mer length; maximum 32   [-k]
int EXCLUDE_PERIODIC_BP = -1; // -1 disables long p=1..3 tract exclusion

static const uint64_t MODIMIZER_MOD = UINT64_C(101);
static const uint64_t MODIMIZER_THRESHOLD = UINT64_C(1);

// ASCII base -> 2-bit value (A=0, C=1, G=2, T=3), matching DAZZ_DB.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

// Bitmask covering the low 2*K bits of a k-mer's 2-bit code.
static uint64_t KMASK;

typedef struct
  { uint64_t code;
    int32_t  rpos;
    uint32_t rid;
  } KposEntry;

typedef char KposEntry_must_be_16_bytes[(sizeof(KposEntry) == 16) ? 1 : -1];

typedef struct
  { KposEntry *entry;
    size_t     len;
    size_t     cap;
  } KposVector;

typedef struct
  { uint64_t periodic_tracts[3];
    uint64_t periodic_tract_bases[3];
    uint64_t periodic_reads;
    uint64_t periodic_seeds_removed;
  } Periodic_Stats;

typedef struct
  { KposEntry *entry;
    size_t     len;
    int        ready;
    Periodic_Stats periodic_stats;
  } BlockResult;

typedef struct
  { DAZZ_DB       *db;
    DAZZ_STUB     *stub;
    Profile_Index *genome;
    BlockResult   *result;
    int            next_block;
    int            end_block;
    pthread_mutex_t lock;
    pthread_cond_t  ready;
    pthread_cond_t  consumed;
  } WorkQueue;

static void push_kpos(KposVector *vector, uint64_t code, int32_t rpos,
                      uint32_t rid)
{
  if (vector->len == vector->cap)
    { size_t next = vector->cap == 0 ? 65536 : 2 * vector->cap;
      KposEntry *entry =
          (KposEntry *)realloc(vector->entry, next * sizeof(KposEntry));
      if (entry == NULL)
        { perror("realloc KPOS block");
          exit(1);
        }
      vector->entry = entry;
      vector->cap = next;
    }

  vector->entry[vector->len].code = code;
  vector->entry[vector->len].rpos = rpos;
  vector->entry[vector->len].rid = rid;
  vector->len += 1;
}

static void periodic_stats_add(Periodic_Stats *dst,
                               const Periodic_Stats *src)
{
  int period;

  dst->periodic_reads += src->periodic_reads;
  dst->periodic_seeds_removed += src->periodic_seeds_removed;
  for (period = 0; period < 3; period++)
    { dst->periodic_tracts[period] += src->periodic_tracts[period];
      dst->periodic_tract_bases[period] +=
          src->periodic_tract_bases[period];
    }
}

static int cmp_kpos(const void *a, const void *b)
{
  const KposEntry *ka = (const KposEntry *)a;
  const KposEntry *kb = (const KposEntry *)b;

  if (ka->code < kb->code) return -1;
  if (ka->code > kb->code) return 1;
  return 0;
}

static int write_kpos_header_for_k(FILE *out, int kmer)
{
  char magic[4] = {'K', 'P', 'O', 'S'};
  uint8_t ver = 2, flags = 0, idx_bytes = sizeof(uint64_t), reserved = 0;
  int32_t Kval = kmer;
  int32_t entry_size = (int32_t)sizeof(KposEntry);
  uint64_t n_entries = 0;

  if (fwrite(magic, 4, 1, out) != 1 ||
      fwrite(&ver, 1, 1, out) != 1 ||
      fwrite(&flags, 1, 1, out) != 1 ||
      fwrite(&idx_bytes, 1, 1, out) != 1 ||
      fwrite(&reserved, 1, 1, out) != 1 ||
      fwrite(&Kval, sizeof(int32_t), 1, out) != 1 ||
      fwrite(&entry_size, sizeof(int32_t), 1, out) != 1 ||
      fwrite(&n_entries, sizeof(uint64_t), 1, out) != 1)
    return -1;

  return 0;
}

static int finish_kpos(FILE *out, uint64_t n_entries)
{
  if (fseek(out, 16, SEEK_SET) != 0 ||
      fwrite(&n_entries, sizeof(uint64_t), 1, out) != 1)
    return -1;

  return 0;
}

static inline void seed_bit_set(uint64_t *bits, size_t pos)
{
  bits[pos >> 6] |= UINT64_C(1) << (pos & 63);
}

// Clear the inclusive bit range [lo, hi] without shifting by 64 bits.
static inline __attribute__((always_inline))
uint64_t seed_bit_clear_range(uint64_t *bits, size_t lo, size_t hi)
{
  size_t first = lo >> 6;
  size_t last = hi >> 6;
  size_t word;
  unsigned int lbit = (unsigned int)(lo & 63);
  unsigned int hbit = (unsigned int)(hi & 63);
  uint64_t left_mask = UINT64_MAX << lbit;
  uint64_t right_mask = UINT64_MAX >> (63U - hbit);
  uint64_t removed = 0;

  if (first == last)
    { uint64_t mask = left_mask & right_mask;
      removed = (uint64_t)__builtin_popcountll(bits[first] & mask);
      bits[first] &= ~mask;
      return removed;
    }

  removed += (uint64_t)__builtin_popcountll(bits[first] & left_mask);
  bits[first] &= ~left_mask;
  for (word = first+1; word < last; word++)
    { removed += (uint64_t)__builtin_popcountll(bits[word]);
      bits[word] = 0;
    }
  removed += (uint64_t)__builtin_popcountll(bits[last] & right_mask);
  bits[last] &= ~right_mask;
  return removed;
}

static inline int periodic_base_code(unsigned char base)
{
  switch (base)
    { case 'A': case 'a': return 0;
      case 'C': case 'c': return 1;
      case 'G': case 'g': return 2;
      case 'T': case 't': return 3;
      default: return -1;
    }
}

// Remove selected K-mers that overlap a long exact period-1, -2, or -3
// tract.  The selected seeds that remain are an exact subset of the baseline
// modimizer result.
static void exclude_long_periodic_seeds(
    const char *bases, int nbase, int kmer, int plen,
    uint64_t *selected, Periodic_Stats *stats)
{
  uint64_t removed = 0;
  int found_tract = 0;
  int period;

  if (EXCLUDE_PERIODIC_BP <= 0 || nbase <= 0 || plen <= 0)
    return;

  for (period = 1; period <= 3; period++)
    { int i = period;

      while (i < nbase)
        { int left_code = periodic_base_code(
              (unsigned char)bases[i-period]);
          int right_code = periodic_base_code((unsigned char)bases[i]);
          int first_match;
          int tract_start;
          int tract_end;
          int raw_span;
          int complete_bp;
          int seed_lo;
          int seed_hi;

          if (left_code < 0 || right_code < 0 || left_code != right_code)
            { i++;
              continue;
            }
          first_match = i;
          do
            { i++;
              if (i >= nbase)
                break;
              left_code = periodic_base_code(
                  (unsigned char)bases[i-period]);
              right_code = periodic_base_code((unsigned char)bases[i]);
            }
          while (left_code >= 0 && right_code >= 0 &&
                 left_code == right_code);

          tract_start = first_match-period;
          tract_end = i;
          raw_span = tract_end-tract_start;
          complete_bp = (raw_span/period)*period;
          if (raw_span < 2*period ||
              complete_bp <= EXCLUDE_PERIODIC_BP)
            continue;
          if (period == 2 &&
              periodic_base_code((unsigned char)bases[tract_start]) ==
              periodic_base_code((unsigned char)bases[tract_start+1]))
            continue;
          if (period == 3)
            { int c0 = periodic_base_code(
                  (unsigned char)bases[tract_start]);
              int c1 = periodic_base_code(
                  (unsigned char)bases[tract_start+1]);
              int c2 = periodic_base_code(
                  (unsigned char)bases[tract_start+2]);
              if (c0 == c1 && c1 == c2)
                continue;
            }

          seed_lo = tract_start-kmer+1;
          if (seed_lo < 0)
            seed_lo = 0;
          seed_hi = tract_end-1;
          if (seed_hi >= plen)
            seed_hi = plen-1;
          if (seed_lo <= seed_hi)
            removed += seed_bit_clear_range(
                selected,(size_t)seed_lo,(size_t)seed_hi);
          found_tract = 1;
          if (stats != NULL)
            { stats->periodic_tracts[period-1] += 1;
              stats->periodic_tract_bases[period-1] +=
                  (uint64_t)raw_span;
            }
        }
    }

  if (stats != NULL)
    { if (found_tract)
        stats->periodic_reads += 1;
      stats->periodic_seeds_removed += removed;
    }
}

// Select every k-mer whose full canonical ntHash satisfies hash % 101 < 1.
static int find_seeds(const char *seq, int plen, int kmer,
                      const uint16 *prof_genome, int read_id,
                      KposVector *out, Periodic_Stats *periodic_stats)
{
  size_t seed_nwords;
  uint64_t *selected;
  uint64_t fhVal = 0;
  uint64_t rhVal = 0;
  uint64_t hash;
  int written = 0;
  int i;

  if (plen <= 0)
    return 0;
  seed_nwords = ((size_t)plen + 63) >> 6;
  selected = (uint64_t *)calloc(seed_nwords, sizeof(uint64_t));
  if (selected == NULL)
    { perror("calloc seed bitset");
      exit(1);
    }

  hash = NTC64_b(seq-kmer+1,kmer,&fhVal,&rhVal);
  if (hash % MODIMIZER_MOD < MODIMIZER_THRESHOLD)
    seed_bit_set(selected,0);
  for (i = 1; i < plen; i++)
    { hash = NTC64_c((unsigned char)seq[i-kmer],
                     (unsigned char)seq[i],kmer,&fhVal,&rhVal);
      if (hash % MODIMIZER_MOD < MODIMIZER_THRESHOLD)
        seed_bit_set(selected,(size_t)i);
    }

  if (EXCLUDE_PERIODIC_BP > 0)
    exclude_long_periodic_seeds(seq-kmer+1,plen+kmer-1,kmer,plen,
                                selected,periodic_stats);

  for (size_t word_index = 0; word_index < seed_nwords; word_index++)
    { uint64_t word = selected[word_index];

      while (word != 0)
        { unsigned int bit = (unsigned int)__builtin_ctzll(word);
          size_t pos = (word_index << 6) + bit;
          const char *kp = seq-kmer+1+pos;
          uint64_t code = 0;
          int x;

          for (x = 0; x < kmer; x++)
            code = (code << 2) | base2bit[(unsigned char)kp[x]];
          code &= KMASK;
          push_kpos(out,code,(int32_t)pos,
                    CN_PACK(read_id,cn_tier_of(prof_genome[pos])));
          written++;
          word &= word-1;
        }
    }

  free(selected);
  return written;
}

static void grow_profile(int need, int *capacity, uint16 **prof_genome)
{
  int next = (int)(1.2 * need) + 1000;
  uint16 *genome =
      (uint16 *)realloc(*prof_genome,(size_t)next*sizeof(uint16));

  if (genome == NULL)
    { perror("realloc profile");
      exit(1);
    }
  *prof_genome = genome;
  *capacity = next;
}

static FILE *open_worker_bases(const DAZZ_DB *db)
{
  size_t len = strlen(db->path);
  char *path = (char *)malloc(len+5);
  FILE *bases;

  if (path == NULL)
    { perror("malloc .bps path");
      exit(1);
    }
  snprintf(path,len+5,"%s.bps",db->path);
  bases = fopen(path,"rb");
  if (bases == NULL)
    { perror(path);
      free(path);
      exit(1);
    }
  free(path);
  return bases;
}

static void *seed_worker(void *arg)
{
  WorkQueue *queue = (WorkQueue *)arg;
  DAZZ_DB db = *queue->db;
  Profile_Reader *Pg_reader;
  char *seq;
  uint16 *prof_genome;
  int pmax = 20000;
  KposVector vector = {0};

  db.bases = open_worker_bases(queue->db);
  Pg_reader = Open_Profile_Reader(queue->genome);
  if (Pg_reader == NULL)
    { fprintf(stderr,"modimizer: cannot open sequential profile reader\n");
      exit(1);
    }
  seq = New_Read_Buffer(&db);
  prof_genome = (uint16 *)malloc((size_t)pmax*sizeof(uint16));
  if (prof_genome == NULL)
    { perror("malloc profile");
      exit(1);
    }

  while (1)
    { int block;
      int first, last, id;
      BlockResult *result;
      Periodic_Stats block_stats = {0};

      pthread_mutex_lock(&queue->lock);
      block = queue->next_block++;
      pthread_mutex_unlock(&queue->lock);
      if (block >= queue->end_block)
        break;

      first = queue->stub->ublocks[block];
      last = queue->stub->ublocks[block+1];
      vector.len = 0;
      for (id = first; id < last; id++)
        { DAZZ_READ *r = db.reads+id;
          int plen = r->rlen-K+1;
          int64 fetched;

          Load_Read(&db,id,seq,2);
          fetched = Fetch_Profile_Sequential(
              Pg_reader,(int64)id,pmax,prof_genome);
          if (fetched > pmax)
            { grow_profile((int)fetched,&pmax,&prof_genome);
              Fetch_Profile_Sequential(
                  Pg_reader,(int64)id,pmax,prof_genome);
            }
          if (plen > 0)
            find_seeds(seq+K-1,plen,K,prof_genome,id-first,
                       &vector,&block_stats);
        }

      if (vector.len > 1)
        qsort(vector.entry,vector.len,sizeof(KposEntry),cmp_kpos);

      result = queue->result+block;
      pthread_mutex_lock(&queue->lock);
      result->entry = vector.entry;
      result->len = vector.len;
      result->periodic_stats = block_stats;
      result->ready = 1;
      pthread_cond_broadcast(&queue->ready);
      while (result->ready)
        pthread_cond_wait(&queue->consumed,&queue->lock);
      pthread_mutex_unlock(&queue->lock);
    }

  free(vector.entry);
  free(prof_genome);
  free(seq-1);
  fclose((FILE *)db.bases);
  Free_Profile_Reader(Pg_reader);
  return NULL;
}

static char *add_suffix(const char *path, const char *suffix)
{
  size_t plen = strlen(path);
  size_t slen = strlen(suffix);
  char *result = (char *)malloc(plen+slen+1);

  if (result == NULL)
    { perror("malloc output path");
      exit(1);
    }
  memcpy(result,path,plen);
  memcpy(result+plen,suffix,slen+1);
  return result;
}

static void print_periodic_stats(const Periodic_Stats *stats, FILE *out,
                                 const char *prefix, int cutoff)
{
  const char *name = prefix == NULL ? "seed" : prefix;

  if (stats == NULL || out == NULL || cutoff <= 0)
    return;
  fprintf(out,
          "PERIODIC_MASK_V1\t%s\tcomplete_motif_bp_gt=%d\t"
          "overlap=any\tstage=post_modimizer_selection\t"
          "hp_tracts=%llu\thp_mask_bases=%llu\t"
          "ds_tracts=%llu\tds_mask_bases=%llu\t"
          "ts_tracts=%llu\tts_mask_bases=%llu\t"
          "reads=%llu\tseeds_removed=%llu\n",
          name,cutoff,
          (unsigned long long)stats->periodic_tracts[0],
          (unsigned long long)stats->periodic_tract_bases[0],
          (unsigned long long)stats->periodic_tracts[1],
          (unsigned long long)stats->periodic_tract_bases[1],
          (unsigned long long)stats->periodic_tracts[2],
          (unsigned long long)stats->periodic_tract_bases[2],
          (unsigned long long)stats->periodic_reads,
          (unsigned long long)stats->periodic_seeds_removed);
}

int main(int argc, char *argv[])
{
  DAZZ_DB _db, *db = &_db;
  int NTHREADS = 1;

  // Optional flags precede the three positional arguments.
  { int i, j = 1;
    for (i = 1; i < argc; i++)
      { if (argv[i][0] == '-')
          { switch (argv[i][1])
              { case 'k': K = atoi(argv[i]+2); break;
                case 'T': NTHREADS = atoi(argv[i]+2); break;
                case 'l':
                  { char *end;
                    long value = strtol(argv[i]+2,&end,10);
                    if (argv[i][2] == '\0' || *end != '\0' ||
                        (value != -1 && value < 1) || value > INT32_MAX)
                      { fprintf(stderr,
                                "%s: invalid -l value '%s' "
                                "(expected -1 or a positive integer)\n",
                                argv[0],argv[i]+2);
                        return 1;
                      }
                    EXCLUDE_PERIODIC_BP = (int)value;
                  }
                  break;
                default:
                  fprintf(stderr,"%s: unknown option '%s'\n",argv[0],argv[i]);
                  return 1;
              }
          }
        else
          argv[j++] = argv[i];
      }
    argc = j;
  }

  if (argc != 4 || K < 1 || K > 32 || NTHREADS < 1)
    { fprintf(stderr,
              "Usage: %s [-T<int>] [-k<32>] [-l<int(-1)>] "
              "<db> <ref.prof> <out.kpos>\n",argv[0]);
      fprintf(stderr,"  -T<N>      : number of DB blocks processed in parallel (default 1)\n");
      fprintf(stderr,"  -k<N>      : k-mer length (<= 32)\n");
      fprintf(stderr,"  -l<N>      : exclude selected K-mers overlapping exact period-1/2/3 tracts whose complete-motif length is >N bp (-1 disables; default -1)\n");
      fprintf(stderr,"  <db>       : whole, partitioned DAZZLER dam/db\n");
      fprintf(stderr,"  <ref.prof> : FastK genome k-mer count profile\n");
      fprintf(stderr,"  <out.kpos> : block-indexed KPOS container; <out.kpos>.idx is also written\n");
      return 1;
    }

  fprintf(stderr,
          "Seed selector: mode=modimizer mod=%llu threshold=%llu "
          "window=disabled thinning=disabled\n",
          (unsigned long long)MODIMIZER_MOD,
          (unsigned long long)MODIMIZER_THRESHOLD);
  if (EXCLUDE_PERIODIC_BP > 0)
    fprintf(stderr,
            "Periodic seed exclusion: periods=1,2,3 "
            "complete_motif_bp>%d overlap=any "
            "stage=post_modimizer_selection\n",
            EXCLUDE_PERIODIC_BP);
  KMASK = K >= 32 ? UINT64_MAX : ((UINT64_C(1) << (2*K))-1);

  int db_type = Open_DB(argv[1],db);
  if (db_type < 0)
    return 1;

  Profile_Index *Pg = Open_Profiles(argv[2]);
  if (Pg == NULL)
    { fprintf(stderr,"%s: Cannot open genome profile '%s'\n",argv[0],argv[2]);
      Close_DB(db);
      return 1;
    }
  if (Pg->kmer != K)
    { fprintf(stderr,"%s: genome profile K=%d does not match K=%d\n",
              argv[0],Pg->kmer,K);
      Free_Profiles(Pg);
      Close_DB(db);
      return 1;
    }
  if (Pg->nreads != db->nreads)
    { fprintf(stderr,
              "%s: profile/DB read-count mismatch: DB=%d genome-profile=%d\n",
              argv[0],db->nreads,Pg->nreads);
      Free_Profiles(Pg);
      Close_DB(db);
      return 1;
    }

  DAZZ_STUB *stub = NULL;
  { char *stub_path;
    const char *stub_suffix = db_type > 0 ? ".dam" : ".db";
    size_t stub_path_len;
    size_t db_arg_len;
    size_t stub_suffix_len;
    int has_stub_suffix;
    int b;

    if (db->part != 0)
      { fprintf(stderr,"%s: requires a whole DB, not block %d\n",
                argv[0],db->part);
        Free_Profiles(Pg);
        Close_DB(db);
        return 1;
      }

    db_arg_len = strlen(argv[1]);
    stub_suffix_len = strlen(stub_suffix);
    has_stub_suffix =
        db_arg_len >= stub_suffix_len &&
        strcmp(argv[1]+db_arg_len-stub_suffix_len,stub_suffix) == 0;
    stub_path_len = db_arg_len+(has_stub_suffix ? 0 : stub_suffix_len)+1;
    stub_path = (char *)malloc(stub_path_len);
    if (stub_path == NULL)
      { perror("malloc DB stub path");
        Free_Profiles(Pg);
        Close_DB(db);
        return 1;
      }
    snprintf(stub_path,stub_path_len,"%s%s",argv[1],
             has_stub_suffix ? "" : stub_suffix);
    stub = Read_DB_Stub(stub_path,DB_STUB_BLOCKS);
    free(stub_path);
    if (stub == NULL || stub->nblocks < 1)
      { fprintf(stderr,"%s: requires a partitioned DB\n",argv[0]);
        if (stub != NULL) Free_DB_Stub(stub);
        Free_Profiles(Pg);
        Close_DB(db);
        return 1;
      }
    if (stub->ublocks[0] != 0 ||
        stub->ublocks[stub->nblocks] != db->nreads)
      { fprintf(stderr,
                "%s: DB block boundaries do not span all %d input reads\n",
                argv[0],db->nreads);
        Free_DB_Stub(stub);
        Free_Profiles(Pg);
        Close_DB(db);
        return 1;
      }
    for (b = 0; b <= stub->nblocks; b++)
      if (stub->ublocks[b] != stub->tblocks[b])
        { fprintf(stderr,
                  "%s: does not support differing trimmed/untrimmed block boundaries\n",
                  argv[0]);
          Free_DB_Stub(stub);
          Free_Profiles(Pg);
          Close_DB(db);
          return 1;
        }
    for (b = 0; b < stub->nblocks; b++)
      if ((uint64_t)(stub->ublocks[b+1]-stub->ublocks[b]) > CN_MAXVAL)
        { fprintf(stderr,
                  "%s: block %d exceeds the copy-tier rid limit (2^29 reads)\n",
                  argv[0],b+1);
          Free_DB_Stub(stub);
          Free_Profiles(Pg);
          Close_DB(db);
          return 1;
        }
  }

  { const char *output_path = argv[3];
    char *index_path = add_suffix(output_path,".idx");
    char *output_tmp = add_suffix(output_path,".tmp");
    char *index_tmp = add_suffix(index_path,".tmp");
    FILE *out = fopen(output_tmp,"wb");
    FILE *idx;
    WorkQueue queue;
    BlockResult *result;
    pthread_t *thread;
    uint64_t total_written = 0;
    Periodic_Stats total_periodic_stats = {0};
    int b, t;

    if (out == NULL)
      { perror(output_tmp);
        exit(1);
      }
    idx = fopen(index_tmp,"wb");
    if (idx == NULL)
      { perror(index_tmp);
        exit(1);
      }
    if (write_kpos_header_for_k(out,K) != 0 ||
        fwrite(&total_written,sizeof(uint64_t),1,idx) != 1)
      { fprintf(stderr,"%s: failed to initialize KPOS output\n",argv[0]);
        exit(1);
      }

    if (NTHREADS > stub->nblocks)
      NTHREADS = stub->nblocks;
    result = (BlockResult *)calloc((size_t)stub->nblocks,sizeof(BlockResult));
    thread = (pthread_t *)malloc((size_t)NTHREADS*sizeof(pthread_t));
    if (result == NULL || thread == NULL)
      { perror("malloc modimizer workers");
        exit(1);
      }

    memset(&queue,0,sizeof(queue));
    queue.db = db;
    queue.stub = stub;
    queue.genome = Pg;
    queue.result = result;
    queue.end_block = stub->nblocks;
    pthread_mutex_init(&queue.lock,NULL);
    pthread_cond_init(&queue.ready,NULL);
    pthread_cond_init(&queue.consumed,NULL);

    for (t = 0; t < NTHREADS; t++)
      if (pthread_create(thread+t,NULL,seed_worker,&queue) != 0)
        { fprintf(stderr,"%s: cannot create worker thread %d\n",argv[0],t+1);
          exit(1);
        }

    for (b = 0; b < stub->nblocks; b++)
      { BlockResult *block = result+b;

        pthread_mutex_lock(&queue.lock);
        while (!block->ready)
          pthread_cond_wait(&queue.ready,&queue.lock);
        pthread_mutex_unlock(&queue.lock);

        periodic_stats_add(&total_periodic_stats,&block->periodic_stats);
        if (block->len > 0 &&
            fwrite(block->entry,sizeof(KposEntry),block->len,out) != block->len)
          { fprintf(stderr,"%s: failed to write block %d to %s\n",
                    argv[0],b+1,output_tmp);
            exit(1);
          }
        total_written += (uint64_t)block->len;
        if (fwrite(&total_written,sizeof(uint64_t),1,idx) != 1)
          { fprintf(stderr,"%s: failed to write block index %s\n",
                    argv[0],index_tmp);
            exit(1);
          }
        fprintf(stderr,"Block %d seeds written: %llu\n",b+1,
                (unsigned long long)block->len);

        pthread_mutex_lock(&queue.lock);
        block->ready = 0;
        pthread_cond_broadcast(&queue.consumed);
        pthread_mutex_unlock(&queue.lock);
      }

    for (t = 0; t < NTHREADS; t++)
      pthread_join(thread[t],NULL);

    if (finish_kpos(out,total_written) != 0 ||
        fclose(out) != 0 || fclose(idx) != 0)
      { fprintf(stderr,"%s: failed to finalize KPOS output\n",argv[0]);
        exit(1);
      }
    if (rename(output_tmp,output_path) != 0)
      { perror(output_path);
        exit(1);
      }
    if (rename(index_tmp,index_path) != 0)
      { perror(index_path);
        exit(1);
      }

    fprintf(stderr,"Seeds written: %llu (%s, %s)\n",
            (unsigned long long)total_written,output_path,index_path);
    print_periodic_stats(&total_periodic_stats,stderr,"standalone_seed",
                         EXCLUDE_PERIODIC_BP);

    pthread_cond_destroy(&queue.consumed);
    pthread_cond_destroy(&queue.ready);
    pthread_mutex_destroy(&queue.lock);
    free(thread);
    free(result);
    free(index_path);
    free(output_tmp);
    free(index_tmp);
  }

  Free_Profiles(Pg);
  Free_DB_Stub(stub);
  Close_DB(db);
  return 0;
}
