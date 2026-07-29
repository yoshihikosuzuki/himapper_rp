#include <stdbool.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <DB.h>
#include "nthash.h"
#include "libfastk.h"
#include "copynum.h"

// #define DEBUG_HASH
// #define DEBUG_NW

int K = 32;          // k-mer length; maximum 32   [-k]
int WSIZE  = 100;    // sliding window width for (count, hash) minimizer   [-w]
int RSPACE = 100;    // minimum spacing between uniformly sparse seeds     [-s]
static const int MOD = 2003;   // hash modulus for tiebreaking (prime > WSIZE recommended)

// ASCII base → 2-bit value  (A=0, C=1, G=2, T=3), matching the DAZZ_DB
// numeric encoding used by himapper's tuple_thread.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

// Bitmask covering the low 2*K bits of a k-mer's 2-bit code; set in main() once K is
// known (K is capped at 32, so K==32 gives the full uint64 mask).
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

static void push_kpos(KposVector *vector, uint64_t code, int32_t rpos, uint32_t rid)
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
      vector->cap   = next;
    }

  vector->entry[vector->len].code = code;
  vector->entry[vector->len].rpos = rpos;
  vector->entry[vector->len].rid  = rid;
  vector->len += 1;
}

static int cmp_kpos(const void *a, const void *b)
{
  const KposEntry *ka = (const KposEntry *)a;
  const KposEntry *kb = (const KposEntry *)b;

  if (ka->code < kb->code) return -1;
  if (ka->code > kb->code) return 1;
  return 0;
}

static int write_kpos_header(FILE *out)
{
  char magic[4] = {'K', 'P', 'O', 'S'};
  uint8_t ver = 2, flags = 0, idx_bytes = sizeof(uint64_t), reserved = 0;
  int32_t Kval = K;
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

// Compare (ca_g, ca_r, ha) < (cb_g, cb_r, hb):
//   primary:   lower genome count  → more unique in genome   → better seed
//   secondary: higher reads count  → better covered in reads → better seed
//   tertiary:  lower hash % MOD    → stable tiebreaker
static inline int tuple_lt(uint16_t ca_g, uint16_t ca_r, uint64_t ha,
                            uint16_t cb_g, uint16_t cb_r, uint64_t hb)
{ if (ca_g != cb_g) return ca_g < cb_g;
  if (ca_r != cb_r) return ca_r > cb_r;   // higher reads count wins
  return ha < hb;
}

// Candidate for uniformly sparse thinning — stores sort-key fields so the
// qsort comparator does not need external array access.
typedef struct {
  int      pos;
  int      nw_val;
  uint16_t g_count;
  uint16_t r_count;
  uint64_t hash;
} Candidate;

// Compare candidates by descending h_nw priority:
//   1. nw      descending (more windows → higher priority)
//   2. g_count ascending  (rarer in genome → better seed)
//   3. r_count descending (more covered by reads → better seed)
//   4. hash    ascending  (stable tiebreaker)
static int cand_cmp(const void *a, const void *b)
{ const Candidate *ca = (const Candidate *)a;
  const Candidate *cb = (const Candidate *)b;
  if (ca->nw_val != cb->nw_val)
    return (cb->nw_val > ca->nw_val) ? 1 : -1;
  if (ca->g_count != cb->g_count)
    return (ca->g_count > cb->g_count) ? 1 : -1;
  if (ca->r_count != cb->r_count)
    return (cb->r_count > ca->r_count) ? 1 : -1;
  if (ca->hash != cb->hash)
    return (ca->hash > cb->hash) ? 1 : -1;
  return 0;
}

// Select seeds using a sliding-window (genome_count, -reads_count, hash) minimizer.
//
// For every window of WSIZE consecutive k-mers the one with the lexicographically
// smallest (prof_genome[i], -prof_reads[i], hVal[i]) tuple is selected.
// The result is the union of all per-window selections.
//
//   seq         – pointer to position (K-1) in the read
//   plen        – number of k-mers (rlen - K + 1)
//   K           – k-mer length
//   prof_genome – FastK uint16 genome k-mer counts, positions 0..plen-1
//   prof_reads  – FastK uint16 reads  k-mer counts, positions 0..plen-1
//   read_id     – 0-based read index
//   out         – block-local KPOS vector
//
// Returns the number of entries appended.
int find_seeds(const char    *seq,
               const int      plen,
               const int      K,
               const uint16  *prof_genome,
               const uint16  *prof_reads,
               const int      read_id,
               KposVector    *out)
{
  if (plen <= 0) return 0;

  // Step 1: compute all canonical hash values.
  uint64_t *hvals = (uint64_t *) malloc((size_t)plen * sizeof(uint64_t));
  if (!hvals) { perror("malloc hvals"); exit(1); }
  { uint64_t fhVal = 0, rhVal = 0;
    hvals[0] = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);
    for (int i = 1; i < plen; i++)
      hvals[i] = NTC64_c((unsigned char)seq[i - K], (unsigned char)seq[i],
                         K, &fhVal, &rhVal);
    // Reduce to [0, MOD) so that hash ties are common — enables robust winnowing
    // to reliably favour the leftmost (earliest) occurrence within a window.
    for (int i = 0; i < plen; i++)
      hvals[i] %= MOD;
  }

  // Step 2: sliding-window minimum with a monotone deque.
  //   The deque stores indices in increasing position order.
  //   Front = index with the minimum (count, hash) in the current window.
  //   We pop from the back any index j whose tuple >= the incoming tuple
  //   (it can never win a future window while the new element is around).
  int  *dq      = (int  *) malloc((size_t)plen * sizeof(int));
  bool *is_seed = (bool *) calloc((size_t)plen,  sizeof(bool));
  int  *nw      = (int  *) calloc((size_t)plen,  sizeof(int));   // # of windows for which position i is minimizer
  if (!dq || !is_seed || !nw) { perror("malloc deque/is_seed/nw"); exit(1); }

  int dq_head = 0, dq_tail = 0;
  for (int i = 0; i < plen; i++)
    { // Robust winnowing: evict back elements whose tuple is STRICTLY GREATER
      // than the current tuple.  On a tie we keep the existing (leftmost) entry,
      // so ties are broken in favour of the element already inside the window.
      if (prof_genome[i] != 0)
        { while (dq_head < dq_tail &&
                 tuple_lt(prof_genome[i],             prof_reads[i],             hvals[i],
                          prof_genome[dq[dq_tail-1]], prof_reads[dq[dq_tail-1]], hvals[dq[dq_tail-1]]))
            dq_tail--;
          dq[dq_tail++] = i;
        }

      // Evict front elements that have slid out of the window.  With g==0 positions
      // skipped the deque can be empty, so guard against dq_head == dq_tail.
      while (dq_head < dq_tail && dq[dq_head] < i - WSIZE + 1)
        dq_head++;

      // Only consider complete windows of exactly WSIZE k-mers (standard minimizer
      // definition).  Emit only if the window holds a g>=1 candidate (deque non-empty);
      // an all-g==0 window produces no seed.
      if (i >= WSIZE - 1 && dq_head < dq_tail)
        { is_seed[dq[dq_head]] = true;
          nw[dq[dq_head]]++;
        }
#ifdef DEBUG_HASH
      if (i >= WSIZE - 1 && dq_head < dq_tail)
      { int win_start = i - WSIZE + 1;
        int min_pos   = dq[dq_head];
        fprintf(stderr, "WINDOW read=%d [%d..%d]  minimizer=@%d g=%u r=%u h=%lu\n",
                read_id, win_start, i, min_pos,
                (unsigned)prof_genome[min_pos], (unsigned)prof_reads[min_pos], hvals[min_pos]);
        for (int j = win_start; j <= i; j++) {
          const char *wkp = seq - K + 1 + j;
          char wkstr[65];
          for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
          wkstr[K] = '\0';
          fprintf(stderr, "  j=%-6d g=%-5u r=%-5u h=%-4lu  %s%s\n",
                  j,
                  (unsigned)prof_genome[j], (unsigned)prof_reads[j], hvals[j],
                  wkstr,
                  (j == min_pos) ? "  <-- MIN" : "");
        }
      }
#endif
    }

  // Step 3: Uniformly sparse thinning.
  // Sort candidates by descending h_nw = (NW, -h_minimax) priority, then
  // greedily accept only those whose position is ≥ RSPACE from every
  // already-accepted seed.
  { int ncand = 0;
    for (int i = 0; i < plen; i++)
      if (is_seed[i]) ncand++;

    if (ncand > 1)
      { Candidate *cands = (Candidate *) malloc((size_t)ncand * sizeof(Candidate));
        if (!cands) { perror("malloc cands"); exit(1); }
        { int ci = 0;
          for (int i = 0; i < plen; i++)
            if (is_seed[i])
              { cands[ci].pos     = i;
                cands[ci].nw_val  = nw[i];
                cands[ci].g_count = prof_genome[i];
                cands[ci].r_count = prof_reads[i];
                cands[ci].hash    = hvals[i];
                ci++;
              }
        }
        qsort(cands, (size_t)ncand, sizeof(Candidate), cand_cmp);

        // Reset is_seed; re-mark only accepted positions.
        memset(is_seed, 0, (size_t)plen * sizeof(bool));

        // Maintain accepted positions in sorted order (by position) for
        // O(log n) nearest-neighbour distance checks.
        int *acc = (int *) malloc((size_t)ncand * sizeof(int));
        if (!acc) { perror("malloc acc"); exit(1); }
        int nacc = 0;

        for (int ci = 0; ci < ncand; ci++)
          { int p = cands[ci].pos;
            // Binary search for insertion point in acc[].
            int lo = 0, hi = nacc;
            while (lo < hi)
              { int mid = (lo + hi) / 2;
                if (acc[mid] < p) lo = mid + 1;
                else              hi = mid;
              }
            // Check left and right neighbours in accepted set.
            int ok = 1;
            if (lo > 0    && p - acc[lo-1] < RSPACE) ok = 0;
            if (lo < nacc && acc[lo] - p   < RSPACE) ok = 0;
            if (ok)
              { memmove(acc + lo + 1, acc + lo,
                        (size_t)(nacc - lo) * sizeof(int));
                acc[lo] = p;
                nacc++;
                is_seed[p] = true;
              }
          }
        free(acc);
        free(cands);
      }
  }

#ifdef DEBUG_NW
  // Print per-position window-count; mark final selected seeds with *****.
  fprintf(stderr, "NW read=%d plen=%d\n", read_id, plen);
  for (int i = 0; i < plen; i++)
    { const char *wkp = seq - K + 1 + i;
      char wkstr[65];
      for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
      wkstr[K] = '\0';
      fprintf(stderr, "  i=%-6d nw=%-5d g=%-5u r=%-5u h=%-4lu  %s  %s%s\n",
              i, nw[i],
              (unsigned)prof_genome[i], (unsigned)prof_reads[i], hvals[i],
              wkstr, nw[i] ? "*" : "", is_seed[i] ? "****" : "");
    }
#endif

  // Step 4: append selected seeds to this block's KPOS run.
  int written = 0;
  for (int i = 0; i < plen; i++)
    { if (!is_seed[i]) continue;
      const char *kp = seq - K + 1 + i;
      uint64_t c = 0;
      for (int x = 0; x < K; x++)
        c = (c << 2) | base2bit[(unsigned char)kp[x]];
      c &= KMASK;
      int32_t rpos = (int32_t)i;
      // Every seed's rid carries the genome copy-tier in its top 3 bits (see copynum.h):
      // the low 29 bits stay ref/read id. Strip the tier with CN_VAL wherever it needs the id.
      // read_id < CN_MAXVAL (2^29) is enforced in main(), so we pack unconditionally.
      uint32_t rid = CN_PACK(read_id, cn_tier_of(prof_genome[i]));
      push_kpos(out, c, rpos, rid);
      written++;
    }

  free(is_seed);
  free(nw);
  free(dq);
  free(hvals);
  return written;
}

typedef struct
  { KposEntry *entry;
    size_t     len;
    int        ready;
  } BlockResult;

typedef struct
  { DAZZ_DB       *db;
    DAZZ_STUB     *stub;
    Profile_Index *genome;
    Profile_Index *reads;
    BlockResult   *result;
    int            next_block;
    pthread_mutex_t lock;
    pthread_cond_t  ready;
    pthread_cond_t  consumed;
  } WorkQueue;

static void grow_profiles(int need, int *capacity,
                          uint16 **prof_genome, uint16 **prof_reads)
{
  int next = (int)(1.2 * need) + 1000;
  uint16 *genome =
      (uint16 *)realloc(*prof_genome, (size_t)next * sizeof(uint16));
  uint16 *reads =
      (uint16 *)realloc(*prof_reads, (size_t)next * sizeof(uint16));

  if (genome == NULL || reads == NULL)
    { perror("realloc profile");
      exit(1);
    }
  *prof_genome = genome;
  *prof_reads  = reads;
  *capacity    = next;
}

static FILE *open_worker_bases(const DAZZ_DB *db)
{
  size_t len = strlen(db->path);
  char *path = (char *)malloc(len + 5);
  FILE *bases;

  if (path == NULL)
    { perror("malloc .bps path");
      exit(1);
    }
  snprintf(path, len + 5, "%s.bps", db->path);
  bases = fopen(path, "rb");
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
  Profile_Index *Pg = Clone_Profiles(queue->genome);
  Profile_Index *Pr = Clone_Profiles(queue->reads);
  char *seq;
  uint16 *prof_genome;
  uint16 *prof_reads;
  int pmax = 20000;
  KposVector vector = {NULL, 0, 0};

  db.bases = open_worker_bases(queue->db);
  seq = New_Read_Buffer(&db);
  prof_genome = (uint16 *)malloc((size_t)pmax * sizeof(uint16));
  prof_reads  = (uint16 *)malloc((size_t)pmax * sizeof(uint16));
  if (prof_genome == NULL || prof_reads == NULL)
    { perror("malloc profile");
      exit(1);
    }

  while (1)
    { int block;
      int first, last, id;
      BlockResult *result;

      pthread_mutex_lock(&queue->lock);
      block = queue->next_block++;
      pthread_mutex_unlock(&queue->lock);
      if (block >= queue->stub->nblocks)
        break;

      first = queue->stub->ublocks[block];
      last  = queue->stub->ublocks[block + 1];
      vector.len = 0;

      for (id = first; id < last; id++)
        { DAZZ_READ *r = db.reads + id;
          int plen = r->rlen - K + 1;
          int64 fetched;

          Load_Read(&db, id, seq, 2);

          fetched = Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
          if (fetched > pmax)
            { grow_profiles((int)fetched, &pmax, &prof_genome, &prof_reads);
              Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
            }

          fetched = Fetch_Profile(Pr, (int64)id, pmax, prof_reads);
          if (fetched > pmax)
            { grow_profiles((int)fetched, &pmax, &prof_genome, &prof_reads);
              Fetch_Profile(Pr, (int64)id, pmax, prof_reads);
            }

          if (plen > 0)
            find_seeds(seq + K - 1, plen, K, prof_genome, prof_reads,
                       id - first, &vector);
        }

      if (vector.len > 1)
        qsort(vector.entry, vector.len, sizeof(KposEntry), cmp_kpos);

      result = queue->result + block;
      pthread_mutex_lock(&queue->lock);
      result->entry = vector.entry;
      result->len   = vector.len;
      result->ready = 1;
      pthread_cond_broadcast(&queue->ready);
      while (result->ready)
        pthread_cond_wait(&queue->consumed, &queue->lock);
      pthread_mutex_unlock(&queue->lock);
    }

  free(vector.entry);
  free(prof_genome);
  free(prof_reads);
  free(seq - 1);
  fclose((FILE *)db.bases);
  Free_Profiles(Pg);
  Free_Profiles(Pr);
  return NULL;
}

static char *add_suffix(const char *path, const char *suffix)
{
  size_t plen = strlen(path);
  size_t slen = strlen(suffix);
  char *result = (char *)malloc(plen + slen + 1);

  if (result == NULL)
    { perror("malloc output path");
      exit(1);
    }
  memcpy(result, path, plen);
  memcpy(result + plen, suffix, slen + 1);
  return result;
}

int main(int argc, char *argv[]) {
  DAZZ_DB _db, *db = &_db;
  int NTHREADS = 1;
  // Optional flags precede the 4 positional args.
  { int i, j = 1;
    for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'k': K      = atoi(argv[i]+2); break;
          case 'w': WSIZE  = atoi(argv[i]+2); break;
          case 's': RSPACE = atoi(argv[i]+2); break;
          case 'T': NTHREADS = atoi(argv[i]+2); break;
          default:
            fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[i]);
            return 1;
        }
      } else
        argv[j++] = argv[i];
    }
    argc = j;
  }
  if (argc != 5 || K < 1 || K > 32 || WSIZE < 1 || RSPACE < 0 ||
      NTHREADS < 1) {
    fprintf(stderr, "Usage: %s [-T<int>] [-k<32>] [-w<100>] [-s<100>] <db> <ref.prof> <ref.reads.prof> <out.kpos>\n", argv[0]);
    fprintf(stderr, "  -T<N>           : number of DB blocks processed in parallel (default 1)\n");
    fprintf(stderr, "  -k<N>           : k-mer length (<= 32)\n");
    fprintf(stderr, "  -w<N>           : minimizer window width  (larger => sparser seeds)\n");
    fprintf(stderr, "  -s<N>           : min spacing of uniform sparse seeds (larger => sparser)\n");
    fprintf(stderr, "  <db>            : whole, partitioned DAZZ reference dam/db\n");
    fprintf(stderr, "  <ref.prof>      : FastK genome k-mer count profile\n");
    fprintf(stderr, "  <ref.reads.prof>: FastK reads  k-mer count profile\n");
    fprintf(stderr, "  <out.kpos>      : block-indexed KPOS container; <out.kpos>.idx is also written\n");
    return 1;
  }
  KMASK = (K >= 32) ? UINT64_MAX : ((UINT64_C(1) << (2 * K)) - 1);
  int db_type = Open_DB(argv[1], db);
  if (db_type < 0)
    return 1;

  // Open genome profile (argv[2]).
  Profile_Index *Pg = Open_Profiles(argv[2]);
  if (Pg == NULL)
    { fprintf(stderr, "%s: Cannot open genome profile '%s'\n", argv[0], argv[2]);
      return 1;
    }
  if (Pg->kmer != K)
    { fprintf(stderr, "%s: genome profile K=%d does not match K=%d\n", argv[0], Pg->kmer, K);
      Free_Profiles(Pg);
      return 1;
    }

  // Open reads profile (argv[3]).
  Profile_Index *Pr = Open_Profiles(argv[3]);
  if (Pr == NULL)
    { fprintf(stderr, "%s: Cannot open reads profile '%s'\n", argv[0], argv[3]);
      Free_Profiles(Pg);
      return 1;
    }
  if (Pr->kmer != K)
    { fprintf(stderr, "%s: reads profile K=%d does not match K=%d\n", argv[0], Pr->kmer, K);
      Free_Profiles(Pg); Free_Profiles(Pr);
      return 1;
    }
  if (Pg->nreads != db->nreads || Pr->nreads != db->nreads)
    { fprintf(stderr,
              "%s: profile/DB read-count mismatch: DB=%d genome-profile=%d reads-profile=%d\n",
              argv[0], db->nreads, Pg->nreads, Pr->nreads);
      Free_Profiles(Pg); Free_Profiles(Pr);
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
      { fprintf(stderr, "%s: requires a whole DB, not block %d\n",
                argv[0], db->part);
        Free_Profiles(Pg); Free_Profiles(Pr);
        Close_DB(db);
        return 1;
      }

    // db->path names the auxiliary-file root and can include a hidden-file
    // component (for example "./.reads").  The DB stub is the user-visible
    // argv[1] path instead, with its .db/.dam suffix added when omitted.
    db_arg_len = strlen(argv[1]);
    stub_suffix_len = strlen(stub_suffix);
    has_stub_suffix =
        db_arg_len >= stub_suffix_len &&
        strcmp(argv[1] + db_arg_len - stub_suffix_len, stub_suffix) == 0;
    stub_path_len =
        db_arg_len + (has_stub_suffix ? 0 : stub_suffix_len) + 1;
    stub_path = (char *)malloc(stub_path_len);
    if (stub_path == NULL)
      { perror("malloc DB stub path");
        Free_Profiles(Pg); Free_Profiles(Pr);
        Close_DB(db);
        return 1;
      }
    snprintf(stub_path, stub_path_len, "%s%s", argv[1],
             has_stub_suffix ? "" : stub_suffix);
    stub = Read_DB_Stub(stub_path, DB_STUB_BLOCKS);
    free(stub_path);
    if (stub == NULL || stub->nblocks < 1)
      { fprintf(stderr, "%s: requires a partitioned DB\n", argv[0]);
        if (stub != NULL) Free_DB_Stub(stub);
        Free_Profiles(Pg); Free_Profiles(Pr);
        Close_DB(db);
        return 1;
      }
    if (stub->ublocks[0] != 0 ||
        stub->ublocks[stub->nblocks] != db->nreads)
      { fprintf(stderr,
                "%s: DB block boundaries do not span all %d input reads\n",
                argv[0], db->nreads);
        Free_DB_Stub(stub);
        Free_Profiles(Pg); Free_Profiles(Pr);
        Close_DB(db);
        return 1;
      }
    for (b = 0; b <= stub->nblocks; b++)
      if (stub->ublocks[b] != stub->tblocks[b])
        { fprintf(stderr,
                  "%s: does not support differing trimmed/untrimmed block boundaries\n",
                  argv[0]);
          Free_DB_Stub(stub);
          Free_Profiles(Pg); Free_Profiles(Pr);
          Close_DB(db);
          return 1;
        }
    for (b = 0; b < stub->nblocks; b++)
      if ((uint64_t)(stub->ublocks[b + 1] - stub->ublocks[b]) > CN_MAXVAL)
        { fprintf(stderr,
                  "%s: block %d exceeds the copy-tier rid limit (2^29 reads)\n",
                  argv[0], b + 1);
          Free_DB_Stub(stub);
          Free_Profiles(Pg); Free_Profiles(Pr);
          Close_DB(db);
          return 1;
        }
  }

  {
    const char *output_path = argv[4];
    char *index_path = add_suffix(output_path, ".idx");
    char *output_tmp = add_suffix(output_path, ".tmp");
    char *index_tmp = add_suffix(index_path, ".tmp");
    FILE *out = fopen(output_tmp, "wb");
    FILE *idx;
    WorkQueue queue;
    BlockResult *result;
    pthread_t *thread;
    uint64_t total_written = 0;
    int b, t;

    if (out == NULL)
      { perror(output_tmp);
        exit(1);
      }
    idx = fopen(index_tmp, "wb");
    if (idx == NULL)
      { perror(index_tmp);
        exit(1);
      }
    if (write_kpos_header(out) != 0 ||
        fwrite(&total_written, sizeof(uint64_t), 1, idx) != 1)
      { fprintf(stderr, "%s: failed to initialize KPOS output\n", argv[0]);
        exit(1);
      }

    if (NTHREADS > stub->nblocks)
      NTHREADS = stub->nblocks;
    result = (BlockResult *)calloc((size_t)stub->nblocks,
                                    sizeof(BlockResult));
    thread = (pthread_t *)malloc((size_t)NTHREADS * sizeof(pthread_t));
    if (result == NULL || thread == NULL)
      { perror("malloc minimaxmer workers");
        exit(1);
      }

    queue.db         = db;
    queue.stub       = stub;
    queue.genome     = Pg;
    queue.reads      = Pr;
    queue.result     = result;
    queue.next_block = 0;
    pthread_mutex_init(&queue.lock, NULL);
    pthread_cond_init(&queue.ready, NULL);
    pthread_cond_init(&queue.consumed, NULL);

    for (t = 0; t < NTHREADS; t++)
      if (pthread_create(thread + t, NULL, seed_worker, &queue) != 0)
        { fprintf(stderr, "%s: cannot create worker thread %d\n",
                  argv[0], t + 1);
          exit(1);
        }

    for (b = 0; b < stub->nblocks; b++)
      { BlockResult *block = result + b;

        pthread_mutex_lock(&queue.lock);
        while (!block->ready)
          pthread_cond_wait(&queue.ready, &queue.lock);
        pthread_mutex_unlock(&queue.lock);

        if (block->len > 0 &&
            fwrite(block->entry, sizeof(KposEntry), block->len, out) !=
                block->len)
          { fprintf(stderr, "%s: failed to write block %d to %s\n",
                    argv[0], b + 1, output_tmp);
            exit(1);
          }
        total_written += (uint64_t)block->len;
        if (fwrite(&total_written, sizeof(uint64_t), 1, idx) != 1)
          { fprintf(stderr, "%s: failed to write block index %s\n",
                    argv[0], index_tmp);
            exit(1);
          }
        fprintf(stderr, "Block %d seeds written: %llu\n", b + 1,
                (unsigned long long)block->len);

        pthread_mutex_lock(&queue.lock);
        block->ready = 0;
        pthread_cond_broadcast(&queue.consumed);
        pthread_mutex_unlock(&queue.lock);
      }

    for (t = 0; t < NTHREADS; t++)
      pthread_join(thread[t], NULL);

    if (finish_kpos(out, total_written) != 0 ||
        fclose(out) != 0 || fclose(idx) != 0)
      { fprintf(stderr, "%s: failed to finalize KPOS output\n", argv[0]);
        exit(1);
      }
    if (rename(output_tmp, output_path) != 0)
      { perror(output_path);
        exit(1);
      }
    if (rename(index_tmp, index_path) != 0)
      { perror(index_path);
        exit(1);
      }

    fprintf(stderr, "Seeds written: %llu (%s, %s)\n",
            (unsigned long long)total_written, output_path, index_path);

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
  Free_Profiles(Pr);
  Free_DB_Stub(stub);
  Close_DB(db);

  return 0;
}
