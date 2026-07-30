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

// Experimental build-time controls.  The defaults preserve the historical
// minimaxmer selector exactly.  They are compile-time rather than run-time
// options so production invocations cannot accidentally change seed semantics.
//
// MINIMAX_HASH_MOD=0 uses the full canonical 64-bit ntHash.  With
// MINIMAX_HASH_REFINE=1, the key is (hash % MINIMAX_HASH_MOD, full_hash);
// this preserves every non-tied legacy ordering and only resolves modulo ties.
// MINIMAX_ROBUST_WINNOW=1 implements Schleimer/Winnowmap robust winnowing:
// retain the previous winner while it remains a minimum, otherwise choose the
// rightmost minimum in the new complete window.
#ifndef MINIMAX_HASH_MOD
#define MINIMAX_HASH_MOD 2003
#endif
#ifndef MINIMAX_HASH_REFINE
#define MINIMAX_HASH_REFINE 0
#endif
#ifndef MINIMAX_ROBUST_WINNOW
#define MINIMAX_ROBUST_WINNOW 0
#endif
#ifndef MINIMAX_TIE_AUDIT
#define MINIMAX_TIE_AUDIT 0
#endif
#if MINIMAX_HASH_REFINE && MINIMAX_HASH_MOD == 0
#error "MINIMAX_HASH_REFINE requires a non-zero MINIMAX_HASH_MOD"
#endif
#if MINIMAX_TIE_AUDIT && \
    (MINIMAX_HASH_MOD != 2003 || MINIMAX_HASH_REFINE || MINIMAX_ROBUST_WINNOW)
#error "MINIMAX_TIE_AUDIT is defined only for the legacy h%2003 selector"
#endif

static const uint64_t MOD = MINIMAX_HASH_MOD;

// ASCII base → 2-bit value  (A=0, C=1, G=2, T=3), matching the DAZZ_DB
// numeric encoding used by himapper's tuple_thread.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

// Bitmask covering the low 2*K bits of a k-mer's 2-bit code; set in main() once K is
// known (K is capped at 32, so K==32 gives the full uint64 mask).
static uint64_t KMASK;

#if MINIMAX_TIE_AUDIT
#define TIE_AUDIT_MAX_W 100

typedef struct
  { uint64_t sequences;
    uint64_t sequences_with_complete_window;
    uint64_t sequences_with_full_tie;
    uint64_t complete_windows;
    uint64_t eligible_windows;
    uint64_t empty_windows;
    uint64_t profile_tie_windows;
    uint64_t full_tie_windows;
    uint64_t full_tie_repeat_only_windows;
    uint64_t full_tie_mod_collision_windows;
    uint64_t full_tie_runs;
    uint64_t full_tied_occurrences;
    uint64_t profile_min_multiplicity[TIE_AUDIT_MAX_W + 1];
    uint64_t full_min_multiplicity[TIE_AUDIT_MAX_W + 1];
  } TieStats;

static void tie_stats_add(TieStats *dst, const TieStats *src)
{
  int i;

  dst->sequences += src->sequences;
  dst->sequences_with_complete_window +=
      src->sequences_with_complete_window;
  dst->sequences_with_full_tie += src->sequences_with_full_tie;
  dst->complete_windows += src->complete_windows;
  dst->eligible_windows += src->eligible_windows;
  dst->empty_windows += src->empty_windows;
  dst->profile_tie_windows += src->profile_tie_windows;
  dst->full_tie_windows += src->full_tie_windows;
  dst->full_tie_repeat_only_windows +=
      src->full_tie_repeat_only_windows;
  dst->full_tie_mod_collision_windows +=
      src->full_tie_mod_collision_windows;
  dst->full_tie_runs += src->full_tie_runs;
  dst->full_tied_occurrences += src->full_tied_occurrences;
  for (i = 0; i <= TIE_AUDIT_MAX_W; i++)
    { dst->profile_min_multiplicity[i] +=
          src->profile_min_multiplicity[i];
      dst->full_min_multiplicity[i] += src->full_min_multiplicity[i];
    }
}
#endif

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
//   tertiary:  lower hash key       → stable tiebreaker
// In refinement builds only, the full hash is the final collision resolver.
static inline int position_tuple_lt(int a, int b,
                                    const uint16_t *prof_genome,
                                    const uint16_t *prof_reads,
                                    const uint64_t *hvals,
                                    const uint64_t *full_hvals)
{ if (prof_genome[a] != prof_genome[b])
    return prof_genome[a] < prof_genome[b];
  if (prof_reads[a] != prof_reads[b])
    return prof_reads[a] > prof_reads[b];   // higher reads count wins
  if (hvals[a] != hvals[b])
    return hvals[a] < hvals[b];
#if MINIMAX_HASH_REFINE
  return full_hvals[a] < full_hvals[b];
#else
  (void)full_hvals;
  return 0;
#endif
}

static inline int position_tuple_equal(int a, int b,
                                       const uint16_t *prof_genome,
                                       const uint16_t *prof_reads,
                                       const uint64_t *hvals,
                                       const uint64_t *full_hvals)
{ if (prof_genome[a] != prof_genome[b] ||
      prof_reads[a] != prof_reads[b] ||
      hvals[a] != hvals[b])
    return 0;
#if MINIMAX_HASH_REFINE
  return full_hvals[a] == full_hvals[b];
#else
  (void)full_hvals;
  return 1;
#endif
}

#if MINIMAX_TIE_AUDIT
static inline int profile_tuple_lt(int a, int b,
                                   const uint16_t *prof_genome,
                                   const uint16_t *prof_reads)
{
  if (prof_genome[a] != prof_genome[b])
    return prof_genome[a] < prof_genome[b];
  return prof_reads[a] > prof_reads[b];
}

static inline int profile_tuple_equal(int a, int b,
                                      const uint16_t *prof_genome,
                                      const uint16_t *prof_reads)
{
  return prof_genome[a] == prof_genome[b] &&
         prof_reads[a] == prof_reads[b];
}

// Exact canonical 2-bit codes distinguish repeated occurrences of the same
// canonical k-mer from different k-mers that collide after h%2003.  K<=32
// makes the representation collision-free.
static void fill_canonical_codes(const char *seq, int plen, int kmer,
                                 uint64_t *codes)
{
  const char *bases = seq - kmer + 1;
  uint64_t forward = 0;
  uint64_t reverse = 0;
  int x;

  for (x = 0; x < kmer; x++)
    { uint64_t b = base2bit[(unsigned char)bases[x]];
      forward = (forward << 2) | b;
      reverse |= (UINT64_C(3) - b) << (2 * x);
    }
  codes[0] = forward < reverse ? forward : reverse;

  for (x = 1; x < plen; x++)
    { uint64_t incoming =
          base2bit[(unsigned char)bases[x + kmer - 1]];
      forward = ((forward << 2) | incoming) & KMASK;
      reverse = (reverse >> 2) |
                ((UINT64_C(3) - incoming) << (2 * (kmer - 1)));
      codes[x] = forward < reverse ? forward : reverse;
    }
}
#endif

// Candidate for uniformly sparse thinning — stores sort-key fields so the
// qsort comparator does not need external array access.
typedef struct {
  int      pos;
  int      nw_val;
  uint16_t g_count;
  uint16_t r_count;
  uint64_t hash;
} Candidate;

#if MINIMAX_HASH_REFINE
// Keep Candidate byte-layout identical to the historical build.  A larger
// element size can change qsort's ordering inside comparator-equal groups and
// would confound collision refinement with an unrelated spacing change.
static _Thread_local const uint64_t *candidate_full_hvals;
#endif

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
#if MINIMAX_HASH_REFINE
  { uint64_t fa = candidate_full_hvals[ca->pos];
    uint64_t fb = candidate_full_hvals[cb->pos];
    if (fa != fb)
      return (fa > fb) ? 1 : -1;
  }
#endif
  return 0;
}

static inline int seed_bit_is_set(const uint64_t *bits, size_t pos)
{
  return (int)((bits[pos >> 6] >> (pos & 63)) & UINT64_C(1));
}

static inline void seed_bit_set(uint64_t *bits, size_t pos)
{
  bits[pos >> 6] |= UINT64_C(1) << (pos & 63);
}

// Set the inclusive bit range [lo, hi].  Keep the first/last-word cases
// separate so no expression ever shifts a uint64_t by 64 bits.
static void seed_bit_set_range(uint64_t *bits, size_t lo, size_t hi)
{
  size_t first = lo >> 6;
  size_t last  = hi >> 6;
  unsigned int lbit = (unsigned int)(lo & 63);
  unsigned int hbit = (unsigned int)(hi & 63);
  uint64_t left_mask  = UINT64_MAX << lbit;
  uint64_t right_mask = UINT64_MAX >> (63U - hbit);

  if (first == last)
    { bits[first] |= left_mask & right_mask;
      return;
    }

  bits[first] |= left_mask;
  if (last > first + 1)
    memset(bits + first + 1, 0xff, (last - first - 1) * sizeof(uint64_t));
  bits[last] |= right_mask;
}

// Count how many complete windows are assigned to each position.  Keeping this
// selector separate from hashing and spacing lets the experiment validate both
// deque policies against a simple O(nW) oracle without changing KPOS output.
static void count_window_minima(const char *seq,
                                int plen,
                                int kmer,
                                const uint16_t *prof_genome,
                                const uint16_t *prof_reads,
                                const uint64_t *hvals,
                                const uint64_t *full_hvals,
                                int read_id,
                                int *dq,
                                int *nw
#if MINIMAX_TIE_AUDIT
                                , const uint64_t *canonical_codes,
                                int *profile_dq,
                                TieStats *tie_stats
#endif
                                )
{
#ifndef DEBUG_HASH
  (void)seq;
  (void)kmer;
  (void)read_id;
#endif

  int dq_head = 0, dq_tail = 0;
#if MINIMAX_TIE_AUDIT
  int profile_head = 0, profile_tail = 0;
  int previous_window_tied = 0;
  int sequence_has_tie = 0;

  tie_stats->sequences++;
  if (plen >= WSIZE)
    tie_stats->sequences_with_complete_window++;
#endif
#if MINIMAX_ROBUST_WINNOW
  int previous_min = -1;
#endif
  for (int i = 0; i < plen; i++)
    {
      if (prof_genome[i] != 0)
        {
#if MINIMAX_TIE_AUDIT
          while (profile_head < profile_tail &&
                 profile_tuple_lt(i, profile_dq[profile_tail-1],
                                  prof_genome, prof_reads))
            profile_tail--;
          profile_dq[profile_tail++] = i;
#endif
#if MINIMAX_ROBUST_WINNOW
          // Keep a rightmost representative in the deque: discard a back
          // element whenever it is >= the incoming tuple.  The previous
          // window winner is tracked separately below, so an equal incoming
          // tuple does not dislodge it while that winner remains active.
          while (dq_head < dq_tail &&
                 !position_tuple_lt(dq[dq_tail-1], i,
                                    prof_genome, prof_reads,
                                    hvals, full_hvals))
            dq_tail--;
#else
          // Historical selector: discard only strictly greater back elements,
          // retaining every equal tuple so the deque front is the leftmost.
          while (dq_head < dq_tail &&
                 position_tuple_lt(i, dq[dq_tail-1],
                                   prof_genome, prof_reads,
                                   hvals, full_hvals))
            dq_tail--;
#endif
          dq[dq_tail++] = i;
        }

      // Evict front elements that have slid out of the window.  With g==0
      // positions skipped the deque can be empty.
      while (dq_head < dq_tail && dq[dq_head] < i - WSIZE + 1)
        dq_head++;
#if MINIMAX_TIE_AUDIT
      while (profile_head < profile_tail &&
             profile_dq[profile_head] < i - WSIZE + 1)
        profile_head++;

      if (i >= WSIZE - 1)
        tie_stats->complete_windows++;
#endif

      // Only complete windows contribute.  An all-g==0 window has no winner.
      if (i >= WSIZE - 1 && dq_head < dq_tail)
        {
          int min_pos = dq[dq_head];
#if MINIMAX_TIE_AUDIT
          int profile_min_pos = profile_dq[profile_head];
          int profile_multiplicity = 1;
          int full_multiplicity = 1;
          int q;

          tie_stats->eligible_windows++;
          for (q = profile_head + 1; q < profile_tail; q++)
            { if (!profile_tuple_equal(profile_dq[q], profile_min_pos,
                                       prof_genome, prof_reads))
                break;
              profile_multiplicity++;
            }
          for (q = dq_head + 1; q < dq_tail; q++)
            { if (!position_tuple_equal(dq[q], min_pos,
                                        prof_genome, prof_reads,
                                        hvals, full_hvals))
                break;
              full_multiplicity++;
            }

          tie_stats->profile_min_multiplicity[profile_multiplicity]++;
          tie_stats->full_min_multiplicity[full_multiplicity]++;
          if (profile_multiplicity > 1)
            tie_stats->profile_tie_windows++;

          if (full_multiplicity > 1)
            { uint64_t first_code = canonical_codes[min_pos];
              int distinct_code_seen = 0;

              tie_stats->full_tie_windows++;
              tie_stats->full_tied_occurrences +=
                  (uint64_t)full_multiplicity;
              sequence_has_tie = 1;
              if (!previous_window_tied)
                tie_stats->full_tie_runs++;

              for (q = dq_head + 1; q < dq_head + full_multiplicity; q++)
                if (canonical_codes[dq[q]] != first_code)
                  { distinct_code_seen = 1;
                    break;
                  }
              if (distinct_code_seen)
                tie_stats->full_tie_mod_collision_windows++;
              else
                tie_stats->full_tie_repeat_only_windows++;
              previous_window_tied = 1;
            }
          else
            previous_window_tied = 0;
#endif
#if MINIMAX_ROBUST_WINNOW
          int win_start = i - WSIZE + 1;
          if (previous_min >= win_start &&
              position_tuple_equal(previous_min, min_pos,
                                   prof_genome, prof_reads,
                                   hvals, full_hvals))
            min_pos = previous_min;
          previous_min = min_pos;
#endif
          nw[min_pos]++;
        }
#if MINIMAX_ROBUST_WINNOW
      else if (i >= WSIZE - 1)
        previous_min = -1;
#endif
#if MINIMAX_TIE_AUDIT
      else if (i >= WSIZE - 1)
        { tie_stats->empty_windows++;
          previous_window_tied = 0;
        }
#endif
#ifdef DEBUG_HASH
      if (i >= WSIZE - 1 && dq_head < dq_tail)
      { int win_start = i - WSIZE + 1;
        int min_pos   = dq[dq_head];
#if MINIMAX_ROBUST_WINNOW
        if (previous_min >= win_start &&
            position_tuple_equal(previous_min, min_pos,
                                 prof_genome, prof_reads,
                                 hvals, full_hvals))
          min_pos = previous_min;
#endif
        fprintf(stderr, "WINDOW read=%d [%d..%d]  minimizer=@%d g=%u r=%u h=%lu\n",
                read_id, win_start, i, min_pos,
                (unsigned)prof_genome[min_pos],
                (unsigned)prof_reads[min_pos], hvals[min_pos]);
        for (int j = win_start; j <= i; j++) {
          const char *wkp = seq - kmer + 1 + j;
          char wkstr[65];
          for (int b2 = 0; b2 < kmer; b2++) wkstr[b2] = wkp[b2];
          wkstr[kmer] = '\0';
          fprintf(stderr, "  j=%-6d g=%-5u r=%-5u h=%-4lu  %s%s\n",
                  j,
                  (unsigned)prof_genome[j],
                  (unsigned)prof_reads[j], hvals[j],
                  wkstr, (j == min_pos) ? "  <-- MIN" : "");
        }
      }
#endif
    }
#if MINIMAX_TIE_AUDIT
  if (sequence_has_tie)
    tie_stats->sequences_with_full_tie++;
#endif
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
               KposVector    *out
#if MINIMAX_TIE_AUDIT
               , TieStats    *tie_stats
#endif
               )
{
  if (plen <= 0) return 0;

  // Step 1: compute all canonical hash values and derive the configured key.
  uint64_t *hvals = (uint64_t *) malloc((size_t)plen * sizeof(uint64_t));
  uint64_t *full_hvals = NULL;
#if MINIMAX_TIE_AUDIT
  uint64_t *canonical_codes =
      (uint64_t *)malloc((size_t)plen * sizeof(uint64_t));
  int *profile_dq = (int *)malloc((size_t)plen * sizeof(int));
#endif
  if (!hvals) { perror("malloc hvals"); exit(1); }
#if MINIMAX_TIE_AUDIT
  if (!canonical_codes || !profile_dq)
    { perror("malloc tie audit arrays");
      exit(1);
    }
  fill_canonical_codes(seq, plen, K, canonical_codes);
#endif
  { uint64_t fhVal = 0, rhVal = 0;
    hvals[0] = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);
    for (int i = 1; i < plen; i++)
      hvals[i] = NTC64_c((unsigned char)seq[i - K], (unsigned char)seq[i],
                         K, &fhVal, &rhVal);
#if MINIMAX_HASH_REFINE
    full_hvals = (uint64_t *)malloc((size_t)plen * sizeof(uint64_t));
    if (!full_hvals) { perror("malloc full_hvals"); exit(1); }
    memcpy(full_hvals, hvals, (size_t)plen * sizeof(uint64_t));
#endif
#if MINIMAX_HASH_MOD > 0
    // Historical builds use MOD=2003.  Larger-modulus and refinement builds
    // are used only by the controlled chr6 ablation.
    for (int i = 0; i < plen; i++)
      hvals[i] %= MOD;
#endif
  }

  // Step 2: sliding-window minimum with a monotone deque.
  //   The deque stores indices in increasing position order.
  //   Front = index with the minimum (count, hash) in the current window.
  //   We pop from the back any index j whose tuple >= the incoming tuple
  //   (it can never win a future window while the new element is around).
  int  *dq      = (int  *) malloc((size_t)plen * sizeof(int));
  int  *nw      = (int  *) calloc((size_t)plen,  sizeof(int));   // # of windows for which position i is minimizer
  if (!dq || !nw) { perror("malloc deque/nw"); exit(1); }

  count_window_minima(seq, plen, K, prof_genome, prof_reads,
                      hvals, full_hvals, read_id, dq, nw
#if MINIMAX_TIE_AUDIT
                      , canonical_codes, profile_dq, tie_stats
#endif
                      );

  // Step 3: Uniformly sparse thinning.
  // Sort candidates by descending h_nw = (NW, -h_minimax) priority, then
  // greedily accept only those whose position is ≥ RSPACE from every
  // already-accepted seed.
  size_t seed_nwords = ((size_t)plen + 63) >> 6;
  uint64_t *selected = (uint64_t *) calloc(seed_nwords, sizeof(uint64_t));
  uint64_t *blocked  = (uint64_t *) calloc(seed_nwords, sizeof(uint64_t));
  if (!selected || !blocked) { perror("calloc seed bitsets"); exit(1); }

  { int ncand = 0;
    for (int i = 0; i < plen; i++)
      if (nw[i] > 0) ncand++;

    if (ncand > 0)
      { Candidate *cands = (Candidate *) malloc((size_t)ncand * sizeof(Candidate));
        if (!cands) { perror("malloc cands"); exit(1); }
        { int ci = 0;
          for (int i = 0; i < plen; i++)
            if (nw[i] > 0)
              { cands[ci].pos     = i;
                cands[ci].nw_val  = nw[i];
                cands[ci].g_count = prof_genome[i];
                cands[ci].r_count = prof_reads[i];
                cands[ci].hash    = hvals[i];
                ci++;
              }
        }
        if (ncand > 1)
          {
#if MINIMAX_HASH_REFINE
            candidate_full_hvals = full_hvals;
#endif
            qsort(cands, (size_t)ncand, sizeof(Candidate), cand_cmp);
#if MINIMAX_HASH_REFINE
            candidate_full_hvals = NULL;
#endif
          }

        for (int ci = 0; ci < ncand; ci++)
          { size_t p = (size_t)cands[ci].pos;
            if (seed_bit_is_set(blocked, p))
              continue;

            seed_bit_set(selected, p);

            // The old test rejects distances strictly less than RSPACE, so an
            // accepted seed blocks radius RSPACE-1.  RSPACE==0 blocks nothing.
            if (RSPACE > 0)
              { size_t radius = (size_t)RSPACE - 1;
                size_t lo = radius > p ? 0 : p - radius;
                size_t remaining = (size_t)plen - 1 - p;
                size_t hi = radius > remaining ? (size_t)plen - 1 : p + radius;
                seed_bit_set_range(blocked, lo, hi);
              }
          }
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
              wkstr, nw[i] ? "*" : "", seed_bit_is_set(selected, (size_t)i) ? "****" : "");
    }
#endif

  // Step 4: append selected seeds to this block's KPOS run.
  int written = 0;
  for (int i = 0; i < plen; i++)
    { if (!seed_bit_is_set(selected, (size_t)i)) continue;
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

  free(blocked);
  free(selected);
  free(nw);
  free(dq);
#if MINIMAX_TIE_AUDIT
  free(profile_dq);
  free(canonical_codes);
#endif
  free(full_hvals);
  free(hvals);
  return written;
}

typedef struct
  { KposEntry *entry;
    size_t     len;
    int        ready;
#if MINIMAX_TIE_AUDIT
    TieStats   tie_stats;
#endif
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
#if MINIMAX_TIE_AUDIT
      TieStats block_tie_stats = {0};
#endif

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
                       id - first, &vector
#if MINIMAX_TIE_AUDIT
                       , &block_tie_stats
#endif
                       );
        }

      if (vector.len > 1)
        qsort(vector.entry, vector.len, sizeof(KposEntry), cmp_kpos);

      result = queue->result + block;
      pthread_mutex_lock(&queue->lock);
      result->entry = vector.entry;
      result->len   = vector.len;
#if MINIMAX_TIE_AUDIT
      result->tie_stats = block_tie_stats;
#endif
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

#if MINIMAX_TIE_AUDIT
static void write_tie_audit(const char *output_path, const TieStats *stats)
{
  char *stats_path = add_suffix(output_path, ".tie.tsv");
  FILE *out = fopen(stats_path, "w");
  int i;

  if (out == NULL)
    { perror(stats_path);
      exit(1);
    }

#define WRITE_TIE_METRIC(name) \
  fprintf(out, #name "\t%llu\n", \
          (unsigned long long)stats->name)

  fprintf(out, "metric\tvalue\n");
  WRITE_TIE_METRIC(sequences);
  WRITE_TIE_METRIC(sequences_with_complete_window);
  WRITE_TIE_METRIC(sequences_with_full_tie);
  WRITE_TIE_METRIC(complete_windows);
  WRITE_TIE_METRIC(eligible_windows);
  WRITE_TIE_METRIC(empty_windows);
  WRITE_TIE_METRIC(profile_tie_windows);
  WRITE_TIE_METRIC(full_tie_windows);
  WRITE_TIE_METRIC(full_tie_repeat_only_windows);
  WRITE_TIE_METRIC(full_tie_mod_collision_windows);
  WRITE_TIE_METRIC(full_tie_runs);
  WRITE_TIE_METRIC(full_tied_occurrences);
#undef WRITE_TIE_METRIC

  for (i = 1; i <= WSIZE; i++)
    { fprintf(out, "profile_min_multiplicity_%d\t%llu\n", i,
              (unsigned long long)stats->profile_min_multiplicity[i]);
      fprintf(out, "full_min_multiplicity_%d\t%llu\n", i,
              (unsigned long long)stats->full_min_multiplicity[i]);
    }

  if (fclose(out) != 0)
    { perror(stats_path);
      exit(1);
    }

  fprintf(stderr,
          "Tie audit: windows=%llu eligible=%llu profile_ties=%llu "
          "full_ties=%llu repeat_only=%llu modulo_collision=%llu "
          "runs=%llu sequences_with_tie=%llu (%s)\n",
          (unsigned long long)stats->complete_windows,
          (unsigned long long)stats->eligible_windows,
          (unsigned long long)stats->profile_tie_windows,
          (unsigned long long)stats->full_tie_windows,
          (unsigned long long)stats->full_tie_repeat_only_windows,
          (unsigned long long)stats->full_tie_mod_collision_windows,
          (unsigned long long)stats->full_tie_runs,
          (unsigned long long)stats->sequences_with_full_tie,
          stats_path);
  free(stats_path);
}
#endif

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
#if MINIMAX_TIE_AUDIT
  if (WSIZE > TIE_AUDIT_MAX_W)
    { fprintf(stderr, "Tie audit supports -w up to %d, observed %d\n",
              TIE_AUDIT_MAX_W, WSIZE);
      return 1;
    }
#endif
#if MINIMAX_TIE_AUDIT
  fprintf(stderr,
          "Seed selector: hash_mod=%llu hash_refine=%d robust_winnow=%d "
          "tie_audit=1\n",
          (unsigned long long)MOD, MINIMAX_HASH_REFINE,
          MINIMAX_ROBUST_WINNOW);
#else
  fprintf(stderr,
          "Seed selector: hash_mod=%llu hash_refine=%d robust_winnow=%d\n",
          (unsigned long long)MOD, MINIMAX_HASH_REFINE,
          MINIMAX_ROBUST_WINNOW);
#endif
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
#if MINIMAX_TIE_AUDIT
    TieStats total_tie_stats = {0};
#endif
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
#if MINIMAX_TIE_AUDIT
        tie_stats_add(&total_tie_stats, &block->tie_stats);
#endif
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
#if MINIMAX_TIE_AUDIT
    write_tie_audit(output_path, &total_tie_stats);
#endif

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
