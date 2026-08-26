// Development-only copy-number-oriented minimizer variant.
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <DB.h>
#include "nthash.h"
#include "libfastk.h"
#include "copynum.h"
#include "minimaxmer.h"

// #define DEBUG_HASH
// #define DEBUG_NW

int K = 32;          // k-mer length; maximum 32   [-k]
int WSIZE  = 100;    // sliding window width for (count, hash) minimizer   [-w]
int RSPACE = 100;    // minimum spacing between uniformly sparse seeds     [-s]
int EXCLUDE_PERIODIC_BP = -1; // -1 disables long p=1..3 tract exclusion

// Build-time selector controls.  Production keeps the historical h%2003 key,
// but uses D-oriented robust winnowing so that tied window minima are selected
// consistently from a sequence and its reverse complement.
//
// MINIMAX_HASH_MOD=0 uses the full canonical 64-bit ntHash.  With
// MINIMAX_HASH_REFINE=1, the key is (hash % MINIMAX_HASH_MOD, full_hash);
// this preserves every non-tied legacy ordering and only resolves modulo ties.
// MINIMAX_ROBUST_WINNOW=1 implements Schleimer/Winnowmap robust winnowing:
// retain the previous winner while it remains a minimum, otherwise choose the
// extreme minimum in the current traversal direction.
// MINIMAX_CANONICAL_WINNOW=1 orients a tied-window run with
// D(W)=#G+#T-#A-#C over its K+W-1 spelled bases.  Negative-D runs are traversed
// left-to-right with rightmost fallback; positive-D runs are traversed
// right-to-left with leftmost fallback.  State is reset when the orientation
// changes.  This makes the window assignments reverse-complement equivariant.
// MINIMAX_CANONICAL_THINNING controls the strand-safe positional-pair policy
// after window counting.  It follows canonical winnowing by default and is
// independently switchable only for controlled performance attribution.
#ifndef MINIMAX_HASH_MOD
#define MINIMAX_HASH_MOD 2003
#endif
#ifndef MINIMAX_HASH_REFINE
#define MINIMAX_HASH_REFINE 0
#endif
#ifndef MINIMAX_ROBUST_WINNOW
#define MINIMAX_ROBUST_WINNOW 1
#endif
#ifndef MINIMAX_CANONICAL_WINNOW
#define MINIMAX_CANONICAL_WINNOW 1
#endif
#ifndef MINIMAX_CANONICAL_THINNING
#define MINIMAX_CANONICAL_THINNING MINIMAX_CANONICAL_WINNOW
#endif
#ifndef MINIMAX_TIE_AUDIT
#define MINIMAX_TIE_AUDIT 0
#endif
#if (MINIMAX_HASH_REFINE != 0 && MINIMAX_HASH_REFINE != 1) || \
    (MINIMAX_ROBUST_WINNOW != 0 && MINIMAX_ROBUST_WINNOW != 1) || \
    (MINIMAX_CANONICAL_WINNOW != 0 && MINIMAX_CANONICAL_WINNOW != 1) || \
    (MINIMAX_CANONICAL_THINNING != 0 && \
     MINIMAX_CANONICAL_THINNING != 1) || \
    (MINIMAX_TIE_AUDIT != 0 && MINIMAX_TIE_AUDIT != 1)
#error "minimaxmer selector controls must be either 0 or 1"
#endif
#if MINIMAX_HASH_REFINE && MINIMAX_HASH_MOD == 0
#error "MINIMAX_HASH_REFINE requires a non-zero MINIMAX_HASH_MOD"
#endif
#if MINIMAX_CANONICAL_WINNOW && !MINIMAX_ROBUST_WINNOW
#error "MINIMAX_CANONICAL_WINNOW requires MINIMAX_ROBUST_WINNOW"
#endif
#if MINIMAX_TIE_AUDIT && \
    (MINIMAX_HASH_MOD != 2003 || MINIMAX_HASH_REFINE || \
     MINIMAX_ROBUST_WINNOW || MINIMAX_CANONICAL_WINNOW)
#error "MINIMAX_TIE_AUDIT is defined only for the legacy h%2003 selector"
#endif

#if MINIMAX_CANONICAL_WINNOW && !MINIMAX_HASH_REFINE && \
    MINIMAX_HASH_MOD > 0 && MINIMAX_HASH_MOD <= 65536
#define MINIMAX_PACK_POSITION_TUPLE 1
#else
#define MINIMAX_PACK_POSITION_TUPLE 0
#endif

static const uint64_t MOD = MINIMAX_HASH_MOD;

#if MINIMAX_PACK_POSITION_TUPLE
typedef uint16_t PositionHash;
#else
typedef uint64_t PositionHash;
#endif

// ASCII base → 2-bit value  (A=0, C=1, G=2, T=3), matching the DAZZ_DB
// numeric encoding used by himapper's tuple_thread.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

#if MINIMAX_CANONICAL_WINNOW
// D(W) contribution and validity tables.  The production path receives
// uppercase A/C/G/T from Load_Read(..., 2); lowercase entries keep the core
// helper safe for direct unit tests.
static const int8_t base_balance[256] = {
    ['A'] = -1, ['C'] = -1, ['G'] = 1, ['T'] = 1,
    ['a'] = -1, ['c'] = -1, ['g'] = 1, ['t'] = 1,
};
static const uint8_t base_is_acgt[256] = {
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1,
    ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1,
};
#endif

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

typedef Minimaxmer_Kpos KposEntry;

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

typedef struct
  { uint64_t periodic_tracts[3];
    uint64_t periodic_tract_bases[3];
    uint64_t periodic_reads;
    uint64_t periodic_seeds_removed;
  } Periodic_Stats;

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

#if MINIMAX_PACK_POSITION_TUPLE
// Production h%2003 tuples fit in one uint64:
//   [ genome count:16 ][ unused:16 ][ hash key:16 ]
// The full-length array stores only the 16-bit hash.  The canonical deque
// caches this packed key for its O(W) live entries, so its hot comparisons
// remain one integer comparison without paying 8 bytes per sequence position.
static inline uint64_t pack_position_tuple(uint16_t genome_count,
                                           uint64_t hash)
{
  return ((uint64_t)genome_count << 32) |
         (hash & UINT64_C(0xffff));
}

static inline uint64_t position_tuple_key(
    int pos,
    const uint16_t *prof_genome,
    const PositionHash *hvals)
{
  return pack_position_tuple(prof_genome[pos], (uint64_t)hvals[pos]);
}
#endif

// Compare (ca_g, ha) < (cb_g, hb):
//   primary:   lower genome count  -> more unique in genome -> better seed
//   secondary: lower hash key      -> stable tiebreaker
// In refinement builds only, the full hash is the final collision resolver.
static inline int position_tuple_lt(int a, int b,
                                    const uint16_t *prof_genome,
                                    const PositionHash *hvals,
                                    const uint64_t *full_hvals)
#if MINIMAX_PACK_POSITION_TUPLE
{
  (void)full_hvals;
  return position_tuple_key(a, prof_genome, hvals) <
         position_tuple_key(b, prof_genome, hvals);
}
#else
{
  if (prof_genome[a] != prof_genome[b])
    return prof_genome[a] < prof_genome[b];
  if (hvals[a] != hvals[b])
    return hvals[a] < hvals[b];
#if MINIMAX_HASH_REFINE
  return full_hvals[a] < full_hvals[b];
#else
  (void)full_hvals;
  return 0;
#endif
}
#endif

static inline int position_tuple_equal(int a, int b,
                                       const uint16_t *prof_genome,
                                       const PositionHash *hvals,
                                       const uint64_t *full_hvals)
#if MINIMAX_PACK_POSITION_TUPLE
{
  (void)full_hvals;
  return position_tuple_key(a, prof_genome, hvals) ==
         position_tuple_key(b, prof_genome, hvals);
}
#else
{
  if (prof_genome[a] != prof_genome[b] ||
      hvals[a] != hvals[b])
    return 0;
#if MINIMAX_HASH_REFINE
  return full_hvals[a] == full_hvals[b];
#else
  (void)full_hvals;
  return 1;
#endif
}
#endif

#if MINIMAX_TIE_AUDIT
static inline int profile_tuple_lt(int a, int b,
                                   const uint16_t *prof_genome)
{
  return prof_genome[a] < prof_genome[b];
}

static inline int profile_tuple_equal(int a, int b,
                                      const uint16_t *prof_genome)
{
  return prof_genome[a] == prof_genome[b];
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
#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
  uint64_t tuple;
#else
  uint16_t g_count;
  uint64_t hash;
#endif
} Candidate;

// Stable radix sorting needs one auxiliary Candidate array.  Production
// workers retain this grow-only buffer across reads so the fast path does not
// add a malloc/free pair to every read; direct find_seeds() callers get the
// same behavior through a call-local scratch object below.
typedef struct {
  Candidate *temporary;
  size_t     capacity;
} CandidateRadixScratch;

static void init_candidate_radix_scratch(CandidateRadixScratch *scratch)
{
  scratch->temporary = NULL;
  scratch->capacity = 0;
}

static void free_candidate_radix_scratch(CandidateRadixScratch *scratch)
{
  free(scratch->temporary);
  scratch->temporary = NULL;
  scratch->capacity = 0;
}

#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
typedef char Candidate_must_be_16_bytes[
    (sizeof(Candidate) == 16) ? 1 : -1];
#endif

#if MINIMAX_HASH_REFINE
// Keep Candidate byte-layout identical to the historical build.  A larger
// element size can change qsort's ordering inside comparator-equal groups and
// would confound collision refinement with an unrelated spacing change.
static _Thread_local const uint64_t *candidate_full_hvals;
#endif

// Compare candidates by descending h_nw priority:
//   1. nw      descending (more windows → higher priority)
//   2. g_count ascending  (rarer in genome → better seed)
//   3. hash    ascending  (stable tiebreaker)
static int cand_cmp(const void *a, const void *b)
{ const Candidate *ca = (const Candidate *)a;
  const Candidate *cb = (const Candidate *)b;
  if (ca->nw_val != cb->nw_val)
    return (cb->nw_val > ca->nw_val) ? 1 : -1;
#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
  if (ca->tuple != cb->tuple)
    return (ca->tuple > cb->tuple) ? 1 : -1;
#else
  if (ca->g_count != cb->g_count)
    return (ca->g_count > cb->g_count) ? 1 : -1;
  if (ca->hash != cb->hash)
    return (ca->hash > cb->hash) ? 1 : -1;
#if MINIMAX_HASH_REFINE
  { uint64_t fa = candidate_full_hvals[ca->pos];
    uint64_t fb = candidate_full_hvals[cb->pos];
    if (fa != fb)
      return (fa > fb) ? 1 : -1;
  }
#endif
#endif
  return 0;
}

#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
static Candidate *reserve_candidate_radix_scratch(
    CandidateRadixScratch *scratch, size_t need)
{
  Candidate *temporary;
  size_t capacity;

  if (need <= scratch->capacity)
    return scratch->temporary;
  if (need > SIZE_MAX / sizeof(Candidate))
    { fprintf(stderr, "candidate radix capacity overflow\n");
      exit(1);
    }
  capacity = scratch->capacity == 0 ? 256 : scratch->capacity;
  while (capacity < need)
    { if (capacity > SIZE_MAX / 2)
        { capacity = need;
          break;
        }
      capacity *= 2;
    }
  temporary = (Candidate *)realloc(
      scratch->temporary, capacity * sizeof(Candidate));
  if (temporary == NULL)
    { perror("realloc candidate radix scratch");
      exit(1);
    }
  scratch->temporary = temporary;
  scratch->capacity = capacity;
  return temporary;
}

// Production Candidate::tuple occupies exactly 48 bits.  Six stable LSD
// byte passes therefore order tuple ascending; the final pass maps nw to its
// reverse byte so larger window counts precede smaller ones.  Input candidates
// are built in position order, hence equal keys remain position-stable.  The
// existing group-local position qsort remains in place as the semantic guard.
static void radix_sort_candidates(Candidate *candidates, size_t length,
                                  CandidateRadixScratch *scratch)
{
  Candidate *source = candidates;
  Candidate *target =
      reserve_candidate_radix_scratch(scratch, length);
  unsigned int pass;

  for (pass = 0; pass < 7; pass++)
    { size_t offset[256] = {0};
      size_t total = 0;
      size_t i;

      if (pass < 6)
        { unsigned int shift = pass * 8;
          for (i = 0; i < length; i++)
            offset[(source[i].tuple >> shift) & UINT64_C(0xff)]++;
        }
      else
        for (i = 0; i < length; i++)
          offset[UINT8_MAX - (unsigned int)source[i].nw_val]++;

      for (i = 0; i < 256; i++)
        { size_t count = offset[i];
          offset[i] = total;
          total += count;
        }

      if (pass < 6)
        { unsigned int shift = pass * 8;
          for (i = 0; i < length; i++)
            { unsigned int key = (unsigned int)
                  ((source[i].tuple >> shift) & UINT64_C(0xff));
              target[offset[key]++] = source[i];
            }
        }
      else
        for (i = 0; i < length; i++)
          { unsigned int key =
                UINT8_MAX - (unsigned int)source[i].nw_val;
            target[offset[key]++] = source[i];
          }

      { Candidate *swap = source;
        source = target;
        target = swap;
      }
    }

  if (source != candidates)
    memcpy(candidates, source, length * sizeof(Candidate));
}
#endif

#if MINIMAX_CANONICAL_THINNING
static int cand_pos_cmp(const void *a, const void *b)
{
  const Candidate *ca = (const Candidate *)a;
  const Candidate *cb = (const Candidate *)b;
  return (ca->pos > cb->pos) - (ca->pos < cb->pos);
}

static int candidate_priority_equal(const Candidate *a, const Candidate *b,
                                    const uint64_t *full_hvals)
{
  if (a->nw_val != b->nw_val)
    return 0;
#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
  (void)full_hvals;
  return a->tuple == b->tuple;
#else
  if (
      a->g_count != b->g_count ||
      a->hash != b->hash)
    return 0;
#if MINIMAX_HASH_REFINE
  return full_hvals[a->pos] == full_hvals[b->pos];
#else
  (void)full_hvals;
  return 1;
#endif
#endif
}
#endif

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
static inline __attribute__((always_inline))
void seed_bit_set_range(uint64_t *bits, size_t lo, size_t hi)
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

// Clear the inclusive bit range [lo, hi].  This mirrors
// seed_bit_set_range(), including its no-shift-by-64 boundary handling.
static inline __attribute__((always_inline))
uint64_t seed_bit_clear_range(uint64_t *bits, size_t lo, size_t hi)
{
  size_t first = lo >> 6;
  size_t last  = hi >> 6;
  size_t word;
  unsigned int lbit = (unsigned int)(lo & 63);
  unsigned int hbit = (unsigned int)(hi & 63);
  uint64_t left_mask  = UINT64_MAX << lbit;
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
// tract.  A p-periodic tract is a maximal run of bases[i]==bases[i-p].
// Its ClassPro-compatible length is the number of complete motif bases,
// p*floor(raw_span/p); the option excludes it only when that value is
// strictly greater than EXCLUDE_PERIODIC_BP.  Period-2 and period-3 tracts
// that are homopolymers are left to the period-1 pass so the diagnostics are
// disjoint.  The selected seeds that remain are an exact subset of the
// baseline thinning result.
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

static inline __attribute__((always_inline))
void accept_sparse_seed(uint64_t *selected, uint64_t *blocked,
                        size_t plen, size_t pos)
{
  seed_bit_set(selected, pos);
  // The historical test rejects distances strictly less than RSPACE, so an
  // accepted seed blocks radius RSPACE-1.  RSPACE==0 blocks nothing.
  if (RSPACE > 0)
    { size_t radius = (size_t)RSPACE - 1;
      size_t lo = radius > pos ? 0 : pos - radius;
      size_t remaining = plen - 1 - pos;
      size_t hi = radius > remaining ? plen - 1 : pos + radius;
      seed_bit_set_range(blocked, lo, hi);
    }
}

#if MINIMAX_CANONICAL_WINNOW
// Return the D(W) balance for bases[start..start+span), or -1 if a direct
// caller supplied a non-ACGT character.  DAZZ Load_Read(..., 2) itself only
// yields uppercase A/C/G/T.
static int initial_window_balance(const char *bases, int start, int span,
                                  int *balance)
{
  int d = 0;
  int j;

  for (j = 0; j < span; j++)
    { unsigned char base = (unsigned char)bases[start + j];
      if (!base_is_acgt[base])
        return -1;
      d += base_balance[base];
    }
  *balance = d;
  return 0;
}

static int roll_window_balance(const char *bases, int start, int span,
                               int *balance)
{
  unsigned char outgoing = (unsigned char)bases[start - 1];
  unsigned char incoming = (unsigned char)bases[start + span - 1];

  if (!base_is_acgt[outgoing] || !base_is_acgt[incoming])
    return -1;
  *balance += base_balance[incoming] - base_balance[outgoing];
  return 0;
}

// D is non-zero whenever the spelled window length is odd.  For even spans,
// an exact comparison against the reverse complement supplies an
// anti-symmetric secondary orientation.  Positive selects the right-to-left
// traversal; zero means that the window is exactly
// self-reverse-complementary.
static int exact_window_orientation(const char *bases, int start, int span,
                                    int balance)
{
  int j;

  if (balance > 0)
    return 1;
  if (balance < 0)
    return -1;

  for (j = 0; j < span; j++)
    { unsigned char left = (unsigned char)bases[start + j];
      unsigned char right =
          (unsigned char)bases[start + span - 1 - j];
      uint8_t forward;
      uint8_t reverse;

      if (!base_is_acgt[left] || !base_is_acgt[right])
        return 2;
      forward = base2bit[left];
      reverse = (uint8_t)(UINT8_C(3) - base2bit[right]);
      if (forward < reverse)
        return 1;
      if (forward > reverse)
        return -1;
    }
  return 0;
}

// Refresh the cached end of the deque's equal-minimum prefix.  The cache is
// invalidated whenever a tail contraction removes its slot, or the front moves
// past it.  New equal minima extend it at insertion time below.  Consequently
// tuple equality is loaded only when a different minimum group reaches the
// front, rather than once per complete window.
typedef struct
{
  int pos;
#if MINIMAX_PACK_POSITION_TUPLE
  uint64_t key;
#endif
} CanonicalMinimum;

typedef struct
{
  void *value;
  int   narrow;
} WindowCounts;

static inline CanonicalMinimum invalid_canonical_minimum(void)
{
  CanonicalMinimum result;

  result.pos = -1;
#if MINIMAX_PACK_POSITION_TUPLE
  result.key = 0;
#endif
  return result;
}

static inline CanonicalMinimum canonical_minimum_at(
    int pos,
    const uint16_t *prof_genome,
    const PositionHash *hvals)
{
  CanonicalMinimum result;

  result.pos = pos;
#if MINIMAX_PACK_POSITION_TUPLE
  result.key = position_tuple_key(pos, prof_genome, hvals);
#else
  (void)prof_genome;
  (void)hvals;
#endif
  return result;
}

static inline int canonical_minimum_lt(
    CanonicalMinimum a,
    CanonicalMinimum b,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals)
{
#if MINIMAX_PACK_POSITION_TUPLE
  (void)prof_genome;
  (void)hvals;
  (void)full_hvals;
  return a.key < b.key;
#else
  return position_tuple_lt(
      a.pos, b.pos, prof_genome, hvals, full_hvals);
#endif
}

static inline int canonical_minimum_equal(
    CanonicalMinimum a,
    CanonicalMinimum b,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals)
{
#if MINIMAX_PACK_POSITION_TUPLE
  (void)prof_genome;
  (void)hvals;
  (void)full_hvals;
  return a.key == b.key;
#else
  return position_tuple_equal(
      a.pos, b.pos, prof_genome, hvals, full_hvals);
#endif
}

static void init_window_counts(WindowCounts *counts, int plen)
{
  size_t element_size;

  counts->narrow = WSIZE <= UINT8_MAX;
  element_size = counts->narrow ? sizeof(uint8_t) : sizeof(int);
  if ((size_t)plen > SIZE_MAX / element_size)
    { fprintf(stderr, "canonical winnowing: window-count size overflow\n");
      exit(1);
    }
  counts->value = calloc((size_t)plen, element_size);
  if (counts->value == NULL && plen > 0)
    { perror("calloc window counts");
      exit(1);
    }
}

static void free_window_counts(WindowCounts *counts)
{
  free(counts->value);
  counts->value = NULL;
  counts->narrow = 0;
}

static inline __attribute__((always_inline))
void increment_window_count(void *value, int narrow, int pos)
{
  if (__builtin_expect(narrow, 1))
    ((uint8_t *)value)[pos]++;
  else
    ((int *)value)[pos]++;
}

static inline unsigned int get_window_count(
    const WindowCounts *counts, int pos)
{
  if (__builtin_expect(counts->narrow, 1))
    return ((const uint8_t *)counts->value)[pos];
  return (unsigned int)((const int *)counts->value)[pos];
}

static CanonicalMinimum cached_rightmost_equal_minimum(
    const CanonicalMinimum *dq,
    size_t mask,
    size_t head,
    size_t tail,
    size_t *equal_end,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals)
{
  size_t q = *equal_end;

  if (q == SIZE_MAX || q < head || q >= tail)
    {
      q = head;
      while (q + 1 < tail &&
             canonical_minimum_equal(
                 dq[(q + 1) & mask], dq[head & mask],
                 prof_genome, hvals, full_hvals))
        q++;
    }
  *equal_end = q;
  return dq[q & mask];
}

typedef struct
{
  CanonicalMinimum *leftmost;
  size_t             capacity;
  size_t             mask;
  CanonicalMinimum *deque;
  size_t             deque_capacity;
  size_t             deque_mask;
} CanonicalWindowScratch;

static size_t canonical_ring_capacity(size_t need)
{
  size_t capacity = 1;

  while (capacity < need)
    {
      if (capacity > SIZE_MAX / 2)
        { fprintf(stderr, "canonical winnowing: ring size overflow\n");
          exit(1);
        }
      capacity <<= 1;
    }
  if (capacity > SIZE_MAX / sizeof(CanonicalMinimum))
    { fprintf(stderr, "canonical winnowing: ring allocation overflow\n");
      exit(1);
    }
  return capacity;
}

static void init_canonical_window_scratch(CanonicalWindowScratch *scratch)
{
  size_t deque_capacity =
      canonical_ring_capacity((size_t)WSIZE + 1);
  size_t capacity = deque_capacity;

  // Most positive-D subruns are short.  A small floor removes their replay
  // without making scratch scale with the sequence length; WSIZE+1 still
  // guarantees enough space for the exact replay fallback.
  if (capacity < 512)
    capacity = 512;
  if (capacity > SIZE_MAX / sizeof(CanonicalMinimum))
    { fprintf(stderr, "canonical winnowing: scratch size overflow\n");
      exit(1);
    }
  scratch->leftmost =
      (CanonicalMinimum *)malloc(capacity * sizeof(CanonicalMinimum));
  scratch->deque =
      (CanonicalMinimum *)malloc(
          deque_capacity * sizeof(CanonicalMinimum));
  if (scratch->leftmost == NULL || scratch->deque == NULL)
    { perror("malloc canonical window scratch");
      exit(1);
    }
  scratch->capacity = capacity;
  scratch->mask = capacity - 1;
  scratch->deque_capacity = deque_capacity;
  scratch->deque_mask = deque_capacity - 1;
}

static void free_canonical_window_scratch(CanonicalWindowScratch *scratch)
{
  free(scratch->leftmost);
  free(scratch->deque);
  scratch->leftmost = NULL;
  scratch->deque = NULL;
  scratch->capacity = 0;
  scratch->mask = 0;
  scratch->deque_capacity = 0;
  scratch->deque_mask = 0;
}

// Assign a positive-D tie subrun in its canonical (right-to-left) traversal.
// The circular monotone deque keeps the leftmost representative of equal
// minima.  right_anchor is the unique winner of the immediately following
// window, or has pos=-1 at an orientation/empty/sequence boundary.
static __attribute__((cold, noinline))
void assign_reverse_tie_subrun(int run_start, int run_end,
                               CanonicalMinimum right_anchor,
                               const uint16_t *prof_genome,
                               const PositionHash *hvals,
                               const uint64_t *full_hvals,
                               CanonicalWindowScratch *scratch,
                               void *nw,
                               int nw_narrow)
{
  size_t head = 0;
  size_t tail = 0;
  CanonicalMinimum previous = right_anchor;
  int high = run_end + WSIZE - 1;
  int pos;

  for (pos = high; pos >= run_start; pos--)
    {
      if (prof_genome[pos] != 0)
        {
          CanonicalMinimum incoming =
              canonical_minimum_at(
                  pos, prof_genome, hvals);

          while (head < tail)
            { CanonicalMinimum back =
                  scratch->leftmost[(tail - 1) & scratch->mask];
              if (canonical_minimum_lt(
                      back, incoming, prof_genome, hvals, full_hvals))
                break;
              tail--;
            }
          scratch->leftmost[tail & scratch->mask] = incoming;
          tail++;
        }

      if (pos <= run_end)
        {
          int window_end = pos + WSIZE - 1;
          CanonicalMinimum minimum;
          CanonicalMinimum winner;

          while (head < tail &&
                 scratch->leftmost[head & scratch->mask].pos > window_end)
            head++;
          if (head == tail)
            { fprintf(stderr,
                      "canonical winnowing: empty deferred tie window at %d\n",
                      pos);
              exit(1);
            }

          minimum = scratch->leftmost[head & scratch->mask];
          if (previous.pos >= pos && previous.pos <= window_end &&
              canonical_minimum_equal(
                  previous, minimum, prof_genome, hvals, full_hvals))
            winner = previous;
          else
            winner = minimum;
          increment_window_count(nw, nw_narrow, winner.pos);
          previous = winner;
        }
    }
}

// The forward scan already knows the leftmost minimum of every positive-D
// window.  Evaluating those positions backwards is exactly the right-to-left
// robust recurrence and avoids rebuilding the deque over the run plus its
// WSIZE-1 halo.
static void assign_buffered_reverse_tie_subrun(
    int run_start,
    size_t run_length,
    CanonicalMinimum right_anchor,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals,
    const CanonicalMinimum *leftmost,
    void *nw,
    int nw_narrow)
{
  CanonicalMinimum previous = right_anchor;
  size_t j = run_length;

  while (j > 0)
    {
      int start;
      CanonicalMinimum minimum;
      CanonicalMinimum winner;

      j--;
      start = run_start + (int)j;
      minimum = leftmost[j];
      if (previous.pos >= start && previous.pos < start + WSIZE &&
          canonical_minimum_equal(
              previous, minimum, prof_genome, hvals, full_hvals))
        winner = previous;
      else
        winner = minimum;
      increment_window_count(nw, nw_narrow, winner.pos);
      previous = winner;
    }
}

static void assign_positive_tie_subrun(
    int run_start,
    int run_end,
    size_t run_length,
    CanonicalMinimum right_anchor,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals,
    CanonicalWindowScratch *scratch,
    void *nw,
    int nw_narrow)
{
  if (run_start < 0)
    return;
  if (run_length != (size_t)(run_end - run_start + 1))
    { fprintf(stderr,
              "canonical winnowing: inconsistent positive-D subrun\n");
      exit(1);
    }

  if (run_length <= scratch->capacity)
    assign_buffered_reverse_tie_subrun(
        run_start, run_length, right_anchor,
        prof_genome, hvals, full_hvals,
        scratch->leftmost, nw, nw_narrow);
  else
    assign_reverse_tie_subrun(
        run_start, run_end, right_anchor,
        prof_genome, hvals, full_hvals,
        scratch, nw, nw_narrow);
}

// D-oriented robust winnowing.  Unique minima are direction-independent.
// Consecutive tied windows are divided at D-sign changes:
//   D < 0: traverse left-to-right, retain the actual previous winner, and
//          choose the rightmost minimum when a new winner is needed.
//   D > 0: mirror the same policy right-to-left with leftmost fallback.
// D==0 uses an exact window-vs-RC comparison.  An exactly self-RC tied window
// has no strand-safe singleton winner, so both extreme minima receive the
// window and directional state is reset.
static inline __attribute__((always_inline))
void count_canonical_robust_minima(const char *seq,
                                   int plen,
                                   int kmer,
                                   const uint16_t *prof_genome,
                                   const PositionHash *hvals,
                                   const uint64_t *full_hvals,
                                   int read_id,
                                   void *nw,
                                   int nw_narrow,
                                   CanonicalWindowScratch *scratch)
{
  const char *bases = seq - kmer + 1;
  CanonicalMinimum *dq = scratch->deque;
  size_t dq_mask = scratch->deque_mask;
  int span;
  size_t dq_head = 0;
  size_t dq_tail = 0;
  size_t equal_end = SIZE_MAX;
  int previous_class = 0;  // 0=reset, 1=unique, -1=negative, 2=positive
  CanonicalMinimum forward_previous = invalid_canonical_minimum();
  int positive_start = -1;
  size_t positive_length = 0;
  int previous_window_tied = 0;
  int balance = 0;
  int i;

  if (plen < WSIZE)
    return;
  span = kmer + WSIZE - 1;

  for (i = 0; i < plen; i++)
    {
      int start;

      if (prof_genome[i] != 0)
        {
          CanonicalMinimum incoming =
              canonical_minimum_at(
                  i, prof_genome, hvals);

          // Retain equal tuples.  The front is therefore the leftmost
          // minimum.  Keep the historical boolean comparator in this hot
          // loop; only an incoming item adjacent to the cached minimum prefix
          // needs a separate equality test.
          while (dq_head < dq_tail &&
                 canonical_minimum_lt(
                     incoming, dq[(dq_tail - 1) & dq_mask],
                     prof_genome, hvals, full_hvals))
            dq_tail--;
          // A tail contraction can invalidate the cached equal-prefix end,
          // including the one-element case where its slot is now overwritten.
          if (equal_end != SIZE_MAX && equal_end >= dq_tail)
            equal_end = SIZE_MAX;
          if (equal_end != SIZE_MAX &&
              dq_head < dq_tail &&
              equal_end == dq_tail - 1 &&
              canonical_minimum_equal(
                  incoming, dq[dq_head & dq_mask],
                  prof_genome, hvals, full_hvals))
            equal_end = dq_tail;
          dq[dq_tail & dq_mask] = incoming;
          dq_tail++;
        }

      start = i - WSIZE + 1;
      while (dq_head < dq_tail &&
             dq[dq_head & dq_mask].pos < start)
        dq_head++;
      if (start < 0)
        continue;

      if (dq_head == dq_tail)
        {
          if (positive_start >= 0)
            assign_positive_tie_subrun(
                positive_start, start - 1, positive_length,
                invalid_canonical_minimum(),
                prof_genome, hvals, full_hvals,
                scratch, nw, nw_narrow);
          positive_start = -1;
          positive_length = 0;
          previous_class = 0;
          forward_previous = invalid_canonical_minimum();
          previous_window_tied = 0;
          continue;
        }

      {
        CanonicalMinimum leftmost = dq[dq_head & dq_mask];
        CanonicalMinimum rightmost =
            cached_rightmost_equal_minimum(
                dq, dq_mask, dq_head, dq_tail, &equal_end,
                prof_genome, hvals, full_hvals);
        int tied = equal_end > dq_head;

        if (!tied)
          {
            if (positive_start >= 0)
              {
                assign_positive_tie_subrun(
                    positive_start, start - 1, positive_length, leftmost,
                    prof_genome, hvals, full_hvals,
                    scratch, nw, nw_narrow);
                positive_start = -1;
                positive_length = 0;
              }
            increment_window_count(nw, nw_narrow, leftmost.pos);
            previous_class = 1;
            forward_previous = leftmost;
            previous_window_tied = 0;
            continue;
          }

        if (previous_window_tied)
          {
            if (roll_window_balance(bases, start, span, &balance) != 0)
              { fprintf(stderr,
                        "canonical winnowing: non-ACGT base in read %d "
                        "window %d\n",
                        read_id, start);
                exit(1);
              }
          }
        else if (initial_window_balance(bases, start, span, &balance) != 0)
          {
            fprintf(stderr,
                    "canonical winnowing: non-ACGT base in read %d window %d\n",
                    read_id, start);
            exit(1);
          }
        previous_window_tied = 1;

        {
          // A sum of an odd number of +/-1 terms cannot be zero.  Production
          // K=32,W=100 has span=131, so it never needs the more expensive
          // exact window-versus-RC fallback.
          int orientation =
              (span & 1) ? (balance > 0 ? 1 : -1) :
              exact_window_orientation(bases, start, span, balance);

          if (orientation == 2)
            { fprintf(stderr,
                      "canonical winnowing: non-ACGT base in read %d "
                      "window %d\n",
                      read_id, start);
              exit(1);
            }

          if (orientation < 0)
            {
              CanonicalMinimum winner;

              if (positive_start >= 0)
                assign_positive_tie_subrun(
                    positive_start, start - 1, positive_length,
                    invalid_canonical_minimum(),
                    prof_genome, hvals, full_hvals,
                    scratch, nw, nw_narrow);
              positive_start = -1;
              positive_length = 0;
              if (previous_class != -1 && previous_class != 1)
                forward_previous = invalid_canonical_minimum();

              if (forward_previous.pos >= start &&
                  forward_previous.pos <= i &&
                  canonical_minimum_equal(
                      forward_previous, leftmost,
                      prof_genome, hvals, full_hvals))
                winner = forward_previous;
              else
                winner = rightmost;
              increment_window_count(nw, nw_narrow, winner.pos);
              forward_previous = winner;
              previous_class = -1;
            }
          else if (orientation > 0)
            {
              if (previous_class != 2)
                { positive_start = start;
                  positive_length = 0;
                }
              if (positive_length < scratch->capacity)
                scratch->leftmost[positive_length] = leftmost;
              positive_length++;
              previous_class = 2;
              forward_previous = invalid_canonical_minimum();
            }
          else
            {
              if (positive_start >= 0)
                assign_positive_tie_subrun(
                    positive_start, start - 1, positive_length,
                    invalid_canonical_minimum(),
                    prof_genome, hvals, full_hvals,
                    scratch, nw, nw_narrow);
              positive_start = -1;
              positive_length = 0;
              increment_window_count(nw, nw_narrow, leftmost.pos);
              if (rightmost.pos != leftmost.pos)
                increment_window_count(nw, nw_narrow, rightmost.pos);
              previous_class = 0;
              forward_previous = invalid_canonical_minimum();
            }
        }
      }
    }

  if (positive_start >= 0)
    assign_positive_tie_subrun(
        positive_start, plen - WSIZE, positive_length,
        invalid_canonical_minimum(),
        prof_genome, hvals, full_hvals,
        scratch, nw, nw_narrow);
}
#endif

#if MINIMAX_CANONICAL_WINNOW
static void count_window_minima_with_scratch(
    const char *seq,
    int plen,
    int kmer,
    const uint16_t *prof_genome,
    const PositionHash *hvals,
    const uint64_t *full_hvals,
    int read_id,
    WindowCounts *nw,
    CanonicalWindowScratch *scratch)
{
  // Dispatch once per sequence.  The always-inlined selector sees a constant
  // storage width in each arm, so assigned windows do not pay this test.
  if (__builtin_expect(nw->narrow, 1))
    count_canonical_robust_minima(
        seq, plen, kmer, prof_genome,
        hvals, full_hvals, read_id, nw->value, 1, scratch);
  else
    count_canonical_robust_minima(
        seq, plen, kmer, prof_genome,
        hvals, full_hvals, read_id, nw->value, 0, scratch);
}
#endif

// Count how many complete windows are assigned to each position.  Keeping this
// selector separate from hashing and spacing lets the experiment validate both
// deque policies against a simple O(nW) oracle without changing KPOS output.
static __attribute__((unused))
void count_window_minima(const char *seq,
                                int plen,
                                int kmer,
                                const uint16_t *prof_genome,
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
#if MINIMAX_CANONICAL_WINNOW
  CanonicalWindowScratch scratch;
  WindowCounts compact_nw;
  int i;
#if MINIMAX_PACK_POSITION_TUPLE
  PositionHash *compact_hvals = NULL;
  const PositionHash *selector_hvals = compact_hvals;

  // This internal test/helper entry receives the configured hash key, not an
  // unmodded full hash.  In this arm MINIMAX_HASH_MOD<=65536, so every key is
  // losslessly representable in 16 bits.
  if (plen > 0)
    {
      compact_hvals =
          (PositionHash *)malloc((size_t)plen * sizeof(PositionHash));
      if (compact_hvals == NULL)
        { perror("malloc compact selector hashes");
          exit(1);
        }
      for (i = 0; i < plen; i++)
        {
          if (hvals[i] > UINT16_MAX)
            { fprintf(stderr,
                      "canonical winnowing: selector hash exceeds 16 bits\n");
              exit(1);
            }
          compact_hvals[i] = (PositionHash)hvals[i];
        }
      selector_hvals = compact_hvals;
    }
#else
  const PositionHash *selector_hvals = hvals;
#endif

  (void)dq;
  init_canonical_window_scratch(&scratch);
  init_window_counts(&compact_nw, plen);
  count_window_minima_with_scratch(
      seq, plen, kmer, prof_genome,
      selector_hvals, full_hvals, read_id, &compact_nw, &scratch);
  for (i = 0; i < plen; i++)
    nw[i] += (int)get_window_count(&compact_nw, i);
  free_window_counts(&compact_nw);
  free_canonical_window_scratch(&scratch);
#if MINIMAX_PACK_POSITION_TUPLE
  free(compact_hvals);
#endif
#else
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
                                  prof_genome))
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
                                    prof_genome, hvals, full_hvals))
            dq_tail--;
#else
          // Historical selector: discard only strictly greater back elements,
          // retaining every equal tuple so the deque front is the leftmost.
          while (dq_head < dq_tail &&
                 position_tuple_lt(i, dq[dq_tail-1],
                                   prof_genome, hvals, full_hvals))
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
                                       prof_genome))
                break;
              profile_multiplicity++;
            }
          for (q = dq_head + 1; q < dq_tail; q++)
            { if (!position_tuple_equal(dq[q], min_pos,
                                        prof_genome, hvals, full_hvals))
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
                                   prof_genome, hvals, full_hvals))
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
                                 prof_genome, hvals, full_hvals))
          min_pos = previous_min;
#endif
        fprintf(stderr, "WINDOW read=%d [%d..%d]  minimizer=@%d g=%u h=%lu\n",
                read_id, win_start, i, min_pos,
                (unsigned)prof_genome[min_pos], hvals[min_pos]);
        for (int j = win_start; j <= i; j++) {
          const char *wkp = seq - kmer + 1 + j;
          char wkstr[65];
          for (int b2 = 0; b2 < kmer; b2++) wkstr[b2] = wkp[b2];
          wkstr[kmer] = '\0';
          fprintf(stderr, "  j=%-6d g=%-5u h=%-4lu  %s%s\n",
                  j,
                  (unsigned)prof_genome[j], hvals[j],
                  wkstr, (j == min_pos) ? "  <-- MIN" : "");
        }
      }
#endif
    }
#if MINIMAX_TIE_AUDIT
  if (sequence_has_tie)
    tie_stats->sequences_with_full_tie++;
#endif
#endif
}

// Select seeds using a sliding-window (genome_count, hash) minimizer.
//
// For every window of WSIZE consecutive k-mers the one with the lexicographically
// smallest (prof_genome[i], hVal[i]) tuple is selected.
// The result is the union of all per-window selections.
//
//   seq         – pointer to position (K-1) in the read
//   plen        – number of k-mers (rlen - K + 1)
//   K           – k-mer length
//   prof_genome – FastK uint16 genome k-mer counts, positions 0..plen-1
//   read_id     – 0-based read index
//   out         – block-local KPOS vector
//
// Returns the number of entries appended.  Production workers pass a reusable
// canonical scratch buffer; the public wrapper below preserves the historical
// find_seeds() interface for direct callers and selector tests.
static int find_seeds_with_scratch(
               const char    *seq,
               const int      plen,
               const int      K,
               const uint16  *prof_genome,
               const int      read_id,
               KposVector    *out,
               Periodic_Stats *periodic_stats,
               CandidateRadixScratch *candidate_scratch
#if MINIMAX_CANONICAL_WINNOW
               , CanonicalWindowScratch *window_scratch
#endif
#if MINIMAX_TIE_AUDIT
               , TieStats    *tie_stats
#endif
               )
{
#if !(MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING)
  (void)candidate_scratch;
#endif
  if (plen <= 0) return 0;

  // Step 1: compute all canonical hash values and derive the configured key.
  PositionHash *hvals =
      (PositionHash *)malloc((size_t)plen * sizeof(PositionHash));
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
#if MINIMAX_PACK_POSITION_TUPLE
    uint64_t hash = NTC64_b(seq - K + 1, K, &fhVal, &rhVal) % MOD;
    hvals[0] = (PositionHash)hash;
    for (int i = 1; i < plen; i++)
      { hash = NTC64_c(
                   (unsigned char)seq[i - K],
                   (unsigned char)seq[i], K, &fhVal, &rhVal) % MOD;
        hvals[i] = (PositionHash)hash;
      }
#else
    hvals[0] = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);
    for (int i = 1; i < plen; i++)
      hvals[i] = NTC64_c(
          (unsigned char)seq[i - K], (unsigned char)seq[i],
          K, &fhVal, &rhVal);
#if MINIMAX_HASH_REFINE
    full_hvals = (uint64_t *)malloc((size_t)plen * sizeof(uint64_t));
    if (!full_hvals) { perror("malloc full_hvals"); exit(1); }
    memcpy(full_hvals, hvals, (size_t)plen * sizeof(PositionHash));
#endif
#if MINIMAX_HASH_MOD > 0
    // Historical builds use MOD=2003.  Larger-modulus and refinement builds
    // are used only by the controlled chr6 ablation.
    for (int i = 0; i < plen; i++)
      hvals[i] %= MOD;
#endif
#endif
  }
  // Step 2: sliding-window minimum with a monotone deque.
  //   The deque stores indices in increasing position order.
  //   Front = index with the minimum (count, hash) in the current window.
  //   Canonical mode retains the equal-minimum prefix so D-oriented robust
  //   traversal can choose either extreme; historical modes use their
  //   original strict/non-strict pop policy.
#if MINIMAX_CANONICAL_WINNOW
  WindowCounts nw;

  init_window_counts(&nw, plen);
  count_window_minima_with_scratch(
      seq, plen, K, prof_genome,
      hvals, full_hvals, read_id, &nw, window_scratch);
#define NW_AT(pos) get_window_count(&nw, (pos))
#else
  int *dq = (int *)malloc((size_t)plen * sizeof(int));
  int *nw = (int *)calloc((size_t)plen, sizeof(int));

  if (!dq || !nw) { perror("malloc deque/nw"); exit(1); }
  count_window_minima(seq, plen, K, prof_genome,
                      hvals, full_hvals, read_id, dq, nw
#if MINIMAX_TIE_AUDIT
                      , canonical_codes, profile_dq, tie_stats
#endif
                      );
#define NW_AT(pos) ((unsigned int)nw[(pos)])
#endif
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
      if (NW_AT(i) > 0) ncand++;

    if (ncand > 0)
      { Candidate *cands = (Candidate *) malloc((size_t)ncand * sizeof(Candidate));
        if (!cands) { perror("malloc cands"); exit(1); }
        { int ci = 0;
          for (int i = 0; i < plen; i++)
            if (NW_AT(i) > 0)
              { cands[ci].pos     = i;
                cands[ci].nw_val  = (int)NW_AT(i);
#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
                cands[ci].tuple =
                    position_tuple_key(i, prof_genome, hvals);
#else
                cands[ci].g_count = prof_genome[i];
                cands[ci].hash    = (uint64_t)hvals[i];
#endif
                ci++;
              }
        }
        if (ncand > 1)
          {
#if MINIMAX_HASH_REFINE
            candidate_full_hvals = full_hvals;
#endif
#if MINIMAX_PACK_POSITION_TUPLE && MINIMAX_CANONICAL_THINNING
            if (WSIZE <= UINT8_MAX)
              radix_sort_candidates(
                  cands, (size_t)ncand, candidate_scratch);
            else
#endif
              qsort(cands, (size_t)ncand, sizeof(Candidate), cand_cmp);
#if MINIMAX_HASH_REFINE
            candidate_full_hvals = NULL;
#endif
          }

#if MINIMAX_CANONICAL_THINNING
        {
          // Within each priority-equal group, order by position and visit the
          // two extremes as a pair.  Keeping the positional key out of the
          // primary qsort is intentional: reads contain many duplicate
          // priorities, and forcing a total order there costs more than these
          // small group-local sorts.  Accept both extremes when mutually
          // compatible, then move inward past candidates blocked by that
          // pair.  Left and right exchange under RC, so every decision is
          // strand-safe.
          int begin = 0;

          while (begin < ncand)
            {
              int end = begin + 1;
              int left;
              int right;

              while (end < ncand &&
                     candidate_priority_equal(cands + begin, cands + end,
                                              full_hvals))
                end++;
              if (end - begin > 1)
                qsort(cands + begin, (size_t)(end - begin),
                      sizeof(Candidate), cand_pos_cmp);
              left = begin;
              right = end - 1;
              while (left <= right)
                {
                  int left_pos;
                  int right_pos;

                  while (left <= right &&
                         seed_bit_is_set(
                             blocked, (size_t)cands[left].pos))
                    left++;
                  if (left > right)
                    break;
                  while (left < right &&
                         seed_bit_is_set(
                             blocked, (size_t)cands[right].pos))
                    right--;

                  left_pos = cands[left].pos;
                  right_pos = cands[right].pos;
                  if (left == right)
                    {
                      accept_sparse_seed(selected, blocked, (size_t)plen,
                                         (size_t)left_pos);
                      break;
                    }

                  if (right_pos - left_pos >= RSPACE)
                    {
                      accept_sparse_seed(selected, blocked, (size_t)plen,
                                         (size_t)left_pos);
                      accept_sparse_seed(selected, blocked, (size_t)plen,
                                         (size_t)right_pos);
                    }
                  left++;
                  right--;
                }
              begin = end;
            }
        }
#else
        for (int ci = 0; ci < ncand; ci++)
          { size_t p = (size_t)cands[ci].pos;
            if (seed_bit_is_set(blocked, p))
              continue;
            accept_sparse_seed(selected, blocked, (size_t)plen, p);
          }
#endif
        free(cands);
      }
  }

  if (EXCLUDE_PERIODIC_BP > 0)
    { exclude_long_periodic_seeds(
          seq-K+1, plen+K-1, K, plen, selected, periodic_stats);
    }

#ifdef DEBUG_NW
  // Print per-position window-count; mark final selected seeds with *****.
  fprintf(stderr, "NW read=%d plen=%d\n", read_id, plen);
  for (int i = 0; i < plen; i++)
    { const char *wkp = seq - K + 1 + i;
      char wkstr[65];
      for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
      wkstr[K] = '\0';
      fprintf(stderr, "  i=%-6d nw=%-5d g=%-5u h=%-4lu  %s  %s%s\n",
              i, (int)NW_AT(i),
              (unsigned)prof_genome[i], (unsigned long)hvals[i],
              wkstr, NW_AT(i) ? "*" : "",
              seed_bit_is_set(selected, (size_t)i) ? "****" : "");
    }
#endif

  // Step 4: append selected seeds to this block's KPOS run.
  int written = 0;
  { size_t word_index;

      // selected only contains valid positions.  Ascending words followed by
      // least-significant-set-bit removal therefore emits the exact same
      // position order as the historical i=0..plen-1 scan.
      for (word_index = 0; word_index < seed_nwords; word_index++)
        { uint64_t word = selected[word_index];

          while (word != 0)
            { unsigned int bit = (unsigned int)__builtin_ctzll(word);
              size_t pos = (word_index << 6) + bit;

              { const char *kp = seq - K + 1 + pos;
                uint64_t code = 0;
                int x;

                for (x = 0; x < K; x++)
                  code = (code << 2) | base2bit[(unsigned char)kp[x]];
                code &= KMASK;
                push_kpos(out, code, (int32_t)pos,
                          CN_PACK(read_id, cn_tier_of(prof_genome[pos])));
              }
              written++;
              word &= word - 1;
            }
        }
    }
  free(blocked);
  free(selected);
#if MINIMAX_CANONICAL_WINNOW
  free_window_counts(&nw);
#else
  free(nw);
  free(dq);
#endif
#if MINIMAX_TIE_AUDIT
  free(profile_dq);
  free(canonical_codes);
#endif
  free(full_hvals);
  free(hvals);
#undef NW_AT
  return written;
}

int find_seeds(const char    *seq,
               const int      plen,
               const int      K,
               const uint16  *prof_genome,
               const int      read_id,
               KposVector    *out
#if MINIMAX_TIE_AUDIT
               , TieStats    *tie_stats
#endif
               )
{
  CandidateRadixScratch candidate_scratch;
#if MINIMAX_CANONICAL_WINNOW
  CanonicalWindowScratch scratch;
  int written;

  init_candidate_radix_scratch(&candidate_scratch);
  init_canonical_window_scratch(&scratch);
  written = find_seeds_with_scratch(
      seq, plen, K, prof_genome, read_id, out, NULL,
      &candidate_scratch, &scratch);
  free_canonical_window_scratch(&scratch);
  free_candidate_radix_scratch(&candidate_scratch);
  return written;
#else
  int written;

  init_candidate_radix_scratch(&candidate_scratch);
  written = find_seeds_with_scratch(
      seq, plen, K, prof_genome, read_id, out, NULL,
      &candidate_scratch
#if MINIMAX_TIE_AUDIT
      , tie_stats
#endif
      );
  free_candidate_radix_scratch(&candidate_scratch);
  return written;
#endif
}

typedef struct
  { KposEntry *entry;
    size_t     len;
    int        ready;
    Periodic_Stats periodic_stats;
#if MINIMAX_TIE_AUDIT
    TieStats   tie_stats;
#endif
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

static void grow_profile(int need, int *capacity, uint16 **prof_genome)
{
  int next = (int)(1.2 * need) + 1000;
  uint16 *genome =
      (uint16 *)realloc(*prof_genome, (size_t)next * sizeof(uint16));

  if (genome == NULL)
    { perror("realloc profile");
      exit(1);
    }
  *prof_genome = genome;
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
  Profile_Index *Pg = NULL;
  char *seq;
  uint16 *prof_genome;
  int pmax = 20000;
  KposVector vector = {0};
  CandidateRadixScratch candidate_scratch;
#if MINIMAX_CANONICAL_WINNOW
  CanonicalWindowScratch window_scratch;
#endif

  db.bases = open_worker_bases(queue->db);
  Pg = Clone_Profiles(queue->genome);
  if (Pg == NULL)
    { fprintf(stderr,"copy_minimizer: cannot clone profile index\n");
      exit(1);
    }
  seq = New_Read_Buffer(&db);
  prof_genome = (uint16 *)malloc((size_t)pmax * sizeof(uint16));
  if (prof_genome == NULL)
    { perror("malloc profile");
      exit(1);
    }
  init_candidate_radix_scratch(&candidate_scratch);
#if MINIMAX_CANONICAL_WINNOW
  init_canonical_window_scratch(&window_scratch);
#endif

  while (1)
    { int block;
      int first, last, id;
      BlockResult *result;
      Periodic_Stats block_stats = {0};
#if MINIMAX_TIE_AUDIT
      TieStats block_tie_stats = {0};
#endif

      pthread_mutex_lock(&queue->lock);
      block = queue->next_block++;
      pthread_mutex_unlock(&queue->lock);
      if (block >= queue->end_block)
        break;

      first = queue->stub->ublocks[block];
      last  = queue->stub->ublocks[block + 1];
      vector.len = 0;

      for (id = first; id < last; id++)
        { int local_id = id-first;
          DAZZ_READ *r = db.reads+id;
          int plen = r->rlen - K + 1;
          int64 fetched;
          Load_Read(&db,id,seq,2);
          fetched = Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
          if (fetched > pmax)
            { grow_profile((int)fetched, &pmax, &prof_genome);
              Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
            }

          if (plen > 0)
            {
              find_seeds_with_scratch(
                  seq + K - 1, plen, K, prof_genome, local_id,
                  &vector, &block_stats,
                  &candidate_scratch
#if MINIMAX_CANONICAL_WINNOW
                  , &window_scratch
#endif
#if MINIMAX_TIE_AUDIT
                  , &block_tie_stats
#endif
                  );
            }
        }

      if (vector.len > 1)
        qsort(vector.entry, vector.len, sizeof(KposEntry), cmp_kpos);

      result = queue->result + block;
      pthread_mutex_lock(&queue->lock);
      result->entry = vector.entry;
      result->len   = vector.len;
      result->periodic_stats = block_stats;
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
  free_candidate_radix_scratch(&candidate_scratch);
#if MINIMAX_CANONICAL_WINNOW
  free_canonical_window_scratch(&window_scratch);
#endif
  free(prof_genome);
  free(seq - 1);
  fclose((FILE *)db.bases);
  Free_Profiles(Pg);
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

static void minimaxmer_print_periodic_stats_impl(
    const Periodic_Stats *stats, FILE *out, const char *prefix,
    int cutoff)
{
  const char *name = prefix == NULL ? "seed" : prefix;

  if (stats == NULL || out == NULL || cutoff <= 0)
    return;
  fprintf(out,
          "PERIODIC_MASK_V1\t%s\tcomplete_motif_bp_gt=%d\t"
          "overlap=any\tstage=post_thinning\t"
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
  // Optional flags precede the 3 positional args.
  { int i, j = 1;
    for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'k': K      = atoi(argv[i]+2); break;
          case 'w': WSIZE  = atoi(argv[i]+2); break;
          case 's': RSPACE = atoi(argv[i]+2); break;
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
            fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[i]);
            return 1;
        }
      } else
        argv[j++] = argv[i];
    }
    argc = j;
  }
  if (argc != 4 || K < 1 || K > 32 || WSIZE < 1 || RSPACE < 0 ||
      NTHREADS < 1) {
    fprintf(stderr, "Usage: %s [-T<int>] [-k<32>] [-w<100>] [-s<100>] [-l<int(-1)>] <db> <ref.prof> <out.kpos>\n", argv[0]);
    fprintf(stderr, "  -T<N>           : number of DB blocks processed in parallel (default 1)\n");
    fprintf(stderr, "  -k<N>           : k-mer length (<= 32)\n");
    fprintf(stderr, "  -w<N>           : minimizer window width  (larger => sparser seeds)\n");
    fprintf(stderr, "  -s<N>           : min spacing of uniform sparse seeds (larger => sparser)\n");
    fprintf(stderr, "  -l<N>           : exclude selected K-mers overlapping exact period-1/2/3 tracts whose complete-motif length is >N bp (-1 disables; default -1)\n");
    fprintf(stderr, "  <db>            : whole, partitioned DAZZLER dam/db\n");
    fprintf(stderr, "  <ref.prof>      : FastK genome k-mer count profile\n");
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
          "Seed selector: mode=copy_minimizer hash_mod=%llu hash_refine=%d robust_winnow=%d "
          "canonical_winnow=%d canonical_thinning=%d packed_tuple=%d "
          "hash_bytes=%zu candidate_bytes=%zu "
          "tie_audit=1\n",
          (unsigned long long)MOD, MINIMAX_HASH_REFINE,
          MINIMAX_ROBUST_WINNOW, MINIMAX_CANONICAL_WINNOW,
          MINIMAX_CANONICAL_THINNING, MINIMAX_PACK_POSITION_TUPLE,
          sizeof(PositionHash), sizeof(Candidate));
#else
  fprintf(stderr,
          "Seed selector: mode=copy_minimizer hash_mod=%llu hash_refine=%d robust_winnow=%d "
          "canonical_winnow=%d canonical_thinning=%d packed_tuple=%d "
          "hash_bytes=%zu candidate_bytes=%zu\n",
          (unsigned long long)MOD, MINIMAX_HASH_REFINE,
          MINIMAX_ROBUST_WINNOW, MINIMAX_CANONICAL_WINNOW,
          MINIMAX_CANONICAL_THINNING, MINIMAX_PACK_POSITION_TUPLE,
          sizeof(PositionHash), sizeof(Candidate));
#endif
  if (EXCLUDE_PERIODIC_BP > 0)
    fprintf(stderr,
            "Periodic seed exclusion: periods=1,2,3 "
            "complete_motif_bp>%d overlap=any stage=post_thinning\n",
            EXCLUDE_PERIODIC_BP);
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

  if (Pg->nreads != db->nreads)
    { fprintf(stderr,
              "%s: profile/DB read-count mismatch: DB=%d genome-profile=%d\n",
              argv[0], db->nreads, Pg->nreads);
      Free_Profiles(Pg);
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
        Free_Profiles(Pg);
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
        Free_Profiles(Pg);
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
        Free_Profiles(Pg);
        Close_DB(db);
        return 1;
      }
    if (stub->ublocks[0] != 0 ||
        stub->ublocks[stub->nblocks] != db->nreads)
      { fprintf(stderr,
                "%s: DB block boundaries do not span all %d input reads\n",
                argv[0], db->nreads);
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
      if ((uint64_t)(stub->ublocks[b + 1] - stub->ublocks[b]) > CN_MAXVAL)
        { fprintf(stderr,
                  "%s: block %d exceeds the copy-tier rid limit (2^29 reads)\n",
                  argv[0], b + 1);
          Free_DB_Stub(stub);
          Free_Profiles(Pg);
          Close_DB(db);
          return 1;
        }
  }

  {
    const char *output_path = argv[3];
    char *index_path = add_suffix(output_path, ".idx");
    char *output_tmp = add_suffix(output_path, ".tmp");
    char *index_tmp = add_suffix(index_path, ".tmp");
    FILE *out = fopen(output_tmp, "wb");
    FILE *idx;
    WorkQueue queue;
    BlockResult *result;
    pthread_t *thread;
    uint64_t total_written = 0;
    Periodic_Stats total_periodic_stats = {0};
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
    if (write_kpos_header_for_k(out,K) != 0 ||
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

    memset(&queue,0,sizeof(queue));
    queue.db         = db;
    queue.stub       = stub;
    queue.genome     = Pg;
    queue.result     = result;
    queue.next_block = 0;
    queue.end_block  = stub->nblocks;
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

        periodic_stats_add(&total_periodic_stats,&block->periodic_stats);
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
    minimaxmer_print_periodic_stats_impl(
        &total_periodic_stats,stderr,"standalone_seed",EXCLUDE_PERIODIC_BP);
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
  Free_DB_Stub(stub);
  Close_DB(db);

  return 0;
}
