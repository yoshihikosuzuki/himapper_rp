#ifndef HIMAPPER_MAPQ_H
#define HIMAPPER_MAPQ_H

#include <stdint.h>
#include <stddef.h>

/* Disabled by default.  The caller emits 255 when MAPQ_MODE is zero. */

extern int MAPQ_MODE;

#define MAPQ_BIN 64

/* Detailed discarded-chain storage is bounded.  Preserve the largest score
   outside that storage, and mark every query bin touched by such a chain.
   bits contains ceil(length/(64*MAPQ_BIN)) words, private to this read. */

typedef struct
  { const uint64_t *bits;
    int             length, max_score;
  } Mapq_Overflow;

typedef struct
  { int      chain_score;
    int      seed_count;
    uint64_t chain_id;
  } Mapq_Evidence;

typedef struct
  { int      qb, qe;
    int      rb, re;
    int      reference, reverse;
    int      chain_score, seed_count;
    uint64_t chain_id;
  } Mapq_Shadow;

#define MAPQ_FRAGMENT_MAX 65536
#define MAPQ_JOIN_MAX 32
#define MAPQ_JOIN_POINTS 8
#define MAPQ_TRACE_MAX_WIDTH 32

typedef struct
  { int query, reference;
  } Mapq_Point;

/* Actual fragment endpoints, sampled across query space; never interpolated. */
typedef struct
  { int        n;
    Mapq_Point points[MAPQ_JOIN_POINTS];
  } Mapq_Join_Geometry;

typedef struct
  { int          n, cap;
    Mapq_Shadow *items;
    Mapq_Join_Geometry *geometry;
  } Mapq_Join_Catalog;

typedef struct
  { uint64_t input, skipped, joins;
    double   cpu_seconds;
    size_t   allocated_bytes;
  } Mapq_Fragment_Stats;

typedef struct Mapq_Fragment_Work Mapq_Fragment_Work;

typedef struct Mapq_Fragment_Pass Mapq_Fragment_Pass;

/* Retain bounded per-pass inputs until alignment query intervals are known.
   Each pass keeps its original DP boundary and its own 65536-input limit. */
typedef struct
  { Mapq_Shadow        *items;
    Mapq_Fragment_Pass *passes;
    size_t              n, cap;
    int                 npasses, pass_cap;
  } Mapq_Fragment_Store;

int Mapq_Fragment_Save(Mapq_Fragment_Work *work, int kmer, int max_gap,
                       int drift_penalty_x100, Mapq_Fragment_Store *store);
void Mapq_Free_Fragment_Store(Mapq_Fragment_Store *store);
size_t Mapq_Fragment_Store_Bytes(const Mapq_Fragment_Store *store);


/* Reuse one scratch per worker, resetting inputs between reads/passes.  Reset
   preserves allocations and cumulative statistics.  Input/skipped count Add
   calls and rejected inputs.  In Assign_Retained, joins counts all evaluated
   multi-fragment DP prefixes; the compatibility Build counts only terminal
   aggregates.  Both include duplicates from the two shifted diagonal bands.
   Evaluation CPU time is measured only when timing is enabled.  Add, Save
   and Build return zero on allocation failure; invalid/over-limit Add inputs
   are skipped successfully.  Save omits passes with fewer than two unique
   fragments, which cannot supply joined evidence.
   The caller must retain its existing individual-shadow/overflow handling. */

Mapq_Fragment_Work *Mapq_New_Fragment_Work(int timing);
void Mapq_Fragment_Reset(Mapq_Fragment_Work *work);
int  Mapq_Fragment_Add(Mapq_Fragment_Work *work, const Mapq_Shadow *fragment);
/* Compatibility catalog for focused diagnostics/tests.  Production uses
   Save/Assign_Retained and does not select a global top-32 joined catalog. */
int  Mapq_Fragment_Build(Mapq_Fragment_Work *work, int kmer, int max_gap,
                         int drift_penalty_x100, Mapq_Join_Catalog *out);
Mapq_Fragment_Stats Mapq_Fragment_Statistics(const Mapq_Fragment_Work *work);
void Mapq_Free_Fragment_Work(Mapq_Fragment_Work *work);
void Mapq_Free_Join_Catalog(Mapq_Join_Catalog *catalog);
size_t Mapq_Join_Bytes(const Mapq_Join_Catalog *catalog);

typedef struct
  { double score, identity;
    int    chain_score, seeds;
    int    secondary, incomplete;
  } Mapq_Input;

typedef struct
  { double align_score;
    int    chain_score;
    int    aligned, n_near;
    int    drop_score;
  } Mapq_Competition;

/* A bounded, uncalibrated confidence estimate.  Scores must exclude the unique
   seed bonus.  Missing main evidence, a tied/better aligned competitor,
   secondary mappings, and invalid/incomplete evidence receive zero.
   identity and unique_ratio must be finite and in [0,1].  An absent aligned
   competitor has aligned=0; an absent unaligned competitor has drop_score=0.
   Valid positive main scores strictly exceeding every aligned competitor
   (or with no aligned competitor) have a floor of 1, including after the
   discarded-chain and joined-fragment caps.  Other results remain in [0,60]. */

int Mapq_Calculate(const Mapq_Input *input, const Mapq_Competition *competition,
                   int min_chain_score, double unique_ratio);

/* Pairwise competition on nonnegative, half-open query intervals.  With
   overlap O and lengths L <= H, require O/L - (L-O)/H > 0.5.  This symmetric
   approximation retains containment, but discounts partly overlapping split
   segments.  It does not use the union of all primary query intervals. */

int Mapq_Query_Overlap(int ab, int ae, int bb, int be);

/* Lower-bound identity estimate from the existing edit count, without CIGAR.
   Invalid/empty spans or negative edit counts return zero. */

double Mapq_Identity(int query_span, int ref_span, int diffs);

/* Keep the strongest source chain; never add duplicated chain support.
   Ties use the larger seed count, then the smaller stable chain ID. */

void Mapq_Combine_Evidence(Mapq_Evidence *dst, const Mapq_Evidence *src);

typedef struct
  { int           qb, qe, rb, re;
    int           reference, reverse;
    int           next, root, incomplete;
    double        score, local_score, identity;
    Mapq_Evidence evidence;
    const uint16_t *trace;
    int           trace_length, trace_spacing;
  } Mapq_Path;

typedef struct Mapq_Work Mapq_Work;

/* Reusable per-worker catalog.  Paths() resets the catalog and returns n
   zeroed entries with next=-1.  Ranks() supplies an index permutation, initially
   0..n-1, which the caller fills in its normal selection order before Assign().
   All root paths, including alternatives omitted from output, must be present.
   next links must be acyclic.  Invalid links/rank permutations fail closed.
   At most 32 discarded-chain summaries are accepted; larger inventories fail
   closed, matching the caller's explicit per-read information limit.
   Additional discarded chains use overflow, or NULL when none were lost.
   Their maximum score caps confidence on any intersecting 64-base bin; it
   may conservatively include unrelated or own-source chains.  A complete
   chain's interval includes gaps between its linked segments.  Invalid
   overflow descriptors fail closed.  Joined fragments impose a separate cap
   after this original calculation.  That cap compares competitor support with
   max(original support, matching own joined support); matching own support
   requires reference/strand agreement and reciprocal 80% query/reference
   coverage plus consistent sampled query/reference correspondence.  Borrowed
   trace pairs hold edit count and reference advance on the query spacing grid.
   Missing/ambiguous correspondence supplies neither own nor rival evidence.
   It never changes the original support penalties or increases Q. */

Mapq_Work *Mapq_New_Work(void);
Mapq_Path *Mapq_Paths(Mapq_Work *work, int n);
int       *Mapq_Ranks(Mapq_Work *work);
void       Mapq_Assign(Mapq_Work *work, int n, const Mapq_Shadow *shadows,
                        int nshadows, int min_chain_score, int repeat_bases,
                        const Mapq_Overflow *overflow, const Mapq_Join_Catalog *joined);

/* Production path: compute the unchanged raw MAPQ, then stream every retained
   multi-fragment DP prefix.  Query applicability precedes score selection.
   Only maximum eligible own/other support is kept per local/complete chain;
   there is no global joined-candidate limit.  Reuse one scratch per reporting
   worker.  Its statistics record streamed prefixes and report-stage CPU time;
   malformed retained inputs fail closed without altering mapping output. */
void Mapq_Assign_Retained(Mapq_Work *work, int n, const Mapq_Shadow *shadows,
                           int nshadows, int min_chain_score, int repeat_bases,
                           const Mapq_Overflow *overflow, const Mapq_Fragment_Store *store,
                           Mapq_Fragment_Work *scratch);


/* Secondary roots receive zero.  For a primary root, use the lower of its
   complete-chain confidence and the emitted segment's local confidence.
   Local duplicates require identical path geometry; complete-chain duplicates
   require the same entire sequence of path geometries.  A shared tail does
   not erase competing heads.  Primary query intervals are the prefixes left
   after earlier primary chains consume shared paths.  Secondary consumption
   and the caller's output filters do not affect confidence assignment. */

int    Mapq_Quality(const Mapq_Work *work, int root, int segment);
size_t Mapq_Work_Bytes(const Mapq_Work *work);
void   Mapq_Free_Work(Mapq_Work *work);

#endif
