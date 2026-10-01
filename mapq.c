#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "mapq.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

int MAPQ_MODE = 0;

static double clip_ratio(double value)
{ if (value < 0.)
    return 0.;
  if (value > 1.)
    return 1.;
  return value;
}

static double mapq_scale(const Mapq_Input *input, double unique_ratio)
{ double support = clip_ratio(input->chain_score / 100.) * unique_ratio;
  double seed_support = clip_ratio(input->seeds / 10.);

  if (support > seed_support)
    support = seed_support;
  return 40. * input->identity * support * log(fmax(input->score,2.));
}

int Mapq_Calculate(const Mapq_Input *input, const Mapq_Competition *competition,
                   int min_chain_score, double unique_ratio)
{ double identity, scale, q;
  double chain_ratio;
  int    second;

  if (!input || !competition)
    return 0;
  if (input->secondary || input->incomplete || input->chain_score <= 0 ||
      input->seeds <= 0 || !isfinite(input->score) || input->score <= 0. ||
      !isfinite(input->identity) || input->identity < 0. || input->identity > 1. ||
      !isfinite(unique_ratio) || unique_ratio < 0. || unique_ratio > 1. ||
      min_chain_score < 0 || competition->chain_score < 0 ||
      competition->drop_score < 0 || competition->n_near < 0)
    return 0;
  if (competition->aligned && (!isfinite(competition->align_score) ||
                               input->score <= competition->align_score))
    return 0;

  identity = input->identity;
  scale = mapq_scale(input,unique_ratio);

  second = competition->chain_score;
  if (second < min_chain_score)
    second = min_chain_score;
  chain_ratio = clip_ratio((double) second / input->chain_score);

  if (competition->aligned)
    { double x   = chain_ratio * clip_ratio(competition->align_score / input->score);
      double cap = 6.02 * identity * identity * (input->score - competition->align_score);

      q = scale * (1. - x*x);
      if (q > cap)
        q = cap;
    }
  else
    q = scale * (1. - chain_ratio);
  q -= 4.343 * log1p((double) competition->n_near);

  if (competition->drop_score > 0)
    { double ratio = clip_ratio((double) competition->drop_score / input->chain_score);
      double cap   = scale * (1. - ratio);

      if (q > cap)
        q = cap;
    }

  if (q <= 0.)
    return 0;
  if (q >= 60.)
    return 60;
  return (int) (q + .499);
}

int Mapq_Query_Overlap(int ab, int ae, int bb, int be)
{ int64_t low, high, overlap;

  if (ab < 0 || bb < 0 || ae <= ab || be <= bb || ae <= bb || be <= ab)
    return 0;
  low     = (int64_t) ae - ab;
  high    = (int64_t) be - bb;
  overlap = (int64_t) (ae < be ? ae : be) - (ab > bb ? ab : bb);
  if (low > high)
    { int64_t swap = low;

      low  = high;
      high = swap;
    }
  return (double) overlap / low - (double) (low-overlap) / high > .5;
}

double Mapq_Identity(int query_span, int ref_span, int diffs)
{ int span;

  if (query_span <= 0 || ref_span <= 0 || diffs < 0)
    return 0.;
  span = query_span > ref_span ? query_span : ref_span;
  if (diffs >= span)
    return 0.;
  return 1. - (double) diffs / span;
}

void Mapq_Combine_Evidence(Mapq_Evidence *dst, const Mapq_Evidence *src)
{ if (src->chain_score > dst->chain_score ||
      (src->chain_score == dst->chain_score &&
       (src->seed_count > dst->seed_count ||
        (src->seed_count == dst->seed_count && src->chain_id < dst->chain_id))))
    *dst = *src;
}

typedef struct
  { int           qb, qe, incomplete;
    int           rb, re, reference, reverse;
    int           primary_qb, primary_qe;
    uint32_t      own_shadows, own_joined, other_joined, near_joined;
    uint32_t      local_own_joined, local_other_joined;
    uint8_t       covered[MAPQ_JOIN_MAX], possible[MAPQ_JOIN_MAX], narrow[MAPQ_JOIN_MAX];
    double        identity;
    Mapq_Evidence evidence;
  } Mapq_Summary;

typedef struct
  { uint8_t covered, possible, narrow, nearby;
  } Mapq_Relation;

typedef struct
  { int own, other;
  } Mapq_Support;

struct Mapq_Work
  { Mapq_Path    *paths;
    Mapq_Summary *summary;
    int          *ranks, *geometry, *suffix, *state, *stack;
    int          *primary, *component_primary, *seen;
    int          *group_quality, *local_quality;
    int          *hash;
    int           n, cap;
    size_t        hash_cap;
    int          *trace_ref;
    size_t       *trace_offset;
    size_t        trace_n, trace_cap;
    Mapq_Relation *joined_local, *joined_group;
    Mapq_Support  *local_support, *group_support;
    int           *local_stamp, *group_stamp, *support_targets;
    size_t         ntargets;
    int            joined_stamp, assigned;
    double         unique_ratio;
  };

static void *mapq_resize(void *ptr, size_t n, size_t size)
{ void *result;

  if (n > SIZE_MAX/size)
    { fputs("MAPQ allocation size overflow\n",stderr);
      exit(1);
    }
  result = realloc(ptr,n*size);
  if (!result)
    { perror("MAPQ allocation");
      exit(1);
    }
  return result;
}

Mapq_Work *Mapq_New_Work(void)
{ Mapq_Work *work = mapq_resize(NULL,1,sizeof(*work));

  memset(work,0,sizeof(*work));
  return work;
}

Mapq_Path *Mapq_Paths(Mapq_Work *work, int n)
{ if (n < 0)
    { work->n = 0;
      return NULL;
    }
  if (n > work->cap)
    { int    cap = n > INT_MAX-64 ? n : n+64;
      size_t hash_cap = 128;

      while (hash_cap < (size_t) cap*2)
        hash_cap *= 2;
      work->paths             = mapq_resize(work->paths,cap,sizeof(*work->paths));
      work->summary           = mapq_resize(work->summary,cap,sizeof(*work->summary));
      work->trace_offset      = mapq_resize(work->trace_offset,cap,sizeof(*work->trace_offset));
      work->joined_local      = mapq_resize(work->joined_local,cap,sizeof(*work->joined_local));
      work->joined_group      = mapq_resize(work->joined_group,cap,sizeof(*work->joined_group));
      work->local_support     = mapq_resize(work->local_support,cap,sizeof(*work->local_support));
      work->group_support     = mapq_resize(work->group_support,cap,sizeof(*work->group_support));
      work->local_stamp       = mapq_resize(work->local_stamp,cap,sizeof(int));
      work->group_stamp       = mapq_resize(work->group_stamp,cap,sizeof(int));
      work->support_targets   = mapq_resize(work->support_targets,(size_t) cap*2,sizeof(int));
      work->ranks             = mapq_resize(work->ranks,cap,sizeof(int));
      work->geometry          = mapq_resize(work->geometry,cap,sizeof(int));
      work->suffix            = mapq_resize(work->suffix,cap,sizeof(int));
      work->state             = mapq_resize(work->state,cap,sizeof(int));
      work->stack             = mapq_resize(work->stack,cap,sizeof(int));
      work->primary           = mapq_resize(work->primary,cap,sizeof(int));
      work->component_primary = mapq_resize(work->component_primary,cap,sizeof(int));
      work->seen              = mapq_resize(work->seen,cap,sizeof(int));
      work->group_quality     = mapq_resize(work->group_quality,cap,sizeof(int));
      work->local_quality     = mapq_resize(work->local_quality,cap,sizeof(int));
      work->hash              = mapq_resize(work->hash,hash_cap,sizeof(int));
      work->cap = cap;
      work->hash_cap = hash_cap;
    }
  work->n = n;
  if (n)
    { memset(work->paths,0,(size_t) n*sizeof(*work->paths));
      memset(work->group_quality,0,(size_t) n*sizeof(int));
      memset(work->local_quality,0,(size_t) n*sizeof(int));
    }
  for (int i = 0; i < n; ++i)
    { work->paths[i].next = -1;
      work->ranks[i] = i;
    }
  return work->paths;
}

int *Mapq_Ranks(Mapq_Work *work)
{ return work->ranks;
}

static uint64_t mapq_hash(uint64_t x)
{ x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

static size_t mapq_geometry_hash(const Mapq_Path *path)
{ uint64_t hash = mapq_hash((uint32_t) path->qb | (uint64_t)(uint32_t) path->qe << 32);

  hash ^= mapq_hash((uint32_t) path->rb | (uint64_t)(uint32_t) path->re << 32);
  hash ^= mapq_hash((uint32_t) path->reference | (uint64_t)(uint32_t) path->reverse << 32);
  return (size_t) hash;
}

static int mapq_same_geometry(const Mapq_Path *a, const Mapq_Path *b)
{ return a->qb == b->qb && a->qe == b->qe && a->rb == b->rb && a->re == b->re &&
         a->reference == b->reference && a->reverse == b->reverse;
}

static int mapq_reciprocal_overlap(int ab, int ae, int bb, int be)
{ int64_t overlap = (int64_t) (ae < be ? ae : be) - (ab > bb ? ab : bb);

  return ae > ab && be > bb && overlap > 0 &&
         overlap*5 >= ((int64_t) ae-ab)*4 && overlap*5 >= ((int64_t) be-bb)*4;
}

/* Sparse traces bound the reference position at a query coordinate.  Build a
   prefix once per inspected path; no trace storage is touched without joined
   candidates on the same reference/strand.  Invalid traces remain unknown. */

static int mapq_trace_prefix(Mapq_Work *work, int at)
{ const Mapq_Path *path = work->paths+at;
  int             count, q, r;
  size_t          need;

  if (work->trace_offset[at] != SIZE_MAX)
    return work->trace_offset[at] != SIZE_MAX-1;
  work->trace_offset[at] = SIZE_MAX-1;
  if (!path->trace || path->trace_spacing <= 0 || path->qb < 0 || path->qe <= path->qb ||
      path->rb < 0 || path->re <= path->rb || path->trace_length <= 0 || path->trace_length%2)
    return 0;
  count = (path->qe-1)/path->trace_spacing-path->qb/path->trace_spacing+1;
  if (count != path->trace_length/2)
    return 0;
  need = work->trace_n+(size_t) count+1;
  if (need < work->trace_n)
    return 0;
  if (need > work->trace_cap)
    { size_t cap = need+need/2+64;

      if (cap < need)
        cap = need;
      work->trace_ref = mapq_resize(work->trace_ref,cap,sizeof(int));
      work->trace_cap = cap;
    }
  q = path->qb;
  r = path->rb;
  work->trace_ref[work->trace_n] = r;
  for (int i = 0; i < count; ++i)
    { int64_t end = ((int64_t) q/path->trace_spacing+1)*path->trace_spacing;
      int     advance = path->trace[2*i+1], edits = path->trace[2*i];
      int64_t delta;

      if (end > path->qe)
        end = path->qe;
      delta = advance-(end-q);
      if (llabs(delta) > edits || (int64_t) r+advance > path->re)
        return 0;
      q = (int) end;
      r += advance;
      work->trace_ref[work->trace_n+i+1] = r;
    }
  if (r != path->re)
    return 0;
  work->trace_offset[at] = work->trace_n;
  work->trace_n = need;
  return 1;
}

static void mapq_trace_interval(const Mapq_Path *path, const int *prefix, int tile,
                                 int query, int64_t *lo, int64_t *hi)
{ int64_t q0 = ((int64_t) path->qb/path->trace_spacing+tile)*path->trace_spacing;
  int64_t q1 = q0+path->trace_spacing;
  int64_t x, length, delta, low, high, edits = path->trace[2*tile];

  if (!tile)
    q0 = path->qb;
  if (q1 > path->qe)
    q1 = path->qe;
  x = query-q0;
  length = (int64_t) prefix[tile+1]-prefix[tile];
  delta = length-(q1-q0);
  /* delta-edits <= 0, so C truncation gives ceil for the lower bound. */
  low  = x+(delta-edits)/2;
  high = x+(delta+edits)/2;
  if (low < 0) low = 0;
  if (high > length) high = length;
  *lo = prefix[tile]+low;
  *hi = prefix[tile]+high;
}

static Mapq_Relation mapq_joined_relation(Mapq_Work *work, int at,
                                           const Mapq_Shadow *item,
                                           const Mapq_Join_Geometry *geometry)
{ const Mapq_Path *path = work->paths+at;
  Mapq_Relation   relation = {0};
  int same_frame = path->reference == item->reference && path->reverse == item->reverse;
  int nearby = same_frame && path->rb < item->re && item->rb < path->re;
  int traced = nearby && mapq_trace_prefix(work,at);

  relation.nearby = (uint8_t) nearby;
  if (!geometry || geometry->n < 2 || geometry->n > MAPQ_JOIN_POINTS)
    return relation;
  for (int j = 0; j < geometry->n; ++j)
    { const Mapq_Point *point = geometry->points+j;
      uint8_t bit = (uint8_t) (1U << j);
      int     tile, count;
      int64_t lo, hi;
      const int *prefix;

      if (point->query < path->qb || point->query > path->qe)
        continue;
      relation.covered |= bit;
      if (!same_frame)
        continue;
      if (!traced)
        { if (nearby)
            relation.possible |= bit;
          continue;
        }
      count = path->trace_length/2;
      tile = point->query/path->trace_spacing-path->qb/path->trace_spacing;
      if (tile == count)
        --tile;
      prefix = work->trace_ref+work->trace_offset[at];
      mapq_trace_interval(path,prefix,tile,point->query,&lo,&hi);
      /* At a grid boundary an indel may belong to either adjacent tile. */
      if (tile > 0 && (int64_t) point->query ==
          ((int64_t) path->qb/path->trace_spacing+tile)*path->trace_spacing)
        { int64_t before_lo, before_hi;

          mapq_trace_interval(path,prefix,tile-1,point->query,&before_lo,&before_hi);
          if (before_lo < lo) lo = before_lo;
          if (before_hi > hi) hi = before_hi;
        }
      if (point->reference >= lo && point->reference <= hi)
        { relation.possible |= bit;
          if (hi-lo <= MAPQ_TRACE_MAX_WIDTH)
            relation.narrow |= bit;
        }
    }
  return relation;
}

static void mapq_joined_path(Mapq_Work *work, int at, const Mapq_Join_Catalog *joined)
{ Mapq_Summary *summary = work->summary+at;

  summary->near_joined = 0;
  memset(summary->covered,0,sizeof(summary->covered));
  memset(summary->possible,0,sizeof(summary->possible));
  memset(summary->narrow,0,sizeof(summary->narrow));
  if (!joined)
    return;
  for (int i = 0; i < joined->n; ++i)
    { Mapq_Relation relation = mapq_joined_relation(work,at,joined->items+i,
                                                   joined->geometry ? joined->geometry+i : NULL);

      if (relation.nearby)
        summary->near_joined |= UINT32_C(1) << i;
      summary->covered[i] = relation.covered;
      summary->possible[i] = relation.possible;
      summary->narrow[i] = relation.narrow;
    }
}

static void mapq_joined_classify(Mapq_Summary *summary, const Mapq_Join_Catalog *joined)
{ summary->own_joined = summary->other_joined = 0;
  if (!joined)
    return;
  for (int i = 0; i < joined->n; ++i)
    { uint32_t bit = UINT32_C(1) << i;
      const Mapq_Join_Geometry *geometry = joined->geometry ? joined->geometry+i : NULL;

      if (!(summary->near_joined & bit) || (summary->covered[i] & ~summary->possible[i]))
        summary->other_joined |= bit;
      else if (geometry && geometry->n >= 2 && geometry->n <= MAPQ_JOIN_POINTS &&
               summary->narrow[i] == (1U << geometry->n)-1)
        summary->own_joined |= bit;
    }
}

/* Evaluate each linked suffix once, including suffixes shared by roots. */

static int mapq_summaries(Mapq_Work *work, const Mapq_Shadow *shadows, int nshadows,
                           const Mapq_Join_Catalog *joined)
{ int n = work->n;

  memset(work->state,0,(size_t) n*sizeof(int));
  memset(work->hash,0,work->hash_cap*sizeof(int));
  for (int i = 0; i < n; ++i)
    if (!work->state[i])
      { int at = i, depth = 0;

        while (at >= 0)
          { if (at >= n || work->state[at] == 1)
              return 0;
            if (work->state[at] == 2)
              break;
            work->state[at] = 1;
            work->stack[depth++] = at;
            at = work->paths[at].next;
          }
        while (depth)
          { Mapq_Path    *path;
            Mapq_Summary *summary;
            size_t        slot;
            int           tail;

            at = work->stack[--depth];
            path = work->paths+at;
            summary = work->summary+at;
            summary->qb         = path->qb;
            summary->qe         = path->qe;
            summary->rb         = path->rb;
            summary->re         = path->re;
            summary->reference  = path->reference;
            summary->reverse    = path->reverse;
            summary->identity   = path->identity;
            summary->incomplete = path->incomplete;
            summary->evidence   = path->evidence;
            summary->own_shadows = 0;
            mapq_joined_path(work,at,joined);
            mapq_joined_classify(summary,joined);
            summary->local_own_joined = summary->own_joined;
            summary->local_other_joined = summary->other_joined;
            for (int j = 0; j < nshadows; ++j)
              if (shadows[j].chain_id == path->evidence.chain_id)
                summary->own_shadows |= UINT32_C(1) << j;
            if (path->qb < 0 || path->qe <= path->qb || path->rb < 0 ||
                path->re <= path->rb || !isfinite(path->identity) ||
                path->identity < 0. || path->identity > 1. || path->next < -1)
              summary->incomplete = 1;
            if (path->next >= 0)
              { Mapq_Summary *tail = work->summary+path->next;

                if (summary->qb > tail->qb)
                  summary->qb = tail->qb;
                if (summary->qe < tail->qe)
                  summary->qe = tail->qe;
                if (summary->rb > tail->rb)
                  summary->rb = tail->rb;
                if (summary->re < tail->re)
                  summary->re = tail->re;
                if (summary->reference != tail->reference || summary->reverse != tail->reverse)
                  summary->reference = -1;
                if (summary->identity > tail->identity)
                  summary->identity = tail->identity;
                summary->incomplete |= tail->incomplete;
                summary->own_shadows |= tail->own_shadows;
                summary->near_joined |= tail->near_joined;
                if (joined)
                  for (int j = 0; j < joined->n; ++j)
                    { summary->covered[j] |= tail->covered[j];
                      summary->possible[j] |= tail->possible[j];
                      summary->narrow[j] |= tail->narrow[j];
                    }
                mapq_joined_classify(summary,joined);
                Mapq_Combine_Evidence(&summary->evidence,&tail->evidence);
              }
            /* Two roots represent the same chain only if every path geometry
               in their complete suffix agrees.  Sharing a tail is insufficient. */
            tail = path->next >= 0 ? work->suffix[path->next] : -1;
            slot = (size_t) mapq_hash((uint32_t) work->geometry[at] |
                                     (uint64_t)(uint32_t) tail << 32) & (work->hash_cap-1);
            while (work->hash[slot])
              { int other = work->hash[slot]-1;
                int next = work->paths[other].next;
                int other_tail = next >= 0 ? work->suffix[next] : -1;

                if (work->geometry[other] == work->geometry[at] && other_tail == tail)
                  break;
                slot = (slot+1) & (work->hash_cap-1);
              }
            if (!work->hash[slot])
              work->hash[slot] = at+1;
            work->suffix[at] = work->hash[slot]-1;
            work->state[at] = 2;
          }
      }
  return 1;
}

/* Primary support counts a source chain only once, using its strongest saved
   score.  Shared tails and duplicate paths must not inflate unique_ratio. */

static size_t mapq_source_slot(const Mapq_Work *work, uint64_t source)
{ size_t slot = (size_t) mapq_hash(source);

  slot &= work->hash_cap-1;
  while (work->hash[slot])
    { int path = work->hash[slot]-1;

      if (work->paths[path].evidence.chain_id == source)
        break;
      slot = (slot+1) & (work->hash_cap-1);
    }
  return slot;
}

static void mapq_competitor(Mapq_Work *work, int stamp, int component,
                             double score, const Mapq_Evidence *evidence,
                             const Mapq_Input *input, Mapq_Competition *competition)
{ if (!competition->aligned || score > competition->align_score)
    competition->align_score = score;
  competition->aligned = 1;
  if (competition->chain_score < evidence->chain_score)
    competition->chain_score = evidence->chain_score;
  if (work->seen[component] != stamp &&
      (evidence->seed_count >= input->seeds || input->score-score <= 4.))
    { work->seen[component] = stamp;
      ++competition->n_near;
    }
}

static int mapq_dropped(int qb, int qe, uint32_t own_shadows,
                         const Mapq_Shadow *shadows, int nshadows)
{ int score = 0;

  for (int i = 0; i < nshadows; ++i)
    { const Mapq_Shadow *shadow = shadows+i;

      if (shadow->chain_score > score && Mapq_Query_Overlap(qb,qe,shadow->qb,shadow->qe) &&
          !(own_shadows & (UINT32_C(1) << i)))
        score = shadow->chain_score;
    }
  return score;
}

static int mapq_joined_score(int qb, int qe, uint32_t other_joined,
                              const Mapq_Join_Catalog *joined)
{ int score = 0;

  if (joined)
    for (int i = 0; i < joined->n; ++i)
      { const Mapq_Shadow *item = joined->items+i;

        if ((other_joined & (UINT32_C(1) << i)) && item->chain_score > score &&
            mapq_reciprocal_overlap(qb,qe,item->qb,item->qe) &&
            Mapq_Query_Overlap(qb,qe,item->qb,item->qe))
          score = item->chain_score;
      }
  return score;
}

/* The new comparison must use like-for-like joined support when it exists on
   the selected locus.  Keep all original alignment scores, support penalties,
   and raw-shadow caps unchanged; only this new cap's denominator can increase. */

static int mapq_own_joined_score(int qb, int qe, int rb, int re, int reference, int reverse,
                                 uint32_t own_joined, const Mapq_Join_Catalog *joined)
{ int score = 0;

  if (joined && reference >= 0)
    for (int i = 0; i < joined->n; ++i)
      { const Mapq_Shadow *item = joined->items+i;

        if ((own_joined & (UINT32_C(1) << i)) && item->chain_score > score &&
            item->reference == reference && item->reverse == reverse &&
            mapq_reciprocal_overlap(qb,qe,item->qb,item->qe) &&
            mapq_reciprocal_overlap(rb,re,item->rb,item->re))
          score = item->chain_score;
      }
  return score;
}

static int mapq_joined_cap(int original, const Mapq_Input *input, double unique_ratio,
                            int own_score, int competitor_score)
{ double q;
  int    denominator = input->chain_score;

  if (!original || competitor_score <= 0)
    return original;
  if (denominator < own_score)
    denominator = own_score;
  q = mapq_scale(input,unique_ratio) *
      (1. - clip_ratio((double) competitor_score / denominator));
  if (q <= 0.)
    return 0;
  if (q >= original)
    return original;
  return (int) (q + .499);
}

/* The maximum score is an upper bound for every discarded chain whose
   detailed record was lost.  Query-bin intersection deliberately errs toward
   including candidates: exact query overlap and source IDs are unavailable. */

static int mapq_overflow_score(int qb, int qe, const Mapq_Overflow *overflow)
{ int first, last;

  if (!overflow || !overflow->max_score)
    return 0;
  if (qb < 0)
    qb = 0;
  if (qe > overflow->length)
    qe = overflow->length;
  if (qb >= qe)
    return 0;
  first = qb/MAPQ_BIN;
  last  = (qe-1)/MAPQ_BIN;
  for (int word = first/64; word <= last/64; ++word)
    { uint64_t mask = UINT64_MAX;

      if (word == first/64)
        mask &= UINT64_MAX << (first%64);
      if (word == last/64)
        mask &= UINT64_MAX >> (63-last%64);
      if (overflow->bits[word] & mask)
        return overflow->max_score;
    }
  return 0;
}

void Mapq_Assign(Mapq_Work *work, int n, const Mapq_Shadow *shadows,
                  int nshadows, int min_chain_score, int repeat_bases,
                  const Mapq_Overflow *overflow, const Mapq_Join_Catalog *joined)
{ double sum = 0., unique_ratio;
  int    nparts = 0;

  if (work) work->assigned = 0;
  if (!work || !work->n)
    return;
  memset(work->group_quality,0,(size_t) work->n*sizeof(int));
  memset(work->local_quality,0,(size_t) work->n*sizeof(int));
  if (n != work->n || nshadows < 0 || nshadows > 32 ||
      (nshadows && !shadows) || repeat_bases < 0 ||
      (joined && (joined->n < 0 || joined->n > MAPQ_JOIN_MAX ||
                  joined->cap < joined->n || (joined->n && !joined->items))) ||
      (overflow && (overflow->max_score < 0 || overflow->length < 0 ||
                    (overflow->max_score && (!overflow->bits || !overflow->length)))))
    return;
  memset(work->seen,0,(size_t) n*sizeof(int));
  for (int i = 0; i < n; ++i)
    { int rank = work->ranks[i];

      if (rank < 0 || rank >= n || work->seen[rank])
        return;
      work->seen[rank] = 1;
    }
  memset(work->hash,0,work->hash_cap*sizeof(int));
  for (int i = 0; i < n; ++i)
    { size_t slot = mapq_geometry_hash(work->paths+i) & (work->hash_cap-1);

      while (work->hash[slot] &&
             !mapq_same_geometry(work->paths+i,work->paths+work->hash[slot]-1))
        slot = (slot+1) & (work->hash_cap-1);
      if (!work->hash[slot])
        work->hash[slot] = i+1;
      work->geometry[i] = work->hash[slot]-1;
    }
  work->trace_n = 0;
  for (int i = 0; i < n; ++i)
    work->trace_offset[i] = SIZE_MAX;
  if (!mapq_summaries(work,shadows,nshadows,joined))
    return;

  memset(work->primary,0,(size_t) n*sizeof(int));
  memset(work->component_primary,0,(size_t) n*sizeof(int));
  memset(work->state,0,(size_t) n*sizeof(int));
  memset(work->hash,0,work->hash_cap*sizeof(int));
  for (int i = 0; i < n; ++i)
    { int           index = work->ranks[i], j;
      Mapq_Path    *path = work->paths+index;
      Mapq_Summary *summary = work->summary+index;
      Mapq_Evidence evidence;
      int           strongest, end;

      if (!path->root || !isfinite(path->score) || path->score < 0. || work->state[index])
        continue;
      /* Confidence uses a primary-only traversal, independent of whether the
         caller emits secondary chains.  An already consumed tail is omitted
         when testing this root's new query interval. */
      summary->primary_qb = path->qb;
      summary->primary_qe = path->qe;
      evidence = path->evidence;
      strongest = end = index;
      for (int p = index; p >= 0 && !work->state[p]; p = work->paths[p].next)
        { Mapq_Path *piece = work->paths+p;

          end = p;
          if (summary->primary_qb > piece->qb)
            summary->primary_qb = piece->qb;
          if (summary->primary_qe < piece->qe)
            summary->primary_qe = piece->qe;
          Mapq_Combine_Evidence(&evidence,&piece->evidence);
          if (evidence.chain_score == piece->evidence.chain_score &&
              evidence.seed_count == piece->evidence.seed_count &&
              evidence.chain_id == piece->evidence.chain_id)
            strongest = p;
        }
      for (j = 0; j < nparts; ++j)
        { Mapq_Summary *previous = work->summary+work->stack[j];

          if ((int64_t) summary->primary_qb < (int64_t) previous->primary_qe-100 &&
              (int64_t) summary->primary_qe > (int64_t) previous->primary_qb+100)
            break;
        }
      if (j == nparts)
        { size_t slot = mapq_source_slot(work,evidence.chain_id);

          work->stack[nparts++] = index;
          work->primary[index] = 1;
          for (int p = index; ; p = work->paths[p].next)
            { work->state[p] = 1;
              work->component_primary[work->geometry[p]] = 1;
              if (p == end)
                break;
            }
          if (evidence.chain_score > 0 &&
              (!work->hash[slot] ||
               work->paths[work->hash[slot]-1].evidence.chain_score < evidence.chain_score))
            { if (work->hash[slot])
                sum -= work->paths[work->hash[slot]-1].evidence.chain_score;
              sum += evidence.chain_score;
              work->hash[slot] = strongest+1;
            }
        }
    }
  unique_ratio = sum > 0. ? sum/(sum+repeat_bases) : 0.;
  work->unique_ratio = unique_ratio;

  memset(work->seen,0,(size_t) n*sizeof(int));
  for (int i = 0; i < n; ++i)
    if (work->primary[i])
      { Mapq_Summary    *summary = work->summary+i;
        Mapq_Input       input;
        Mapq_Competition competition = {0};
        int              component = work->suffix[i], lost_score, quality, own_score;

        input = (Mapq_Input) {work->paths[i].score,summary->identity,
                            summary->evidence.chain_score,summary->evidence.seed_count,
                            0,summary->incomplete};
        for (int j = 0; j < n; ++j)
          if (work->paths[j].root && work->suffix[j] != component &&
              Mapq_Query_Overlap(summary->qb,summary->qe,work->summary[j].qb,work->summary[j].qe))
            { if (!isfinite(work->paths[j].score))
                input.incomplete = 1;
              else
                mapq_competitor(work,i+1,work->suffix[j],work->paths[j].score,
                               &work->summary[j].evidence,&input,&competition);
            }
        competition.drop_score = mapq_dropped(summary->qb,summary->qe,summary->own_shadows,
                                             shadows,nshadows);
        lost_score = mapq_overflow_score(summary->qb,summary->qe,overflow);
        if (competition.drop_score < lost_score)
          competition.drop_score = lost_score;
        quality = Mapq_Calculate(&input,&competition,min_chain_score,unique_ratio);
        lost_score = mapq_joined_score(summary->qb,summary->qe,summary->other_joined,joined);
        own_score = mapq_own_joined_score(summary->qb,summary->qe,summary->rb,summary->re,
                                          summary->reference,summary->reverse,
                                          summary->own_joined,joined);
        work->group_quality[i] = mapq_joined_cap(quality,&input,unique_ratio,own_score,lost_score);
      }

  memset(work->seen,0,(size_t) n*sizeof(int));
  for (int i = 0; i < n; ++i)
    if (work->component_primary[work->geometry[i]])
      { Mapq_Path       *path = work->paths+i;
        Mapq_Input       input;
        Mapq_Competition competition = {0};
        int              component = work->geometry[i], lost_score, quality, own_score;
        uint32_t         own_shadows = 0;

        input = (Mapq_Input) {path->local_score,path->identity,
                            path->evidence.chain_score,path->evidence.seed_count,
                            0,path->incomplete || work->summary[i].incomplete};
        for (int j = 0; j < n; ++j)
          if (work->geometry[j] != component &&
              Mapq_Query_Overlap(path->qb,path->qe,work->paths[j].qb,work->paths[j].qe))
            { if (!isfinite(work->paths[j].local_score))
                input.incomplete = 1;
              else
                mapq_competitor(work,i+1,work->geometry[j],work->paths[j].local_score,
                               &work->paths[j].evidence,&input,&competition);
            }
        for (int j = 0; j < nshadows; ++j)
          if (shadows[j].chain_id == path->evidence.chain_id)
            own_shadows |= UINT32_C(1) << j;
        competition.drop_score = mapq_dropped(path->qb,path->qe,own_shadows,shadows,nshadows);
        lost_score = mapq_overflow_score(path->qb,path->qe,overflow);
        if (competition.drop_score < lost_score)
          competition.drop_score = lost_score;
        quality = Mapq_Calculate(&input,&competition,min_chain_score,unique_ratio);
        lost_score = mapq_joined_score(path->qb,path->qe,work->summary[i].local_other_joined,joined);
        own_score = mapq_own_joined_score(path->qb,path->qe,path->rb,path->re,
                                          path->reference,path->reverse,
                                          work->summary[i].local_own_joined,joined);
        work->local_quality[i] = mapq_joined_cap(quality,&input,unique_ratio,own_score,lost_score);
      }
  work->assigned = 1;
}

int Mapq_Quality(const Mapq_Work *work, int root, int segment)
{ int group, local;

  if (!work || root < 0 || root >= work->n || segment < 0 || segment >= work->n)
    return 0;
  group = work->group_quality[root];
  local = work->local_quality[segment];
  return group < local ? group : local;
}

size_t Mapq_Work_Bytes(const Mapq_Work *work)
{ if (!work)
    return 0;
  return sizeof(*work) + (size_t) work->cap*(sizeof(Mapq_Path)+sizeof(Mapq_Summary)+14*sizeof(int)+sizeof(size_t)
                        +2*sizeof(Mapq_Relation)+2*sizeof(Mapq_Support))
                      + work->hash_cap*sizeof(int) + work->trace_cap*sizeof(int);
}

void Mapq_Free_Work(Mapq_Work *work)
{ if (!work)
    return;
  free(work->paths);
  free(work->summary);
  free(work->ranks);
  free(work->geometry);
  free(work->suffix);
  free(work->state);
  free(work->stack);
  free(work->primary);
  free(work->component_primary);
  free(work->seen);
  free(work->group_quality);
  free(work->local_quality);
  free(work->joined_local);
  free(work->joined_group);
  free(work->local_support);
  free(work->group_support);
  free(work->local_stamp);
  free(work->group_stamp);
  free(work->support_targets);
  free(work->trace_ref);
  free(work->trace_offset);
  free(work->hash);
  free(work);
}

typedef struct
  { Mapq_Shadow fragment;
    Mapq_Join_Geometry geometry;
    int64_t     band;
    int         score, seeds, first, predecessor, parts, used;
  } Mapq_Fragment_Node;

struct Mapq_Fragment_Pass
  { size_t offset;
    int    count, kmer, max_gap, drift_penalty_x100;
  };

struct Mapq_Fragment_Work
  { Mapq_Shadow         *input;
    Mapq_Fragment_Node  *nodes;
    int                  n, cap, node_cap, timing;
    Mapq_Fragment_Stats  stats;
  };

static double mapq_thread_seconds(void)
{ struct timespec now;

  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID,&now))
    return 0.;
  return now.tv_sec + now.tv_nsec*1e-9;
}

Mapq_Fragment_Work *Mapq_New_Fragment_Work(int timing)
{ Mapq_Fragment_Work *work = calloc(1,sizeof(*work));

  if (work)
    work->timing = timing != 0;
  return work;
}

void Mapq_Fragment_Reset(Mapq_Fragment_Work *work)
{ if (work)
    work->n = 0;
}

int Mapq_Fragment_Add(Mapq_Fragment_Work *work, const Mapq_Shadow *fragment)
{ if (!work)
    { errno = EINVAL;
      return 0;
    }
  ++work->stats.input;
  if (!fragment || fragment->qb < 0 || fragment->qe <= fragment->qb ||
      fragment->rb < 0 || fragment->re <= fragment->rb || fragment->reference < 0 ||
      (fragment->reverse != 0 && fragment->reverse != 1) || fragment->chain_score <= 0 ||
      fragment->seed_count <= 0 || !fragment->chain_id || work->n == MAPQ_FRAGMENT_MAX)
    { ++work->stats.skipped;
      return 1;
    }
  if (work->n == work->cap)
    { int          cap = work->cap ? work->cap*2 : 64;
      Mapq_Shadow *input = realloc(work->input,(size_t) cap*sizeof(*input));

      if (!input)
        return 0;
      work->input = input;
      work->cap = cap;
    }
  work->input[work->n++] = *fragment;
  return 1;
}

static int mapq_compare_int(int a, int b)
{ return (a > b) - (a < b);
}

static int mapq_fragment_geometry(const void *left, const void *right)
{ const Mapq_Shadow *a = left, *b = right;
  int cmp;

  if ((cmp = mapq_compare_int(a->reference,b->reference))) return cmp;
  if ((cmp = mapq_compare_int(a->reverse,b->reverse))) return cmp;
  if ((cmp = mapq_compare_int(a->qb,b->qb))) return cmp;
  if ((cmp = mapq_compare_int(a->qe,b->qe))) return cmp;
  if ((cmp = mapq_compare_int(a->rb,b->rb))) return cmp;
  if ((cmp = mapq_compare_int(a->re,b->re))) return cmp;
  if ((cmp = mapq_compare_int(b->chain_score,a->chain_score))) return cmp;
  if ((cmp = mapq_compare_int(b->seed_count,a->seed_count))) return cmp;
  return (a->chain_id > b->chain_id) - (a->chain_id < b->chain_id);
}

static int mapq_fragment_id(const void *left, const void *right)
{ const Mapq_Shadow *a = left, *b = right;
  int cmp = (a->chain_id > b->chain_id) - (a->chain_id < b->chain_id);

  if (cmp) return cmp;
  if ((cmp = mapq_compare_int(b->chain_score,a->chain_score))) return cmp;
  if ((cmp = mapq_compare_int(b->seed_count,a->seed_count))) return cmp;
  return mapq_fragment_geometry(a,b);
}

static int mapq_fragment_same(const Mapq_Shadow *a, const Mapq_Shadow *b)
{ return a->reference == b->reference && a->reverse == b->reverse &&
         a->qb == b->qb && a->qe == b->qe && a->rb == b->rb && a->re == b->re;
}

static int mapq_node_compare(const void *left, const void *right)
{ const Mapq_Fragment_Node *a = left, *b = right;
  int cmp;

  if ((cmp = mapq_compare_int(a->fragment.reference,b->fragment.reference))) return cmp;
  if ((cmp = mapq_compare_int(a->fragment.reverse,b->fragment.reverse))) return cmp;
  if (a->band != b->band) return (a->band > b->band) - (a->band < b->band);
  return mapq_fragment_geometry(&a->fragment,&b->fragment);
}

static int mapq_join_compare(const void *left, const void *right)
{ const Mapq_Shadow *a = left, *b = right;
  int cmp;
  int64_t a_span = (int64_t) a->qe-a->qb + (int64_t) a->re-a->rb;
  int64_t b_span = (int64_t) b->qe-b->qb + (int64_t) b->re-b->rb;

  if ((cmp = mapq_compare_int(b->chain_score,a->chain_score))) return cmp;
  if (a_span != b_span) return (b_span > a_span) - (b_span < a_span);
  if ((cmp = mapq_compare_int(b->seed_count,a->seed_count))) return cmp;
  return mapq_fragment_geometry(a,b);
}

/* Compress only actual boundary points.  First/last always survive; interior
   samples are spread across query coordinates after the best predecessor has
   been chosen.  This work is done once per node, not once per DP comparison. */

static void mapq_join_points(Mapq_Join_Geometry *geometry,
                              const Mapq_Join_Geometry *previous, const Mapq_Shadow *fragment)
{ Mapq_Point points[MAPQ_JOIN_POINTS+2];
  int count = previous ? previous->n : 0, kept = 0;

  if (previous)
    memcpy(points,previous->points,(size_t) count*sizeof(*points));
  points[count++] = (Mapq_Point) {fragment->qb,fragment->rb};
  points[count++] = (Mapq_Point) {fragment->qe,fragment->re};
  for (int i = 0; i < count; ++i)
    if (!kept || points[i].query != points[kept-1].query ||
                 points[i].reference != points[kept-1].reference)
      points[kept++] = points[i];
  if (kept <= MAPQ_JOIN_POINTS)
    { geometry->n = kept;
      memcpy(geometry->points,points,(size_t) kept*sizeof(*points));
    }
  else
    { int last = 0;

      geometry->n = MAPQ_JOIN_POINTS;
      geometry->points[0] = points[0];
      for (int i = 1; i < MAPQ_JOIN_POINTS-1; ++i)
        { int best = last+1;
          int limit = kept-(MAPQ_JOIN_POINTS-i);
          int64_t target = (int64_t) points[0].query*(MAPQ_JOIN_POINTS-1-i)
                         + (int64_t) points[kept-1].query*i;
          int64_t distance = llabs((int64_t) points[best].query*(MAPQ_JOIN_POINTS-1)-target);

          for (int j = best+1; j <= limit; ++j)
            { int64_t delta = llabs((int64_t) points[j].query*(MAPQ_JOIN_POINTS-1)-target);

              if (delta < distance)
                { best = j;
                  distance = delta;
                }
            }
          geometry->points[i] = points[best];
          last = best;
        }
      geometry->points[MAPQ_JOIN_POINTS-1] = points[kept-1];
    }
}

/* A shared bounding box alone cannot identify a repeat copy.  Require two
   exact sampled correspondences spanning >=80% of the shorter query range.
   This collapses the same chain found in both shifted bands without merging
   nearby parallel placements.  Unproven duplicates can occupy separate slots. */

static int mapq_join_same(const Mapq_Shadow *a, const Mapq_Join_Geometry *ag,
                            const Mapq_Shadow *b, const Mapq_Join_Geometry *bg)
{ int first = INT_MAX, last = -1;

  if (a->reference != b->reference || a->reverse != b->reverse ||
      !mapq_reciprocal_overlap(a->qb,a->qe,b->qb,b->qe) ||
      !mapq_reciprocal_overlap(a->rb,a->re,b->rb,b->re))
    return 0;
  for (int i = 0; i < ag->n; ++i)
    for (int j = 0; j < bg->n; ++j)
      if (ag->points[i].query == bg->points[j].query &&
          ag->points[i].reference == bg->points[j].reference)
        { int q = ag->points[i].query;

          if (q < first) first = q;
          if (q > last) last = q;
        }
  return last > first && (int64_t) (last-first)*5 >=
         (int64_t) ((a->qe-a->qb < b->qe-b->qb) ? a->qe-a->qb : b->qe-b->qb)*4;
}

static int mapq_join_merge(Mapq_Join_Catalog *out, const Mapq_Shadow *joined,
                            const Mapq_Join_Geometry *geometry)
{ int at = -1;

  for (int i = 0; i < out->n; ++i)
    if (mapq_join_same(out->items+i,out->geometry+i,joined,geometry))
      { if (mapq_join_compare(joined,out->items+i) >= 0)
          return 1;
        at = i;
        break;
      }
  if (at < 0)
    { if (out->n == MAPQ_JOIN_MAX)
        { if (mapq_join_compare(joined,out->items+out->n-1) >= 0)
            return 1;
          at = out->n-1;
        }
      else
        { if (out->n == out->cap)
            { int cap = out->cap ? out->cap*2 : 4;
              Mapq_Shadow *items;
              Mapq_Join_Geometry *points;

              if (cap > MAPQ_JOIN_MAX)
                cap = MAPQ_JOIN_MAX;
              items = realloc(out->items,(size_t) cap*sizeof(*items));
              if (!items)
                return 0;
              out->items = items;
              points = realloc(out->geometry,(size_t) cap*sizeof(*points));
              if (!points)
                return 0;
              out->geometry = points;
              out->cap = cap;
            }
          at = out->n++;
        }
    }
  /* Every insertion/replacement improves this slot, so it moves only left. */
  while (at > 0 && mapq_join_compare(joined,out->items+at-1) < 0)
    { out->items[at] = out->items[at-1];
      out->geometry[at] = out->geometry[at-1];
      --at;
    }
  out->items[at] = *joined;
  out->geometry[at] = *geometry;
  return 1;
}

static int64_t mapq_diagonal_band(int64_t diagonal, int64_t offset, int64_t width)
{ int64_t value = diagonal-offset;
  int64_t band = value/width;

  if (value < 0 && value%width)
    --band;
  return band;
}

static void mapq_fragment_unique(Mapq_Fragment_Work *work)
{
  /* A raw chain or an identical geometric fragment can contribute only once.
     Keep its strongest record before either band is evaluated. */
  qsort(work->input,work->n,sizeof(*work->input),mapq_fragment_id);
  { int kept = 0;

    for (int i = 0; i < work->n; ++i)
      if (!kept || work->input[i].chain_id != work->input[kept-1].chain_id)
        work->input[kept++] = work->input[i];
    work->n = kept;
  }
  qsort(work->input,work->n,sizeof(*work->input),mapq_fragment_geometry);
  { int kept = 0;

    for (int i = 0; i < work->n; ++i)
      if (!kept || !mapq_fragment_same(work->input+i,work->input+kept-1))
        work->input[kept++] = work->input[i];
    work->n = kept;
  }

}

int Mapq_Fragment_Save(Mapq_Fragment_Work *work, int kmer, int max_gap,
                        int drift_penalty_x100, Mapq_Fragment_Store *store)
{ size_t need;

  if (!work || !store || kmer <= 0 || max_gap < 0 || drift_penalty_x100 < 0 ||
      store->n > store->cap || store->npasses < 0 || store->pass_cap < store->npasses ||
      (store->cap && !store->items) || (store->pass_cap && !store->passes))
    { errno = EINVAL;
      return 0;
    }
  if (work->n < 2)
    return 1;
  mapq_fragment_unique(work);
  if (work->n < 2)
    return 1;
  if (store->n > SIZE_MAX-(size_t) work->n || store->npasses == INT_MAX)
    { errno = ENOMEM;
      return 0;
    }
  need = store->n+(size_t) work->n;
  if (need > store->cap)
    { size_t cap = need+need/2+4;
      Mapq_Shadow *items;

      if (cap < need || cap > SIZE_MAX/sizeof(*items))
        cap = need;
      if (cap > SIZE_MAX/sizeof(*items))
        { errno = ENOMEM;
          return 0;
        }
      items = realloc(store->items,cap*sizeof(*items));
      if (!items)
        return 0;
      store->items = items;
      store->cap = cap;
    }
  if (store->npasses == store->pass_cap)
    { int cap = store->npasses > (INT_MAX-4)/2 ? INT_MAX : store->npasses*2+4;
      Mapq_Fragment_Pass *passes;

      if ((size_t) cap > SIZE_MAX/sizeof(*passes))
        { errno = ENOMEM;
          return 0;
        }
      passes = realloc(store->passes,(size_t) cap*sizeof(*passes));
      if (!passes)
        return 0;
      store->passes = passes;
      store->pass_cap = cap;
    }
  memcpy(store->items+store->n,work->input,(size_t) work->n*sizeof(*work->input));
  store->passes[store->npasses++] = (Mapq_Fragment_Pass)
      {store->n,work->n,kmer,max_gap,drift_penalty_x100};
  store->n = need;
  return 1;
}

void Mapq_Free_Fragment_Store(Mapq_Fragment_Store *store)
{ if (!store)
    return;
  free(store->items);
  free(store->passes);
  memset(store,0,sizeof(*store));
}

size_t Mapq_Fragment_Store_Bytes(const Mapq_Fragment_Store *store)
{ return store ? store->cap*sizeof(*store->items) +
                 (size_t) store->pass_cap*sizeof(*store->passes) : 0;
}

typedef int (*Mapq_Fragment_Emit)(void *, const Mapq_Shadow *, const Mapq_Join_Geometry *);

static int mapq_fragment_evaluate(Mapq_Fragment_Work *work, int kmer, int max_gap,
                                   int drift_penalty_x100, int terminal_only,
                                   const Mapq_Shadow *input, int n,
                                   Mapq_Fragment_Emit emit, void *context)
{ int    result = 1;
  double started;

  if (!work || !emit || kmer <= 0 || max_gap < 0 || drift_penalty_x100 < 0)
    { errno = EINVAL;
      return 0;
    }
  if (!input)
    n = work->n;
  if (n < 2)
    return 1;
  started = work->timing ? mapq_thread_seconds() : 0.;
  if (n > work->node_cap)
    { Mapq_Fragment_Node *nodes = realloc(work->nodes,(size_t) n*sizeof(*nodes));

      if (!nodes)
        { result = 0;
          goto done;
        }
      work->nodes = nodes;
      work->node_cap = n;
    }

  if (!input)
    { mapq_fragment_unique(work);
      input = work->input;
      n = work->n;
    }

  for (int phase = 0; phase < 2; ++phase)
    { int     count = 0, first = 0;
      int64_t width = (int64_t) kmer*2, offset = (int64_t) phase*kmer;

      for (int i = 0; i < n; ++i)
        { const Mapq_Shadow *fragment = input+i;
          int64_t begin = mapq_diagonal_band((int64_t) fragment->rb-fragment->qb,offset,width);
          int64_t end   = mapq_diagonal_band((int64_t) fragment->re-fragment->qe,offset,width);

          if (begin == end)
            { Mapq_Fragment_Node *node = work->nodes+count++;

              node->fragment = *fragment;
              node->band = begin;
            }
        }
      qsort(work->nodes,count,sizeof(*work->nodes),mapq_node_compare);
      for (int i = 0; i < count; ++i)
        { Mapq_Fragment_Node *node = work->nodes+i;
          Mapq_Shadow *fragment = &node->fragment;

          if (!i || fragment->reference != work->nodes[i-1].fragment.reference ||
              fragment->reverse != work->nodes[i-1].fragment.reverse ||
              node->band != work->nodes[i-1].band)
            first = i;
          node->score = fragment->chain_score;
          node->seeds = fragment->seed_count;
          node->first = i;
          node->predecessor = -1;
          node->parts = 1;
          node->used = 0;
          for (int j = i-1; j >= first && i-j <= 64; --j)
            { Mapq_Fragment_Node *previous = work->nodes+j;
              Mapq_Shadow *last = &previous->fragment;
              int64_t dq, dr, drift, score, bonus, seeds;

              if (fragment->qb < last->qe || fragment->rb < last->re)
                continue;
              dq = (int64_t) fragment->qb+kmer-last->qe;
              dr = (int64_t) fragment->rb+kmer-last->re;
              if (dq < 0 || dr < 0 || dq > max_gap || dr > max_gap)
                continue;
              drift = dq > dr ? dq-dr : dr-dq;
              bonus = dq < dr ? dq : dr;
              if (bonus > kmer)
                bonus = kmer;
              score = (int64_t) previous->score+fragment->chain_score+bonus-kmer
                    - (int64_t) drift_penalty_x100*drift/100;
              if (score > INT_MAX)
                score = INT_MAX;
              if (score < node->score ||
                  (score == node->score && previous->parts+1 <= node->parts))
                continue;
              seeds = (int64_t) previous->seeds+fragment->seed_count;
              node->score = (int) score;
              node->seeds = seeds > INT_MAX ? INT_MAX : (int) seeds;
              node->first = previous->first;
              node->predecessor = j;
              node->parts = previous->parts+1;
            }
          mapq_join_points(&node->geometry,
                           node->predecessor >= 0 ? &work->nodes[node->predecessor].geometry : NULL,
                           fragment);
        }
      for (int i = 0; i < count; ++i)
        if (work->nodes[i].predecessor >= 0)
          work->nodes[work->nodes[i].predecessor].used = 1;
      for (int i = 0; i < count; ++i)
        { Mapq_Fragment_Node *node = work->nodes+i;

          if (node->parts >= 2 && (!terminal_only || !node->used))
            { Mapq_Shadow joined = node->fragment;
              const Mapq_Shadow *begin = &work->nodes[node->first].fragment;

              joined.qb = begin->qb;
              joined.rb = begin->rb;
              joined.chain_score = node->score;
              joined.seed_count = node->seeds;
              joined.chain_id = 0;
              ++work->stats.joins;
              if (!emit(context,&joined,&node->geometry))
                { result = 0;
                  goto done;
                }
            }
        }
    }
done:
  if (work->timing)
    work->stats.cpu_seconds += mapq_thread_seconds()-started;
  return result;
}

/* Compatibility materialization uses the same DP engine; production streams
   every prefix directly to query-specific support accumulators instead. */

static int mapq_emit_catalog(void *context, const Mapq_Shadow *item,
                              const Mapq_Join_Geometry *geometry)
{ return mapq_join_merge(context,item,geometry);
}

int Mapq_Fragment_Build(Mapq_Fragment_Work *work, int kmer, int max_gap,
                         int drift_penalty_x100, Mapq_Join_Catalog *out)
{ if (!out || out->n < 0 || out->n > MAPQ_JOIN_MAX || out->cap < out->n ||
      out->cap > MAPQ_JOIN_MAX || (out->cap && (!out->items || !out->geometry)))
    { errno = EINVAL;
      return 0;
    }
  return mapq_fragment_evaluate(work,kmer,max_gap,drift_penalty_x100,1,NULL,0,
                                 mapq_emit_catalog,out);
}

static Mapq_Relation mapq_retained_local(Mapq_Work *work, int at,
                                          const Mapq_Shadow *item,
                                          const Mapq_Join_Geometry *geometry)
{ if (work->local_stamp[at] != work->joined_stamp)
    { work->joined_local[at] = mapq_joined_relation(work,at,item,geometry);
      work->local_stamp[at] = work->joined_stamp;
    }
  return work->joined_local[at];
}

/* Cache suffix unions within this streamed candidate.  Shared suffixes are
   visited once even when many roots converge on them.  An anchor explained by
   any suffix path must not become DIFFERENT because another path disagrees. */

static Mapq_Relation mapq_retained_group(Mapq_Work *work, int root,
                                          const Mapq_Shadow *item,
                                          const Mapq_Join_Geometry *geometry)
{ int at = root, depth = 0;

  while (at >= 0 && work->group_stamp[at] != work->joined_stamp)
    { work->stack[depth++] = at;
      at = work->paths[at].next;
    }
  while (depth)
    { Mapq_Relation relation;
      int next;

      at = work->stack[--depth];
      relation = mapq_retained_local(work,at,item,geometry);
      next = work->paths[at].next;
      if (next >= 0)
        { const Mapq_Relation *tail = work->joined_group+next;

          relation.covered  |= tail->covered;
          relation.possible |= tail->possible;
          relation.narrow   |= tail->narrow;
          relation.nearby   |= tail->nearby;
        }
      work->joined_group[at] = relation;
      work->group_stamp[at] = work->joined_stamp;
    }
  return work->joined_group[root];
}

static int mapq_emit_support(void *context, const Mapq_Shadow *item,
                              const Mapq_Join_Geometry *geometry)
{ Mapq_Work *work = context;

  if (work->joined_stamp == INT_MAX)
    { memset(work->local_stamp,0,(size_t) work->n*sizeof(int));
      memset(work->group_stamp,0,(size_t) work->n*sizeof(int));
      work->joined_stamp = 0;
    }
  ++work->joined_stamp;
  for (size_t k = 0; k < work->ntargets; ++k)
    { int target = work->support_targets[k];
      int group = target < 0;
      int at = group ? ~target : target;
      Mapq_Support *support = group ? work->group_support+at : work->local_support+at;
      Mapq_Relation relation;
      int qb, qe, rb, re, reference, reverse, original;

      if (group)
        { const Mapq_Summary *summary = work->summary+at;

          qb = summary->qb; qe = summary->qe; rb = summary->rb; re = summary->re;
          reference = summary->reference; reverse = summary->reverse;
          original = summary->evidence.chain_score;
        }
      else
        { const Mapq_Path *path = work->paths+at;

          qb = path->qb; qe = path->qe; rb = path->rb; re = path->re;
          reference = path->reference; reverse = path->reverse;
          original = path->evidence.chain_score;
        }
      /* Eligibility comes before retaining a score.  Evidence from an
         unrelated query range never displaces an applicable candidate. */
      if (!mapq_reciprocal_overlap(qb,qe,item->qb,item->qe) ||
          !Mapq_Query_Overlap(qb,qe,item->qb,item->qe))
        continue;
      if (item->chain_score <= support->other &&
          (item->chain_score <= support->own || item->chain_score <= original))
        continue;
      /* A homogeneous group's union cannot overlap this reference span if
         its bounding box does not; no trace inspection is needed here. */
      if (reference >= 0 && (reference != item->reference || reverse != item->reverse ||
                            rb >= item->re || item->rb >= re))
        { if (support->other < item->chain_score)
            support->other = item->chain_score;
          continue;
        }
      relation = group ? mapq_retained_group(work,at,item,geometry) :
                         mapq_retained_local(work,at,item,geometry);
      if (!relation.nearby || (relation.covered & ~relation.possible))
        { if (support->other < item->chain_score)
            support->other = item->chain_score;
        }
      else if (relation.narrow == (1U << geometry->n)-1 &&
               reference >= 0 && reference == item->reference && reverse == item->reverse &&
               mapq_reciprocal_overlap(rb,re,item->rb,item->re))
        { if (support->own < item->chain_score)
            support->own = item->chain_score;
        }
      /* UNKNOWN supplies neither a new rival cap nor own support. */
    }
  return 1;
}

static int mapq_valid_store(const Mapq_Fragment_Store *store)
{ size_t offset = 0;

  if (!store)
    return 1;
  if (store->n > store->cap || store->npasses < 0 || store->pass_cap < store->npasses ||
      (store->cap && !store->items) || (store->pass_cap && !store->passes))
    return 0;
  for (int p = 0; p < store->npasses; ++p)
    { const Mapq_Fragment_Pass *pass = store->passes+p;

      if (pass->offset != offset || pass->count <= 0 || pass->count > MAPQ_FRAGMENT_MAX ||
          pass->kmer <= 0 || pass->max_gap < 0 || pass->drift_penalty_x100 < 0 ||
          offset > store->n || (size_t) pass->count > store->n-offset)
        return 0;
      for (int i = 0; i < pass->count; ++i)
        { const Mapq_Shadow *fragment = store->items+offset+i;

          if (fragment->qb < 0 || fragment->qe <= fragment->qb ||
              fragment->rb < 0 || fragment->re <= fragment->rb || fragment->reference < 0 ||
              (fragment->reverse != 0 && fragment->reverse != 1) || fragment->chain_score <= 0 ||
              fragment->seed_count <= 0 || !fragment->chain_id)
            return 0;
        }
      offset += pass->count;
    }
  return offset == store->n;
}

void Mapq_Assign_Retained(Mapq_Work *work, int n, const Mapq_Shadow *shadows,
                           int nshadows, int min_chain_score, int repeat_bases,
                           const Mapq_Overflow *overflow, const Mapq_Fragment_Store *store,
                           Mapq_Fragment_Work *scratch)
{ Mapq_Assign(work,n,shadows,nshadows,min_chain_score,repeat_bases,overflow,NULL);
  if (!work || !work->assigned || !work->n)
    return;
  if (!mapq_valid_store(store) || (store && store->n && !scratch))
    { memset(work->group_quality,0,(size_t) work->n*sizeof(int));
      memset(work->local_quality,0,(size_t) work->n*sizeof(int));
      work->assigned = 0;
      return;
    }
  if (!store || !store->n)
    return;
  memset(work->group_support,0,(size_t) n*sizeof(*work->group_support));
  memset(work->local_support,0,(size_t) n*sizeof(*work->local_support));
  memset(work->group_stamp,0,(size_t) n*sizeof(int));
  memset(work->local_stamp,0,(size_t) n*sizeof(int));
  work->joined_stamp = 0;
  work->ntargets = 0;
  for (int i = 0; i < n; ++i)
    { if (work->group_quality[i]) work->support_targets[work->ntargets++] = ~i;
      if (work->local_quality[i]) work->support_targets[work->ntargets++] = i;
    }
  if (!work->ntargets)
    return;
  for (int p = 0; p < store->npasses; ++p)
    { const Mapq_Fragment_Pass *pass = store->passes+p;

      if (pass->count < 2)
        continue;
      if (!mapq_fragment_evaluate(scratch,pass->kmer,pass->max_gap,pass->drift_penalty_x100,
                                   0,store->items+pass->offset,pass->count,mapq_emit_support,work))
        { perror("MAPQ retained-fragment evaluation");
          exit(1);
        }
    }
  for (size_t k = 0; k < work->ntargets; ++k)
    { int target = work->support_targets[k];
      int group = target < 0;
      int at = group ? ~target : target;
      const Mapq_Support *support = group ? work->group_support+at : work->local_support+at;
      Mapq_Input input;
      int *quality;

      if (group)
        { const Mapq_Summary *summary = work->summary+at;

          input = (Mapq_Input) {work->paths[at].score,summary->identity,
                               summary->evidence.chain_score,summary->evidence.seed_count,
                               0,summary->incomplete};
          quality = work->group_quality+at;
        }
      else
        { const Mapq_Path *path = work->paths+at;

          input = (Mapq_Input) {path->local_score,path->identity,
                               path->evidence.chain_score,path->evidence.seed_count,
                               0,path->incomplete || work->summary[at].incomplete};
          quality = work->local_quality+at;
        }
      *quality = mapq_joined_cap(*quality,&input,work->unique_ratio,support->own,support->other);
    }
}

Mapq_Fragment_Stats Mapq_Fragment_Statistics(const Mapq_Fragment_Work *work)
{ Mapq_Fragment_Stats stats = {0};

  if (work)
    { stats = work->stats;
      stats.allocated_bytes = sizeof(*work) + (size_t) work->cap*sizeof(*work->input)
                            + (size_t) work->node_cap*sizeof(*work->nodes);
    }
  return stats;
}

void Mapq_Free_Fragment_Work(Mapq_Fragment_Work *work)
{ if (!work)
    return;
  free(work->input);
  free(work->nodes);
  free(work);
}

void Mapq_Free_Join_Catalog(Mapq_Join_Catalog *catalog)
{ if (!catalog)
    return;
  free(catalog->items);
  free(catalog->geometry);
  memset(catalog,0,sizeof(*catalog));
}

size_t Mapq_Join_Bytes(const Mapq_Join_Catalog *catalog)
{ return catalog ? (size_t) catalog->cap*(sizeof(*catalog->items)+sizeof(*catalog->geometry)) : 0;
}
