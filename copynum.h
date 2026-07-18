/*******************************************************************************************
 *
 *  copynum.h — copy-number-aware anchoring support (shared by minimaxmer & himapper).
 *
 *  minimax-mer already uses the genome copy number (#G, from FastK) when SELECTING
 *  reference seeds, but until now that information was discarded afterwards: the
 *  .kpos stored only (code, rpos, rid) and himapper weighted every seed hit equally
 *  during chaining / alignment selection.  In homologous / repetitive regions this
 *  lets many shared (copy>=2) seeds drown out the few distinguishing (copy-1) seeds,
 *  so a read can be anchored to the wrong haplotype copy.
 *
 *  To carry the copy number downstream WITHOUT growing any hot data structure
 *  (KmerPos and SeedPair stay 16 bytes — minimax-mer's low-memory footprint is a
 *  headline feature), we pack a 3-bit genome copy *tier* into the spare top bits of
 *  two int fields that are guaranteed to have room:
 *
 *     - reference contig id (rid):  KmerPos.read (uint32) / SeedPair.bread.  Real
 *       assemblies have far fewer than 2^29 contigs, so bits 29..31 are free.
 *     - read position (apos):       SeedPair.apos.  Read length < 2^24 (and, for the
 *       radix sort key to stay below byte 7, < 2^24 = 16Mbp — guarded at runtime).
 *
 *  The low 29 bits hold the real value; bits 29..31 hold the 3-bit tier (0..7).  Tier 0
 *  means "0 copies / unknown": it is emitted only for a genome count of 0 (which cannot
 *  happen for a seed that actually matches the reference, so it flags an anomaly) or by
 *  untagged / legacy indices, and contributes no copy evidence — old indices keep working.
 *
 ********************************************************************************************/

#ifndef _COPYNUM_H
#define _COPYNUM_H

#include <stdint.h>

#define CN_TIER_SHIFT 29
#define CN_LOW_MASK   ((1u << CN_TIER_SHIFT) - 1u)        //  0x1FFFFFFF
#define CN_MAXVAL     (1u << CN_TIER_SHIFT)               //  2^29: values must be below this
#define CN_MAXLEN     (1 << 24)                           //  2^24: read length must be below this
                                                          //  (so the radix sort key stays in bytes 4..6)

#define CN_TIER(x)    ((unsigned)(x) >> CN_TIER_SHIFT)            //  extract tier (0..7)
#define CN_VAL(x)     ((int)((unsigned)(x) & CN_LOW_MASK))       //  extract real value
#define CN_PACK(v, t) ((uint32_t)(((unsigned)(v) & CN_LOW_MASK) | ((unsigned)(t) << CN_TIER_SHIFT)))

// Quantize a genome k-mer copy number to a tier.
static inline int cn_tier_of(uint32_t g) {
  if (g == 0)   return 0;   // 0 copies: cannot happen for a seed that matches the ref => a bug
  if (g == 1)   return 1;   // haploid, unique
  if (g == 2)   return 2;   // diploid
  if (g == 3)   return 3;
  if (g == 4)   return 4;
  if (g <= 10)  return 5;   // copy 5-10
  if (g <= 100) return 6;   // copy 11-100
  return 7;                 // copy >100
}

#endif  //  _COPYNUM_H
