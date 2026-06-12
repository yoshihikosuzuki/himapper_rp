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
 *     - reference contig id (rid):  KmerPos.read / SeedPair.bread.  Real assemblies
 *       have far fewer than 2^29 contigs, so bits 29..31 are free.
 *     - read position (apos):       SeedPair.apos.  Read length < 2^29 (and, for the
 *       radix sort key to stay below byte 7, < 2^24 = 16Mbp — guarded at runtime).
 *
 *  The low 29 bits hold the real value; bits 29..31 hold the tier.  Tier 0 means
 *  "unknown" (e.g. a .kpos produced before this feature) and contributes no copy
 *  evidence, so old indices keep working unchanged.
 *
 ********************************************************************************************/

#ifndef _COPYNUM_H
#define _COPYNUM_H

#include <stdint.h>

#define CN_TIER_SHIFT 29
#define CN_LOW_MASK   ((1u << CN_TIER_SHIFT) - 1u)        //  0x1FFFFFFF
#define CN_MAXVAL     (1 << CN_TIER_SHIFT)                //  2^29: values must be below this
#define CN_MAXLEN     (1 << 24)                           //  2^24: read length must be below this
                                                          //  (so the radix sort key stays in bytes 4..6)

#define CN_TIER(x)    ((unsigned)(x) >> CN_TIER_SHIFT)            //  extract tier (0..6)
#define CN_VAL(x)     ((int)((unsigned)(x) & CN_LOW_MASK))       //  extract real value
#define CN_PACK(v, t) ((int)(((unsigned)(v) & CN_LOW_MASK) | ((unsigned)(t) << CN_TIER_SHIFT)))

#define CN_NTIER 7    //  tiers 0..6 (0 = unknown)

//  Quantize a genome k-mer copy number to a tier.  Lower tier = rarer = more
//  informative for haplotype/locus discrimination.
//     1 = unique (copy 1)      2 = copy 2        3 = copy 3-4
//     4 = copy 5-8             5 = copy 9-16     6 = copy >16
static inline int cn_tier_of(uint32_t g) {
  if (g <= 1)  return 1;
  if (g == 2)  return 2;
  if (g <= 4)  return 3;
  if (g <= 8)  return 4;
  if (g <= 16) return 5;
  return 6;
}

#endif  //  _COPYNUM_H
