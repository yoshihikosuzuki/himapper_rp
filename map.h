/*******************************************************************************************
 *
 *  Filter interface for the dazzler.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#ifndef _FILTER

#define _FILTER

#include <DB.h>
#include "align.h"

extern int    VERBOSE;
extern int    PROFILE;
extern int    SPACING;
extern double BEST_TIE;
extern int    PRIMARY_ONLY;
extern int    JUNCTION_PEN;   //  collinear-gap penalties (default on; -G disables)
extern char  *SORT_PATH;

extern uint64 MEM_LIMIT;
extern uint64 MEM_PHYSICAL;

int Set_Filter_Params(int kmer, int suppress, int nthreads); 

void *Sort_Kmers(DAZZ_DB *block, int *len);

typedef struct Kpos_Index Kpos_Index;

// Open a block-indexed KPOS container.  The sidecar path is <path>.idx and
// contains nblocks+1 cumulative uint64 entry offsets.
Kpos_Index *Open_Kpos_Index(const char *path, int nblocks);
void Close_Kpos_Index(Kpos_Index *index);

// Load one 1-based DB block.  Entries on disk are code-sorted within each
// block.  comp=0 preserves that order; comp=1 reverse-complements and re-sorts.
// The returned array has the same +2-sentinel contract as Sort_Kmers.
void *Load_Kpos_Block(Kpos_Index *index, int block_no, int comp,
                      DAZZ_DB *block, int *len);

// Memory retained by the caller in addition to the current read/reference
// block pair (for example the whole reference DB and cached KPOS arrays).
void Set_Filter_Extra_Memory(uint64 bytes);

// Reset the per-read candidate state before comparing one reads DB block
// against the reference blocks.
void Begin_Mapping_Block(DAZZ_DB *block);

void Match_Filter(DAZZ_DB *ablock, DAZZ_DB *bblock,
                  void *atable, int alen, void *btable, int blen, int comp, int start);

#define FLAG_DOA    0x1
#define FLAG_DOB    0x2

void Reporter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
              Align_Spec *asettings, int mflag);

void Clean_Exit(int val);

#endif
