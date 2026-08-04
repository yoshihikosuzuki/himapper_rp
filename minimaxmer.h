#ifndef HIMAPPER_MINIMAXMER_H
#define HIMAPPER_MINIMAXMER_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <DB.h>

typedef struct
  { uint64_t code;
    int32_t  rpos;
    uint32_t rid;
  } Minimaxmer_Kpos;

typedef struct
  { Minimaxmer_Kpos *entry;
    DAZZ_DB          *loaded_db;
    size_t           len;
    size_t           capacity;
    uint64_t         resident_bytes;
    int              block;
    int              streaming;
  } Minimaxmer_Block;

/*
 * Compute-phase times are sums of per-worker monotonic elapsed time.  They
 * describe where seed-generation work is spent, not critical-path wall time.
 * kpos_write_ns is instead control-thread serial elapsed time.
 */
typedef struct
  { uint64_t block_load_ns;
    uint64_t sequence_read_ns;
    uint64_t profile_fetch_ns;
    uint64_t hash_ns;
    uint64_t window_ns;
    uint64_t candidate_build_ns;
    uint64_t candidate_sort_ns;
    uint64_t thinning_ns;
    uint64_t code_emit_ns;
    uint64_t block_sort_ns;
    uint64_t kpos_write_ns;
    uint64_t generated_blocks;
    uint64_t generated_reads;
    uint64_t generated_seeds;
    /* Appended to preserve the offsets of the existing phase-stat prefix. */
    uint64_t profile_open_ns;
    uint64_t profile_seek_ns;
    uint64_t profile_read_ns;
    uint64_t profile_decode_control_ns;
    uint64_t profile_unclassified_ns;
    uint64_t profile_fetch_calls;
    uint64_t profile_open_calls;
    uint64_t profile_seek_calls;
    uint64_t profile_refill_calls;
    uint64_t profile_read_calls;
    uint64_t profile_read_bytes;
    /* Periodic-seed exclusion counters are collected whenever enabled. */
    uint64_t periodic_filter_ns;
    uint64_t periodic_tracts[3];
    uint64_t periodic_tract_bases[3];
    uint64_t periodic_reads;
    uint64_t periodic_seeds_removed;
  } Minimaxmer_Phase_Stats;

typedef struct Minimaxmer_Context Minimaxmer_Context;
typedef struct Minimaxmer_Kpos_Writer Minimaxmer_Kpos_Writer;

/*
 * Open the whole, partitioned DB and the two FastK profiles used by seed
 * selection.  Selector configuration is process-global, so at most one
 * Minimaxmer_Context may be open at a time.
 */
Minimaxmer_Context *Minimaxmer_Open(const char *db_path,
                                    const char *genome_profile,
                                    const char *reads_profile,
                                    int kmer, int window, int spacing,
                                    int phase_timing);
void Minimaxmer_Close(Minimaxmer_Context *context);
int Minimaxmer_Block_Count(const Minimaxmer_Context *context);

/*
 * Generate a bounded batch without overlapping mapping.  first_block is
 * 1-based and count must be positive.  The returned blocks are caller-owned.
 */
int Minimaxmer_Generate_Batch(Minimaxmer_Context *context,
                              int first_block, int count, int nthreads,
                              Minimaxmer_Block *blocks);

void Minimaxmer_Free_Block(Minimaxmer_Block *block);

/* Retain each decoded DB block with its KPOS block for direct mapping. */
int Minimaxmer_Set_Retain_Loaded_DB(Minimaxmer_Context *context,
                                    int enabled);
DAZZ_DB *Minimaxmer_Block_DB(Minimaxmer_Block *block);

/*
 * Streaming mode keeps at most nthreads block vectors in flight.  Next must
 * be consumed in DB-block order and Release must follow before reusing it.
 */
int Minimaxmer_Start_Stream(Minimaxmer_Context *context, int nthreads);
int Minimaxmer_Next_Block(Minimaxmer_Context *context,
                          Minimaxmer_Block *block);
void Minimaxmer_Release_Block(Minimaxmer_Context *context,
                              Minimaxmer_Block *block);
int Minimaxmer_Stop_Stream(Minimaxmer_Context *context);

/* Decoded DB plus KPOS capacity retained by streaming seed workers. */
uint64_t Minimaxmer_Resident_Bytes(Minimaxmer_Context *context);

/* Open_DB uses internal static path buffers; serialize calls across teams. */
int Minimaxmer_Open_DB_Serialized(char *path, DAZZ_DB *db);

const Minimaxmer_Phase_Stats *Minimaxmer_Get_Phase_Stats(
    const Minimaxmer_Context *context);
void Minimaxmer_Print_Phase_Stats(const Minimaxmer_Context *context,
                                  FILE *out, const char *prefix);

/* Optional KPOS tee.  Benchmarks can leave this unopened to avoid all KPOS I/O. */
Minimaxmer_Kpos_Writer *Minimaxmer_Open_Kpos_Writer(const char *path,
                                                     int kmer);
int Minimaxmer_Write_Kpos_Block(Minimaxmer_Kpos_Writer *writer,
                                const Minimaxmer_Block *block,
                                Minimaxmer_Context *context);
int Minimaxmer_Close_Kpos_Writer(Minimaxmer_Kpos_Writer *writer,
                                 Minimaxmer_Context *context);

#endif
