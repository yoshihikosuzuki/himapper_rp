#ifndef HIMAPPER_MINIMAXMER_H
#define HIMAPPER_MINIMAXMER_H

#include <stddef.h>
#include <stdint.h>
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
    uint64_t         resident_bytes;
    int              block;
  } Minimaxmer_Block;

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
                                    int kmer, int window, int spacing);
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

/* Optional KPOS tee.  Benchmarks can leave this unopened to avoid all KPOS I/O. */
Minimaxmer_Kpos_Writer *Minimaxmer_Open_Kpos_Writer(const char *path,
                                                     int kmer);
int Minimaxmer_Write_Kpos_Block(Minimaxmer_Kpos_Writer *writer,
                                const Minimaxmer_Block *block);
int Minimaxmer_Close_Kpos_Writer(Minimaxmer_Kpos_Writer *writer);

#endif
