#ifndef HIMAPPER_PAF_H
#define HIMAPPER_PAF_H

#include "align.h"

extern char *PAF_PATH, *PAF_REF_FAI, *PAF_READ_FAI;
extern int   PAF_HISTORY, PAF_KEEP_LAS, PAF_STATS;

typedef struct
  { char    *name;
    int64_t  length;
    int      offset;
  } Paf_Name;

typedef struct Paf_Source Paf_Source;

Paf_Source *Paf_Open_Source(DAZZ_DB *db, const char *stub, int dam, const char *fai);
Paf_Name *Paf_Read_Names(Paf_Source *src, DAZZ_DB *db);
void Paf_Free_Names(Paf_Name *names, int n);
void Paf_Close_Source(Paf_Source *src);

void Paf_Set_Names(Paf_Name *reads, Paf_Name *reference);
void Paf_Open_Output(void);
void Paf_Close_Output(void);

typedef struct
  { uint64_t key;
    int      ab, bb;
    Cigar    cigar;
  } Cigar_Entry;

typedef struct
  { Cigar_Entry *entry;
    size_t       n, cap, bytes, peak;
  } Cigar_Cache;

void Cigar_Cache_Clear(Cigar_Cache *cache);
void Cigar_Cache_Free(Cigar_Cache *cache);
void Cigar_Cache_Put(Cigar_Cache *cache, uint64_t key, const Path *p, const Cigar *c);
const Cigar *Cigar_Cache_Get(Cigar_Cache *cache, uint64_t key);
void Cigar_Cache_Fuse(Cigar_Cache *cache, uint64_t key, const Path *a,
                       const Path *b, int q, int r);

typedef struct
  { Work_Data *trace_work;
    Cigar      cigar;
    char      *text;
    size_t     text_cap, peak;
    uint64_t   local_calls, selected, fallback, candidate_ops, output_ops;
    uint64_t   exact_tiles, traced_tiles, exact_bases, traced_bases;
    double     local_cpu, finish_cpu, report_cpu;
  } Paf_Work;

double Paf_Thread_Seconds(void);
double Paf_Wall_Seconds(void);

void Paf_Emit(Paf_Work *w, const Overlap *ovl, int read_index,
              const char *query, const char *reverse_query, int qlen,
              const char *reference, int rlen,
              int spacing, const Cigar *saved, int best, int chain, int segment, int mapq);
void Paf_Free_Work(Paf_Work *w);

#endif
