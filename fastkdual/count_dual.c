/*********************************************************************************************\
 *
 *  Phase 2 for FastKdual.
 *
 *  The ordinary FastK implementation has one inverse-profile channel selected by
 *  PRO_TABLE.  FastKdual needs both channels, so it builds one inverse list whose records
 *  share the position index and contain both self and relative counts, then completes the
 *  two profile channels sequentially.  The existing FastK implementation is included under a private
 *  Sorting name so that its well-tested list, table, profile, and sort helpers can be reused
 *  without changing count.c or the FastK executable.
 *
 *********************************************************************************************/

#define Sorting FastK_Single_Sorting_Not_Used
#include "count.c"
#undef Sorting

#include <errno.h>
#include <time.h>

#include "FastKdual.h"

static int DUAL_CMER_WORD;

static void dual_read_exact(int fd, void *buffer, size_t length, char *name)
{
  uint8  *cursor = (uint8 *) buffer;
  size_t  done = 0;

  while (done < length)
    { ssize_t count = read(fd,cursor+done,length-done);
      if (count > 0)
        done += (size_t) count;
      else if (count < 0 && errno == EINTR)
        continue;
      else
        { fprintf(stderr,"%s: Cannot read complete record from %s: %s\n",
                           Prog_Name,name,
                           count == 0 ? "unexpected end of file" : strerror(errno));
          Clean_Exit(1);
        }
    }
}

static void dual_write_exact(int fd, const void *buffer, size_t length, char *name)
{
  const uint8 *cursor = (const uint8 *) buffer;
  size_t       done = 0;

  while (done < length)
    { ssize_t count = write(fd,cursor+done,length-done);
      if (count > 0)
        done += (size_t) count;
      else if (count < 0 && errno == EINTR)
        continue;
      else
        { fprintf(stderr,"%s: Cannot write complete record to %s: %s\n",
                           Prog_Name,name,
                           count == 0 ? "zero-length write" : strerror(errno));
          Clean_Exit(1);
        }
    }
}

static void *dual_realloc(void *buffer, int64 size, char *message)
{
  void *resized = Realloc(buffer,size,message);

  if (resized == NULL)
    Clean_Exit(1);
  return resized;
}

typedef struct
  { struct timespec wall;
    struct timespec cpu;
  } Dual_Time;

static Dual_Time dual_time_now(void)
{
  Dual_Time now;

  clock_gettime(CLOCK_MONOTONIC,&now.wall);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&now.cpu);
  return now;
}

static void dual_time_add(Dual_Time start, double *wall, double *cpu)
{
  Dual_Time now = dual_time_now();

  *wall += (now.wall.tv_sec-start.wall.tv_sec)
         + (now.wall.tv_nsec-start.wall.tv_nsec)/1.e9;
  *cpu  += (now.cpu.tv_sec-start.cpu.tv_sec)
         + (now.cpu.tv_nsec-start.cpu.tv_nsec)/1.e9;
}

static void dual_thread_create(THREAD *thread, void *(*worker)(void *), void *arg)
{
  int error = pthread_create(thread,NULL,worker,arg);

  if (error != 0)
    { fprintf(stderr,"%s: Cannot create worker thread: %s\n",
                     Prog_Name,strerror(error));
      Clean_Exit(1);
    }
}

static void dual_thread_join(THREAD thread)
{
  int error = pthread_join(thread,NULL);

  if (error != 0)
    { fprintf(stderr,"%s: Cannot join worker thread: %s\n",
                     Prog_Name,strerror(error));
      Clean_Exit(1);
    }
}

static void unlink_if_present(char *path)
{
  if (unlink(path) != 0 && errno != ENOENT)
    { fprintf(stderr,"%s: Cannot remove temporary file %s: %s\n",
                     Prog_Name,path,strerror(errno));
      Clean_Exit(1);
    }
}

static void remove_partition_table(char *fname, char *root, int part)
{
  int t;

  sprintf(fname,"%s/%s.U%d.ktab",SORT_PATH,root,part);
  unlink_if_present(fname);
  for (t = 1; t <= NTHREADS; t++)
    { sprintf(fname,"%s/.%s.U%d.ktab.%d",SORT_PATH,root,part,t);
      unlink_if_present(fname);
    }
}

/*
 * table_write_thread writes the bucket prefix back into the first byte of each group leader.
 * That is harmless in ordinary self mode, but prefix zero then erases a group marker needed
 * by the relative merge.  Emit the same table record without modifying the shared k-mer sort.
 */
static void *dual_table_write_thread(void *arg)
{
  Twrite_Arg *data  = (Twrite_Arg *) arg;
  int         beg   = data->beg;
  int         end   = data->end;
  int64      *part  = data->parts;
  int         kfile = data->kfile;

  uint8  bufr[0x10000];
  uint8 *fill, *bend;
  int    x, ct;
  uint8 *kptr, *lptr, *kend;
  int64  tmer;

  fill = bufr;
  bend = bufr+(0x10000-TMER_WORD);
  tmer = 0;
  kptr = data->sort+data->off;

  for (x = beg; x < end; x++)
    for (kend = kptr+part[x]; kptr < kend; kptr = lptr)
      { ct = *((uint16 *) (kptr+KMER_BYTES));
        lptr = kptr+KMER_WORD;
        while (*lptr == 0)
          lptr += KMER_WORD;

        if (ct >= DO_TABLE)
          { if (fill >= bend)
              { dual_write_exact(kfile,bufr,(size_t)(fill-bufr),data->kname);
                fill = bufr;
              }
            *fill++ = x;
            memcpy(fill,kptr+1,TMER_WORD-1);
            fill += TMER_WORD-1;
            tmer += 1;
          }
      }

  if (fill > bufr)
    dual_write_exact(kfile,bufr,(size_t)(fill-bufr),data->kname);

  data->tmers = tmer;
  return NULL;
}

typedef struct
  { uint8       *sort;
    int64       *parts;
    int          beg;
    int          end;
    int64        off;
    uint8       *fours[256];
    Kmer_Stream *stm;
  } Dual_Clist_Arg;

/*
 * Merge the weighted k-mer list with the relative table once and emit one inverse record
 * containing both counts:
 *
 *   self:uint16, relative:uint16, position-index-without-its-LSD-byte
 *
 * The omitted position byte is represented by the 256 pre-partitioned fill ranges, exactly
 * as in cmer_list_thread/cmer_merge_thread.
 */
static void *dual_cmer_merge_thread(void *arg)
{
  Dual_Clist_Arg *data  = (Dual_Clist_Arg *) arg;
  int             beg   = data->beg;
  int             end   = data->end;
  int64          *part  = data->parts;
  uint8         **fours = data->fours;
  Kmer_Stream    *S     = data->stm;
  int             hbyte = S->hbyte;

  int KM1 = KMER_WORD-1;

  int    x, d, k, c;
  int    first;
  uint8 *kptr, *lptr, *kend;
  uint8 *fill;
  int    ct;

  first = beg;
  while (first < end && part[first] == 0)
    first += 1;
  if (first >= end)
    return NULL;

  if (beg > 0)
    S = Clone_Kmer_Stream(S);

  kptr = data->sort+data->off;
  *kptr = first;
  GoTo_Kmer_Entry(S,kptr);

  for (x = beg; x < end; x++)
    { while (S->cpre < x)
        Next_Kmer_Entry(S);

      for (kend = kptr+part[x]; kptr < kend; )
        { lptr = kptr+KMER_WORD;
          while (*lptr == 0)
            lptr += KMER_WORD;

          while (1)
            { if (S->cpre > x)
                { c = 1;
                  break;
                }
              c = memcmp(S->csuf,kptr+1,hbyte);
              if (c >= 0)
                break;
              Next_Kmer_Entry(S);
            }

          if (c != 0)
            ct = 0;
          else
            { ct = *((uint16 *) (S->csuf+hbyte));
              Next_Kmer_Entry(S);
            }

          while (kptr < lptr)
            { d = kptr[KM1];
              fill = fours[d];
              *((uint16 *) fill) = *((uint16 *) (kptr+KMER_BYTES));
              fill += 2;
              *((uint16 *) fill) = ct;
              fill += 2;
              for (k = TMER_WORD; k < KM1; k++)
                *fill++ = kptr[k];
              fours[d] = fill;
              kptr += KMER_WORD;
            }
        }
    }

  if (beg > 0)
    Free_Kmer_Stream(S);
  return NULL;
}

typedef struct
  { uint8 *sort;
    int64 *parts;
    int    beg;
    int    end;
    int64  off;
    uint8 *counts;
    int64  cnts;
    int    count_offset;
    uint8 *prol;
    uint8 *fill;
  } Dual_Plist_Arg;

/*
 * Dual-channel variant of profile_list_thread.  The original helper clears the last-fragment
 * bit in the sorted super-mer array.  Here that bit is masked only while copying the run id,
 * leaving the common sorted input immutable for the second channel.
 */
static void *dual_profile_list_thread(void *arg)
{
  Dual_Plist_Arg *data = (Dual_Plist_Arg *) arg;
  int        beg       = data->beg;
  int        end       = data->end;
  int64     *part      = data->parts;

  uint8 *cnts, *fill, *prof;
  uint8 *sptr, *send;
  int    i, x;
  uint8 *asp;
  uint8  b[2*MAX_SUPER+2];
  uint64 pidx;

  int STOT = SMER_BYTES + SLEN_BYTES;

  int sln = 0;
  int len = 0;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  uint8 *sb = ((uint8 *) &sln) - 1;
  uint8 *lb = ((uint8 *) &len) - 1;
#else
  uint8 *sb = ((uint8 *) &sln) + (sizeof(int)-SLEN_BYTES);
  uint8 *lb = ((uint8 *) &len) + (sizeof(int)-PLEN_BYTES);
#endif

  prof = data->prol;
  pidx = CMER_WORD*data->cnts;
  cnts = data->counts+DUAL_CMER_WORD*data->cnts+data->count_offset;
  sptr = data->sort+data->off;
  fill = data->fill;

  for (x = beg; x < end; x++)
    for (send = sptr+part[x]; sptr < send; )
      { asp = sptr+SMER_BYTES;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        for (i = SLEN_BYTES; i > 0; i--)
          sb[i] = *asp++;
#else
        for (i = 0; i < SLEN_BYTES; i++)
          sb[i] = *asp++;
#endif

        { int j, p, c, d;
          int run;
          uint8 *db = (uint8 *) &d;

          p = *((uint16 *) cnts);
          cnts += DUAL_CMER_WORD;
          *((uint16 *) b) = p;
          len = 2;
          run = 0;

          for (j = 1; j <= sln; j++)
            { c = *((uint16 *) cnts);
              cnts += DUAL_CMER_WORD;

              if (c == p)
                { if (run > 0)
                    { if (run >= 63)
                        { b[len++] = run;
                          run = 1;
                        }
                      else
                        run += 1;
                    }
                  else
                    run = 1;
                }
              else
                { if (run > 0)
                    { b[len++] = run;
                      run = 0;
                    }
                  d = c-p;
                  if (abs(d) < 32)
                    b[len++] = 0x40 | (d & 0x3f);
                  else
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
                    { b[len++] = db[1] | 0x80;
                      b[len++] = db[0];
#else
                    { b[len++] = db[0] | 0x80;
                      b[len++] = db[1];
#endif
                    }
                }
              p = c;
            }

          if (run > 0)
            b[len++] = run;
          if (len > 2)
            *((uint16 *) (b+len)) = p;
        }

        do
          { if ((sptr[STOT] & 0x80) != 0)
              *((uint64 *) fill) = (pidx<<1) | 0x1;
            else
              *((uint64 *) fill) = (pidx<<1);
            fill += sizeof(uint64);

            *fill++ = sptr[STOT] & 0x7f;
            for (i = STOT+1; i < SMER_WORD; i++)
              *fill++ = sptr[i];
            sptr += SMER_WORD;
          }
        while (*sptr == 0);

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        for (i = PLEN_BYTES; i > 0; i--)
          prof[pidx++] = lb[i];
#else
        for (i = 0; i < PLEN_BYTES; i++)
          prof[pidx++] = lb[i];
#endif
        if (len > 2)
          len += 2;
        for (i = 0; i < len; i++)
          prof[pidx++] = b[i];
      }

  return NULL;
}

/*
 * Complete one inverse-list channel.
 */
static void finish_profile_channel(uint8 *dual_counts, uint8 *profile_buffer,
                                   int count_offset, uint8 *s_sort, int64 nmers,
                                   int64 *Sparts,
                                   Klist_Arg *parmk, Slist_Arg *parms,
                                   Dual_Plist_Arg *parmp, Pwrite_Arg *parmw,
                                   THREAD *threads, char *fname, char *channel_root,
                                   int part)
{
  uint8 *p_sort, *a_sort;
  int    t;

  p_sort = Malloc(nmers*PROF_BYTES*2,"Allocating dual profile link array");
  if (p_sort == NULL)
    Clean_Exit(1);

  for (t = 0; t < NTHREADS; t++)
    { parmp[t].sort  = s_sort;
      parmp[t].parts = Sparts;
      parmp[t].beg   = parmk[t].beg;
      parmp[t].end   = parmk[t].end;
      parmp[t].off   = parmk[t].off;
      parmp[t].counts = dual_counts;
      parmp[t].cnts  = parmk[t].kidx;
      parmp[t].count_offset = count_offset;
      parmp[t].prol  = profile_buffer;
      parmp[t].fill  = p_sort + (parmk[t].off / SMER_WORD) * PROF_BYTES;
    }

  for (t = 1; t < NTHREADS; t++)
    dual_thread_create(threads+t,dual_profile_list_thread,parmp+t);
  dual_profile_list_thread(parmp);
  for (t = 1; t < NTHREADS; t++)
    dual_thread_join(threads[t]);

  { int i, x;
    int bytes[PROF_BYTES+1];

    x = 0;
    for (i = PROF_BYTES-1; i >= (int) sizeof(uint64); i--)
      bytes[x++] = i;
    bytes[x] = -1;

    a_sort = LSD_Sort(nmers,p_sort,p_sort+nmers*PROF_BYTES,PROF_BYTES,bytes);
  }

  { int64 o;

    sprintf(fname,"%s/%s.%d.P",SORT_PATH,channel_root,part);
    o = 0;
    for (t = 0; t < ITHREADS; t++)
      { parmw[t].sort  = a_sort;
        parmw[t].beg   = o;
        o += parms[t].nmers;
        parmw[t].end   = o;
        parmw[t].nidxs = parms[t].nidxs;
        parmw[t].prol  = profile_buffer;
        parmw[t].root  = fname;
        parmw[t].wch   = t;
      }
  }

  for (t = 1; t < ITHREADS; t++)
    dual_thread_create(threads+t,profile_write_thread,parmw+t);
  profile_write_thread(parmw);
  for (t = 1; t < ITHREADS; t++)
    dual_thread_join(threads[t]);

  free(p_sort);
}

void Sorting_Dual(char *path, char *root, char *relative_root)
{
  char  *fname;
  int64  counts[0x8000];
  int64  max_inst;
  int   *reload;
  double common_wall = 0., common_cpu = 0.;
  double self_wall = 0., self_cpu = 0.;
  double relative_wall = 0., relative_cpu = 0.;

  if (VERBOSE)
    { fprintf(stderr,"\nPhase 2: Sorting once and producing self + relative profiles in %d blocks\n\n",
                     NPARTS);
      fflush(stderr);
    }

  fname = Malloc(2*(strlen(SORT_PATH) + strlen(path) + strlen(root)
                    + strlen(relative_root)) + 200,"File name buffer");
  if (fname == NULL)
    Clean_Exit(1);

  { int  i, n;
    int  l0, l1, l2, l3;
    int *s;

    reload = Malloc(sizeof(int)*IO_UBITS*MAX_SUPER,"Allocating reload table");
    if (reload == NULL)
      Clean_Exit(1);

    for (i = 0; i < 0x8000; i++)
      counts[i] = 0;
    max_inst = 0;

    s = reload;
    for (i = 0; i < IO_UBITS; i++)
      { Fixed_Reload[IO_UBITS-i] = (i + SLEN_BITS - 1) / IO_UBITS;
        if (RUN_BITS < 17)
          Runer_Reload[IO_UBITS-i] = (i + 18) / IO_UBITS;
        else
          Runer_Reload[IO_UBITS-i] = (i + RUN_BITS - 1) / IO_UBITS;
        Super_Reload[IO_UBITS-i] = s;
        for (n = 0; n < MAX_SUPER; n++)
          *s++ = (i + 2*(n+KMER) - 1) / IO_UBITS;
      }

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
  }

  { Slist_Arg  *parms = Malloc(sizeof(Slist_Arg)*ITHREADS,"Allocating sort controls");
    Klist_Arg  *parmk = Malloc(sizeof(Klist_Arg)*NTHREADS,"Allocating sort controls");
    Dual_Clist_Arg *parmc =
      Malloc(sizeof(Dual_Clist_Arg)*NTHREADS,"Allocating dual inverse controls");
    Twrite_Arg *parmt = Malloc(sizeof(Twrite_Arg)*NTHREADS,"Allocating sort controls");
    Dual_Plist_Arg *parmp =
      Malloc(sizeof(Dual_Plist_Arg)*NTHREADS,"Allocating dual profile controls");
    Pwrite_Arg *parmw = Malloc(sizeof(Pwrite_Arg)*ITHREADS,"Allocating sort controls");

    int   *Table_Split = Malloc(sizeof(int)*NTHREADS,"Allocating sort controls");
    int64 *Sparts      = Malloc(sizeof(int64)*256,"Allocating sort controls");
    int64 *Kparts      = Malloc(sizeof(int64)*256,"Allocating sort controls");
    Range *Panels      = Malloc(sizeof(Range)*NTHREADS,"Allocating sort controls");
    int64 *Wkmers      = Malloc(sizeof(int64)*NPARTS,"Allocating sort controls");
    int64 *Ukmers      = Malloc(sizeof(int64)*NPARTS,"Allocating sort controls");
    THREAD *threads    = Malloc(sizeof(THREAD)*NTHREADS,"Allocating sort controls");

    uint8 *s_sort;
    int64  tmers;
    int    t, p;

    if (parms == NULL || parmk == NULL || parmc == NULL || parmt == NULL ||
        parmp == NULL || parmw == NULL || Table_Split == NULL || Sparts == NULL ||
        Kparts == NULL || Panels == NULL || Wkmers == NULL || Ukmers == NULL ||
        threads == NULL)
      Clean_Exit(1);

    KMER_WORD += KMAX_BYTES;
    CMER_WORD  = KMAX_BYTES+1;
    DUAL_CMER_WORD = KMAX_BYTES+3;
    RUN_BYTES  = (RUN_BITS+7) >> 3;
    SMER_WORD += RUN_BYTES;
    PROF_BYTES = RUN_BYTES + sizeof(uint64);

    s_sort = Malloc((NMAX+1)*SMER_WORD+1,"Allocating super-mer sort array");
    if (s_sort == NULL)
      Clean_Exit(1);
    *s_sort++ = 0;

    tmers = 0;
    for (p = 0; p < NPARTS; p++)
      { uint8 *k_sort, *dual_inverse, *dual_counts, *profile_buffer;
        int64  kmers, nmers, skmers;
        Dual_Time phase_time = dual_time_now();

        kmers = 0;
        nmers = 0;
        for (t = 0; t < ITHREADS; t++)
          { int64 k, n;
            int   f;

            sprintf(fname,"%s/%s.%d.T%d",SORT_PATH,root,p,t);
            f = open(fname,O_RDONLY);
            if (f < 0)
              { fprintf(stderr,"\n%s: File %s should exist but doesn't?\n",Prog_Name,fname);
                Clean_Exit(1);
              }

            parms[t].tfile = f;
            parms[t].nidxs = NUM_RID[t];
            dual_read_exact(f,&k,sizeof(int64),fname);
            dual_read_exact(f,&n,sizeof(int64),fname);
            parms[t].nmers = n;
            kmers += k;
            nmers += n;
            dual_read_exact(f,Panels[t].khist,sizeof(int64)*256,fname);
          }

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Sorting super-mers     ",p+1);
            fflush(stderr);
          }

        { int64 o, x;
          int   j;

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < ITHREADS; t++)
              { parms[t].fours[j] = s_sort + o*SMER_WORD;
                o += Panels[t].khist[j];
              }

          o = 0;
          for (t = 0; t < ITHREADS; t++)
            { x = parms[t].nidxs;
              parms[t].nbase = o;
              o += x;
            }
        }

        for (t = 1; t < ITHREADS; t++)
          dual_thread_create(threads+t,supermer_list_thread,parms+t);
        supermer_list_thread(parms);
        for (t = 1; t < ITHREADS; t++)
          dual_thread_join(threads[t]);

        for (t = 0; t < ITHREADS; t++)
          { close(parms[t].tfile);
            sprintf(fname,"%s/%s.%d.T%d",SORT_PATH,root,p,t);
            unlink(fname);
          }

        { uint8 *o, *x;
          int    j;

          o = s_sort;
          for (j = 0; j < 256; j++)
            { x = parms[ITHREADS-1].fours[j];
              Sparts[j] = x-o;
              o = x;
            }
        }

        Supermer_Sort(s_sort,nmers,SMER_WORD,SMER_BYTES+SLEN_BYTES,
                      Sparts,NTHREADS,Panels);

        { int64 o, x;
          int   j;

          o = 0;
          for (t = 0; t < NTHREADS; t++)
            { parmk[t].kidx = o;
              for (j = 0; j < 256; j++)
                o += Panels[t].khist[j];
            }

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < NTHREADS; t++)
              { x = Panels[t].khist[j];
                Panels[t].khist[j] = o;
                o += x;
              }
          skmers = o;
        }

        Wkmers[p] = skmers;
        Ukmers[p] = kmers;
        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Sorting weighted k-mers",p+1);
            fflush(stderr);
          }

        k_sort = Malloc(skmers*KMER_WORD+1,"Allocating weighted k-mer list");
        dual_inverse = Malloc(skmers*DUAL_CMER_WORD+1,"Allocating dual-count inverse list");
        if (k_sort == NULL || dual_inverse == NULL)
          Clean_Exit(1);

        for (t = 0; t < NTHREADS; t++)
          { int j;

            parmk[t].sort  = s_sort;
            parmk[t].parts = Sparts;
            parmk[t].beg   = Panels[t].beg;
            parmk[t].end   = Panels[t].end;
            parmk[t].off   = Panels[t].off;
            for (j = 0; j < 256; j++)
              parmk[t].fours[j] = k_sort + Panels[t].khist[j]*KMER_WORD;
          }

        for (t = 1; t < NTHREADS; t++)
          dual_thread_create(threads+t,kmer_list_thread,parmk+t);
        kmer_list_thread(parmk);
        for (t = 1; t < NTHREADS; t++)
          dual_thread_join(threads[t]);

        { uint8 *o, *x;
          int    j;

          o = k_sort;
          for (j = 0; j < 256; j++)
            { x = parmk[NTHREADS-1].fours[j];
              Kparts[j] = x-o;
              o = x;
            }
        }

        Weighted_Kmer_Sort(k_sort,skmers,KMER_WORD,KMER_BYTES,
                           Kparts,NTHREADS,Panels);
        dual_time_add(phase_time,&common_wall,&common_cpu);

        { int64 *ncnt;
          int    x;

          for (t = 0; t < NTHREADS; t++)
            { ncnt = Panels[t].count;
              for (x = 1; x < 0x8000; x++)
                counts[x] += ncnt[x];
              max_inst += Panels[t].max_inst + parmk[t].overflow;
            }
        }

        if (DO_TABLE > 0)
          { phase_time = dual_time_now();

            if (p == 0)
              for (t = 0; t < NTHREADS; t++)
                { Table_Split[t] = Panels[t].beg;
                  parmt[t].off = Panels[t].off;
                }
            else
              { int64 off;
                int   beg;

                off = 0;
                beg = 0;
                for (t = 0; t < NTHREADS; t++)
                  { while (beg < Table_Split[t])
                      { off += Kparts[beg];
                        beg += 1;
                      }
                    parmt[t].off = off;
                  }
              }

            for (t = 0; t < NTHREADS; t++)
              { parmt[t].sort  = k_sort;
                parmt[t].parts = Kparts;
                parmt[t].beg   = Table_Split[t];
                parmt[t].end   = (t < NTHREADS-1 ? Table_Split[t+1] : 256);
                sprintf(fname,"%s/%s.%d.L%d",SORT_PATH,root,p,t);
                parmt[t].kfile = open(fname,O_WRONLY|O_CREAT|O_TRUNC,
                                      S_IRWXU|S_IRWXG|S_IRWXO);
                parmt[t].kname = Strdup(fname,"Allocating stream name");
                if (parmt[t].kfile < 0 || parmt[t].kname == NULL)
                  Clean_Exit(1);
              }

            for (t = 1; t < NTHREADS; t++)
              dual_thread_create(threads+t,dual_table_write_thread,parmt+t);
            dual_table_write_thread(parmt);
            for (t = 1; t < NTHREADS; t++)
              dual_thread_join(threads[t]);

            for (t = 0; t < NTHREADS; t++)
              tmers += parmt[t].tmers;

            if (p == NPARTS-1)
              { if (tmers > 0x4000000ll && KMER >= 12)
                  IDX_BYTES = 3;
                else if (tmers >= 0x40000ll && KMER >= 8)
                  IDX_BYTES = 2;
                else
                  IDX_BYTES = 1;
              }

            for (t = 0; t < NTHREADS; t++)
              { free(parmt[t].kname);
                close(parmt[t].kfile);
              }

            dual_time_add(phase_time,&self_wall,&self_cpu);
          }

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Building dual-count inverse list",p+1);
            fflush(stderr);
          }

        { int64 o, x;
          int   j;

          o = 0;
          for (j = 0; j < 256; j++)
            for (t = 0; t < NTHREADS; t++)
              { x = Panels[t].khist[j];
                Panels[t].khist[j] = o;
                o += x;
              }
        }

        for (t = 0; t < NTHREADS; t++)
          { int j;

            parmc[t].sort  = k_sort;
            parmc[t].parts = Kparts;
            parmc[t].beg   = Panels[t].beg;
            parmc[t].end   = Panels[t].end;
            parmc[t].off   = Panels[t].off;
            for (j = 0; j < 256; j++)
              parmc[t].fours[j] = dual_inverse + Panels[t].khist[j]*DUAL_CMER_WORD;
          }

        sprintf(fname,"%s/%s.U%d",SORT_PATH,root,p);
        parmc[0].stm = Open_Kmer_Stream(fname);
        if (parmc[0].stm == NULL)
          { fprintf(stderr,"\n%s: Table %s should exist but doesn't?\n",Prog_Name,fname);
            Clean_Exit(1);
          }
        for (t = 1; t < NTHREADS; t++)
          parmc[t].stm = parmc[0].stm;

        phase_time = dual_time_now();
        for (t = 1; t < NTHREADS; t++)
          dual_thread_create(threads+t,dual_cmer_merge_thread,parmc+t);
        dual_cmer_merge_thread(parmc);
        for (t = 1; t < NTHREADS; t++)
          dual_thread_join(threads[t]);

        Free_Kmer_Stream(parmc[0].stm);
        remove_partition_table(fname,root,p);
        dual_time_add(phase_time,&relative_wall,&relative_cpu);

        phase_time = dual_time_now();
        { int i, x;
          int bytes[KMAX_BYTES];

          x = 0;
          for (i = DUAL_CMER_WORD-1; i >= 4; i--)
            bytes[x++] = i;
          bytes[x] = -1;

          dual_counts = LSD_Sort(skmers,dual_inverse,k_sort,DUAL_CMER_WORD,bytes);
        }

        /*
         * Keep only the sorted dual record and one CMER-sized profile workspace.  Reallocating
         * both backing arrays here is what removes the large resident k-mer scratch before
         * either PROF_BYTES link array is allocated.
         */
        if (dual_counts == dual_inverse)
          { dual_inverse = dual_realloc(dual_inverse,skmers*DUAL_CMER_WORD+1,
                                        "Pruning sorted dual-count list");
            dual_counts = dual_inverse;
            k_sort = dual_realloc(k_sort,skmers*CMER_WORD+1,
                                  "Pruning profile workspace");
            profile_buffer = k_sort;
          }
        else
          { k_sort = dual_realloc(k_sort,skmers*DUAL_CMER_WORD+1,
                                  "Pruning sorted dual-count list");
            dual_counts = k_sort;
            dual_inverse = dual_realloc(dual_inverse,skmers*CMER_WORD+1,
                                        "Pruning profile workspace");
            profile_buffer = dual_inverse;
          }
        dual_time_add(phase_time,&relative_wall,&relative_cpu);

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Writing self profile fragments",p+1);
            fflush(stderr);
          }
        phase_time = dual_time_now();
        finish_profile_channel(dual_counts,profile_buffer,0,s_sort,nmers,Sparts,
                               parmk,parms,parmp,parmw,threads,fname,root,p);
        dual_time_add(phase_time,&self_wall,&self_cpu);

        if (VERBOSE)
          { fprintf(stderr,"\r  Processing block %d: Writing relative profile fragments",p+1);
            fflush(stderr);
          }
        phase_time = dual_time_now();
        finish_profile_channel(dual_counts,profile_buffer,2,s_sort,nmers,Sparts,
                               parmk,parms,parmp,parmw,threads,fname,relative_root,p);
        dual_time_add(phase_time,&relative_wall,&relative_cpu);

        free(dual_inverse);
        free(k_sort);
      }

    free(s_sort-1);

    if (VERBOSE)
      { int64 wtot, utot;
        double psav;
        int p;

        wtot = utot = 0;
        psav = 0.;
        for (p = 0; p < NPARTS; p++)
          { wtot += Wkmers[p];
            utot += Ukmers[p];
            if (Wkmers[p] > 0 && (1.*Ukmers[p])/Wkmers[p] > psav)
              psav = (1.*Ukmers[p])/Wkmers[p];
          }
        fprintf(stderr,"\r                                               \r");
        fprintf(stderr,"  Dual weighted k-mers: ");
        Print_Number(wtot,0,stderr);
        fprintf(stderr," (unweighted ");
        Print_Number(utot,0,stderr);
        fprintf(stderr,", max savings %.1f)\n",psav);
        fprintf(stderr,"  Dual phase-2 cumulative resources:\n");
        fprintf(stderr,"    common super-mer + weighted-k-mer sort: %.3fs CPU, %.3fs wall\n",
                       common_cpu,common_wall);
        fprintf(stderr,"    self table + profile fragments:         %.3fs CPU, %.3fs wall\n",
                       self_cpu,self_wall);
        fprintf(stderr,"    relative join + dual inverse + fragments: %.3fs CPU, %.3fs wall\n",
                       relative_cpu,relative_wall);
      }

    free(threads);
    free(Ukmers);
    free(Wkmers);
    free(Panels);
    free(Kparts);
    free(Sparts);
    free(Table_Split);
    free(parmw);
    free(parmp);
    free(parmt);
    free(parmc);
    free(parmk);
    free(parms);
  }

  { int i, f;

    sprintf(fname,"%s/%s.hist",path,root);
    f = open(fname,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
    if (f < 0)
      Clean_Exit(1);
    FastKdual_Mark_Histogram_Touched();
    dual_write_exact(f,&KMER,sizeof(int),fname);
    i = 1;
    dual_write_exact(f,&i,sizeof(int),fname);
    i = 0x7fff;
    dual_write_exact(f,&i,sizeof(int),fname);
    dual_write_exact(f,counts+1,sizeof(int64),fname);
    dual_write_exact(f,&max_inst,sizeof(int64),fname);
    dual_write_exact(f,counts+1,0x7fff*sizeof(int64),fname);
    close(f);
  }

  free(reload);
  free(fname);
}
