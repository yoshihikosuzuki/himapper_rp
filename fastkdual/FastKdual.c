/********************************************************************************************
 *
 *  FastKdual: build a self table/profile and a relative profile from one input partition
 *  and one super-mer/weighted-k-mer sort.
 *
 *  The original FastK executable is intentionally untouched.  FastKdual has its own CLI and
 *  phase-2 implementation; it reuses the stable input, table distribution, table merge, and
 *  profile merge modules.
 *
 *********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>

#include "libfastk.h"
#include "FastKdual.h"

// The reused FastK profile writer reserves one complete super-mer record in a
// fixed 64-KiB buffer.  This bound keeps that record-size arithmetic in range.
#define FASTKDUAL_MAX_KMER 32761

static char *Usage[] =
  { "[-k<int(40)>] [-t<int(1)>] -p:<table>[.ktab] -R<relative_output_root>",
    "[-c] [-bc<int>] [-v] [-N<self_output_root>] [-P<dir($TMPDIR)>]",
    "[-M<int(12)>] [-T<int(4)>] <source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz] ..."
  };

int    VERBOSE;
int    NTHREADS;
int    ITHREADS;
int64  SORT_MEMORY;
char  *SORT_PATH;

int    KMER;
int    DO_TABLE;
int    DO_PROFILE;
Kmer_Stream *PRO_TABLE;
char        *PRO_NAME;
int    BC_PREFIX;
int    COMPRESS;

int    NPARTS;
int    SMER;
int64  KMAX;
int64  NMAX;
int64  RMAX;

int    MOD_LEN;
int    MOD_MSK;
int    MAX_SUPER;

int    SLEN_BITS;
uint64 SLEN_BIT_MASK;
int    RUN_BITS;
int    RUN_BYTES;
int    SLEN_BYTES;
int    PLEN_BYTES;
int    PROF_BYTES;
int    IDX_BYTES;
int    SLEN_BYTE_MASK;
int    KMER_BYTES;
int    SMER_BYTES;
int    KMAX_BYTES;
int    KMER_WORD;
int    SMER_WORD;
int    TMER_WORD;
int    CMER_WORD;

char *DUAL_ROOT;

static char            *SELF_PATH;
static char            *SELF_ROOT;
static char            *REL_PATH;
static Input_Partition *IOPACK;

static struct rusage   Itime;
static struct timespec Iwall;
static struct rusage   Mtime;
static struct timespec Mwall;
static pthread_mutex_t Cleanup_Lock = PTHREAD_MUTEX_INITIALIZER;
static int             SELF_HISTOGRAM_TOUCHED;
static int             SELF_TABLE_TOUCHED;
static int             SELF_PROFILE_TOUCHED;
static int             RELATIVE_PROFILE_TOUCHED;

void FastKdual_Mark_Histogram_Touched(void)
{
  SELF_HISTOGRAM_TOUCHED = 1;
}

static void start_time(void)
{
  getrusage(RUSAGE_SELF,&Itime);
  clock_gettime(CLOCK_MONOTONIC,&Iwall);
  Mtime = Itime;
  Mwall = Iwall;
}

static void time_to(FILE *f, int all)
{
  struct rusage now, *beg;
  struct timespec wall, *wbeg;
  double cpu, elapsed;

  getrusage(RUSAGE_SELF,&now);
  clock_gettime(CLOCK_MONOTONIC,&wall);
  beg  = all ? &Itime : &Mtime;
  wbeg = all ? &Iwall : &Mwall;

  cpu = (now.ru_utime.tv_sec-beg->ru_utime.tv_sec)
      + (now.ru_utime.tv_usec-beg->ru_utime.tv_usec)/1.e6
      + (now.ru_stime.tv_sec-beg->ru_stime.tv_sec)
      + (now.ru_stime.tv_usec-beg->ru_stime.tv_usec)/1.e6;
  elapsed = (wall.tv_sec-wbeg->tv_sec) + (wall.tv_nsec-wbeg->tv_nsec)/1.e9;

  fprintf(f,"\n%s resources: %.3fs CPU, %.3fs wall, %.1f%%",
            all ? "Total" : "  Phase",cpu,elapsed,elapsed > 0. ? 100.*cpu/elapsed : 0.);
  if (all)
    fprintf(f,", %.3f GB peak RSS",now.ru_maxrss/1.e6);
  fprintf(f,"\n");

  Mtime = now;
  Mwall = wall;
}

static void cleanup_unlink(char *directory, char *format, ...)
{
  char    leaf[PATH_MAX];
  char    path[PATH_MAX];
  va_list args;
  int     length;

  if (directory == NULL || format == NULL)
    return;

  va_start(args,format);
  length = vsnprintf(leaf,sizeof(leaf),format,args);
  va_end(args);
  if (length < 0 || (size_t) length >= sizeof(leaf))
    return;

  length = snprintf(path,sizeof(path),"%s/%s",directory,leaf);
  if (length < 0 || (size_t) length >= sizeof(path))
    return;

  unlink(path);
}

static void cleanup_output_artifacts(void)
{
  int t;

  if (SELF_PATH != NULL && SELF_ROOT != NULL)
    { if (SELF_HISTOGRAM_TOUCHED)
        cleanup_unlink(SELF_PATH,"%s.hist",SELF_ROOT);
      if (SELF_TABLE_TOUCHED)
        cleanup_unlink(SELF_PATH,"%s.ktab",SELF_ROOT);
      if (SELF_PROFILE_TOUCHED)
        cleanup_unlink(SELF_PATH,"%s.prof",SELF_ROOT);
      for (t = 1; t <= NTHREADS; t++)
        { if (SELF_TABLE_TOUCHED)
            cleanup_unlink(SELF_PATH,".%s.ktab.%d",SELF_ROOT,t);
          if (SELF_PROFILE_TOUCHED)
            { cleanup_unlink(SELF_PATH,".%s.pidx.%d",SELF_ROOT,t);
              cleanup_unlink(SELF_PATH,".%s.prof.%d",SELF_ROOT,t);
            }
        }
    }

  if (RELATIVE_PROFILE_TOUCHED && REL_PATH != NULL && DUAL_ROOT != NULL)
    { cleanup_unlink(REL_PATH,"%s.prof",DUAL_ROOT);
      for (t = 1; t <= NTHREADS; t++)
        { cleanup_unlink(REL_PATH,".%s.pidx.%d",DUAL_ROOT,t);
          cleanup_unlink(REL_PATH,".%s.prof.%d",DUAL_ROOT,t);
        }
    }
}

static void cleanup_sort_root(char *root, int remove_table_parts)
{
  int p, t, panel;
  int input_threads;

  if (SORT_PATH == NULL || root == NULL)
    return;

  input_threads = ITHREADS;
  if (input_threads <= 0 || input_threads > NTHREADS)
    input_threads = NTHREADS;

  for (t = 0; t < input_threads; t++)
    cleanup_unlink(SORT_PATH,"%s.NS.T%d",root,t);

  for (p = 0; p < NPARTS; p++)
    { for (t = 0; t < input_threads; t++)
        { if (remove_table_parts)
            cleanup_unlink(SORT_PATH,"%s.%d.T%d",root,p,t);
          for (panel = 0; panel < NPANELS; panel++)
            cleanup_unlink(SORT_PATH,"%s.%d.P%d.%d",root,p,t,panel);
        }

      if (remove_table_parts)
        { for (t = 0; t < NTHREADS; t++)
            cleanup_unlink(SORT_PATH,"%s.%d.L%d",root,p,t);
          cleanup_unlink(SORT_PATH,"%s.U%d.ktab",root,p);
          for (t = 1; t <= NTHREADS; t++)
            cleanup_unlink(SORT_PATH,".%s.U%d.ktab.%d",root,p,t);
        }
    }
}

void Clean_Exit(int status)
{
  char message[128];
  int  length;

  /*
   * Keep the lock through _exit: if several workers report an error, only the first one
   * performs cleanup and the others cannot terminate the process underneath it.
   */
  if (pthread_mutex_lock(&Cleanup_Lock) != 0)
    _exit(status == 0 ? 1 : status);

  length = snprintf(message,sizeof(message),"\n*** FastKdual error exit %d; cleaning partial outputs ***\n",
                    status);
  if (length > 0)
    { if ((size_t) length >= sizeof(message))
        length = sizeof(message)-1;
      write(STDERR_FILENO,message,(size_t) length);
    }

  cleanup_output_artifacts();
  cleanup_sort_root(SELF_ROOT,1);
  cleanup_sort_root(DUAL_ROOT,0);

  _exit(status == 0 ? 1 : status);
}

static void usage(void)
{
  fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
  fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
  fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  -p:<table>  reference table used for the relative profile (required)\n");
  fprintf(stderr,"  -R<root>    output root for the relative profile (required)\n");
  fprintf(stderr,"  -N<root>    output root for self histogram, table, and profile\n");
  fprintf(stderr,"  -t<int>     self table cutoff (default 1)\n");
  fprintf(stderr,"  -T<int>     number of threads\n");
  fprintf(stderr,"  -P<dir>     directory for temporary partitions and sorts\n");
  fprintf(stderr,"  -M<int>     GB available to sorting\n");
  fprintf(stderr,"  -k<int>     k-mer size\n");
  fprintf(stderr,"  -bc<int>    ignore this prefix length in every sequence\n");
  fprintf(stderr,"  -c          homopolymer-compress input\n");
  fprintf(stderr,"  -v          verbose progress and resource reporting\n\n");
}

static void check_directory(char *path, char *option)
{
  DIR *dirp = opendir(path);
  if (dirp == NULL)
    { fprintf(stderr,"%s: %s cannot open directory %s: %s\n",
                     Prog_Name,option,path,strerror(errno));
      exit(1);
    }
  closedir(dirp);
}

static int safe_path_character(unsigned char c)
{
  return ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
          (c >= '0' && c <= '9') ||
          c == '/' || c == '.' || c == '_' || c == '-' ||
          c == '+' || c == ':' || c == '=' || c == ',' || c == '@' || c == '%');
}

static void check_safe_path(char *path, char *label)
{
  unsigned char *scan = (unsigned char *) path;

  if (path == NULL || path[0] == '\0')
    { fprintf(stderr,"%s: %s path is empty\n",Prog_Name,label);
      exit(1);
    }

  for (; *scan != '\0'; scan++)
    if (! safe_path_character(*scan))
      { if (*scan >= 0x20 && *scan <= 0x7e)
          fprintf(stderr,
                  "%s: %s path '%s' contains unsupported shell-sensitive character '%c'\n",
                  Prog_Name,label,path,*scan);
        else
          fprintf(stderr,
                  "%s: %s path contains unsupported whitespace/control byte 0x%02x\n",
                  Prog_Name,label,*scan);
        exit(1);
      }
}

static void check_safe_current_directory(void)
{
  char *cwd = getcwd(NULL,0);

  if (cwd == NULL)
    { fprintf(stderr,"%s: Cannot determine current directory\n",Prog_Name);
      exit(1);
    }
  check_safe_path(cwd,"current-directory");
  free(cwd);
}

static char *absolute_sort_path(char *path)
{
  char *cwd, *result;

  if (path[0] == '/')
    return Strdup(path,"Allocating sort path");

  cwd = getcwd(NULL,0);
  if (cwd == NULL)
    { fprintf(stderr,"%s: Cannot determine current directory\n",Prog_Name);
      exit(1);
    }

  if (strcmp(path,".") == 0)
    result = Strdup(cwd,"Allocating sort path");
  else if (strncmp(path,"./",2) == 0)
    result = Strdup(Catenate(cwd,"/",path+2,""),"Allocating sort path");
  else
    result = Strdup(Catenate(cwd,"/",path,""),"Allocating sort path");
  free(cwd);
  return result;
}

static int copy_file(char *source, char *target)
{
  int in, out;
  char *buf;
  ssize_t n;

  in = open(source,O_RDONLY);
  if (in < 0)
    return -1;
  out = open(target,O_WRONLY|O_CREAT|O_EXCL,0600);
  if (out < 0)
    { close(in);
      return -1;
    }

  buf = Malloc(1<<20,"Allocating NS copy buffer");
  if (buf == NULL)
    { close(in);
      close(out);
      return -1;
    }

  while ((n = read(in,buf,1<<20)) > 0)
    { ssize_t done = 0;
      while (done < n)
        { ssize_t w = write(out,buf+done,n-done);
          if (w <= 0)
            { free(buf);
              close(in);
              close(out);
              unlink(target);
              return -1;
            }
          done += w;
        }
    }

  free(buf);
  close(in);
  if (close(out) != 0 || n < 0)
    { unlink(target);
      return -1;
    }
  return 0;
}

static void clone_invalid_streams(char *self_root, char *relative_root)
{
  char *source, *target;
  int   t;
  size_t cap = 2*strlen(SORT_PATH)+strlen(self_root)+strlen(relative_root)+100;

  source = Malloc(cap,"Allocating NS source name");
  target = Malloc(cap,"Allocating NS target name");
  if (source == NULL || target == NULL)
    Clean_Exit(1);

  for (t = 0; t < ITHREADS; t++)
    { sprintf(source,"%s/%s.NS.T%d",SORT_PATH,self_root,t);
      sprintf(target,"%s/%s.NS.T%d",SORT_PATH,relative_root,t);

      if (unlink(target) != 0 && errno != ENOENT)
        { fprintf(stderr,"%s: Cannot replace %s: %s\n",Prog_Name,target,strerror(errno));
          Clean_Exit(1);
        }

      if (link(source,target) != 0)
        if (copy_file(source,target) != 0)
          { fprintf(stderr,"%s: Cannot clone %s as %s: %s\n",
                           Prog_Name,source,target,strerror(errno));
            Clean_Exit(1);
          }
    }

  free(target);
  free(source);
}

static int64 parse_integer_option(char option, char *text, int64 minimum, int64 maximum)
{
  char      *end;
  long long  value;

  errno = 0;
  value = strtoll(text,&end,10);
  if (text[0] == '\0' || *end != '\0' || errno == ERANGE ||
      value < minimum || value > maximum)
    { fprintf(stderr,
              "%s: -%c value '%s' is outside the supported range [%lld,%lld]\n",
              Prog_Name,option,text,(long long) minimum,(long long) maximum);
      exit(1);
    }
  return (int64) value;
}

static void ensure_file_limit(int parts)
{
  struct rlimit rlp;
  uint64 nfiles = ((uint64) parts+3)*(uint64) NTHREADS + 8;

  if (getrlimit(RLIMIT_NOFILE,&rlp) != 0)
    { fprintf(stderr,"%s: Cannot inspect open-file limit: %s\n",
                     Prog_Name,strerror(errno));
      Clean_Exit(1);
    }
  if (nfiles > rlp.rlim_max)
    { fprintf(stderr,"%s: Cannot open %llu files simultaneously\n",
                     Prog_Name,(unsigned long long) nfiles);
      Clean_Exit(1);
    }
  if (nfiles > rlp.rlim_cur)
    { rlp.rlim_cur = nfiles;
      if (setrlimit(RLIMIT_NOFILE,&rlp) != 0)
        { fprintf(stderr,"%s: Cannot raise open-file limit: %s\n",
                         Prog_Name,strerror(errno));
          Clean_Exit(1);
        }
    }
}

int main(int argc, char *argv[])
{
  char *out_name = NULL;
  char *relative_name = NULL;
  int   i, j;

  start_time();
  Prog_Name = Strdup("FastKdual","");

  KMER        = 40;
  SORT_MEMORY = 12000000000ll;
  NTHREADS    = 4;
  SORT_PATH   = getenv("TMPDIR");
  if (SORT_PATH == NULL)
    SORT_PATH = "/tmp";
  DO_TABLE    = 1;
  DO_PROFILE  = 1;
  PRO_TABLE   = NULL;
  PRO_NAME    = NULL;
  BC_PREFIX   = 0;
  COMPRESS    = 0;

  j = 1;
  for (i = 1; i < argc; i++)
    if (argv[i][0] != '-')
      argv[j++] = argv[i];
    else
      { int64 value;

        if (strcmp(argv[i],"-v") == 0)
          VERBOSE = 1;
        else if (strcmp(argv[i],"-c") == 0)
          COMPRESS = 1;
        else if (strncmp(argv[i],"-bc",3) == 0)
          { value = parse_integer_option('b',argv[i]+3,0,INT_MAX);
            BC_PREFIX = value;
          }
        else if (strncmp(argv[i],"-p:",3) == 0 && argv[i][3] != '\0')
          { if (PRO_NAME != NULL)
              { fprintf(stderr,"%s: -p:<table> may be specified only once\n",Prog_Name);
                exit(1);
              }
            PRO_NAME = argv[i]+3;
          }
        else if (strncmp(argv[i],"-R",2) == 0 && argv[i][2] != '\0')
          relative_name = argv[i]+2;
        else if (strncmp(argv[i],"-N",2) == 0 && argv[i][2] != '\0')
          out_name = argv[i]+2;
        else if (strncmp(argv[i],"-P",2) == 0 && argv[i][2] != '\0')
          SORT_PATH = argv[i]+2;
        else if ((argv[i][1] == 'k' || argv[i][1] == 't' ||
                  argv[i][1] == 'M' || argv[i][1] == 'T') &&
                 argv[i][2] != '\0')
          { if (argv[i][1] == 'M')
              value = parse_integer_option('M',argv[i]+2,1,INT64_MAX/1000000000ll);
            else
              value = parse_integer_option(argv[i][1],argv[i]+2,1,
                                           argv[i][1] == 'k' ? FASTKDUAL_MAX_KMER
                                                             : INT_MAX);
            if (argv[i][1] == 'k')
              KMER = value;
            else if (argv[i][1] == 't')
              DO_TABLE = value;
            else if (argv[i][1] == 'M')
              SORT_MEMORY = value*1000000000ll;
            else
              NTHREADS = value;
          }
        else
          { fprintf(stderr,"%s: illegal option %s\n",Prog_Name,argv[i]);
            usage();
            exit(1);
          }
      }
  argc = j;

  if (argc < 2 || PRO_NAME == NULL || relative_name == NULL)
    { usage();
      exit(1);
    }

  /*
   * Vendored FastK io/table/merge modules contain shell-based legacy cleanup/decompression.
   * Validate every user-controlled path before Open_Kmer_Stream or Partition_Input can reach
   * those modules.
   */
  check_safe_current_directory();
  check_safe_path(SORT_PATH,"-P");
  check_safe_path(PRO_NAME,"-p");
  check_safe_path(relative_name,"-R");
  if (out_name != NULL)
    check_safe_path(out_name,"-N");
  for (i = 1; i < argc; i++)
    check_safe_path(argv[i],"source");

  // Partition_Input allocates NTHREADS-sized control arrays, so reject a
  // thread count that cannot support even one sort partition before entering it.
  ensure_file_limit(1);

  PRO_TABLE = Open_Kmer_Stream(PRO_NAME);
  if (PRO_TABLE == NULL)
    { fprintf(stderr,"%s: Cannot open FastK table %s\n",Prog_Name,PRO_NAME);
      exit(1);
    }
  if (PRO_TABLE->kmer != KMER)
    { fprintf(stderr,"%s: reference table k-mer size (%d) != requested k-mer size (%d)\n",
                     Prog_Name,PRO_TABLE->kmer,KMER);
      exit(1);
    }

  SORT_PATH = absolute_sort_path(SORT_PATH);
  check_safe_path(SORT_PATH,"absolute -P");
  check_directory(SORT_PATH,"-P");

  SELF_PATH = SELF_ROOT = REL_PATH = DUAL_ROOT = NULL;
  IOPACK = Partition_Input(argc,argv);
  if (IOPACK == NULL)
    exit(1);

  if (out_name == NULL)
    { SELF_ROOT = First_Root(IOPACK);
      SELF_PATH = First_Pwd(IOPACK);
    }
  else
    { SELF_ROOT = Root(out_name,"");
      SELF_PATH = PathTo(out_name);
    }
  DUAL_ROOT = Root(relative_name,"");
  REL_PATH  = PathTo(relative_name);

  check_safe_path(SELF_PATH,"self output directory");
  check_safe_path(SELF_ROOT,"self output root");
  check_safe_path(REL_PATH,"relative output directory");
  check_safe_path(DUAL_ROOT,"relative output root");
  check_directory(SELF_PATH,"-N");
  check_directory(REL_PATH,"-R");
  if (strcmp(SELF_ROOT,DUAL_ROOT) == 0)
    { fprintf(stderr,
              "%s: self and relative output basenames must differ because they share -P (%s)\n",
              Prog_Name,SELF_ROOT);
      exit(1);
    }

  { DATA_BLOCK *block;
    int64       gsize;
    int         rsize, val;

    if (VERBOSE)
      fprintf(stderr,"\nDetermining minimizer scheme & partition for %s\n",SELF_ROOT);

    block = Get_First_Block(IOPACK,1000000000);
    KMER_BYTES = (KMER*2+7) >> 3;
    rsize = KMER_BYTES+2;
    gsize = block->totlen - KMER*block->nreads;
    if (gsize < block->totlen/3)
      fprintf(stderr,"%s: Warning: sequences average less than 1.5x k-mer size\n",Prog_Name);
    gsize = gsize*block->ratio*rsize;
    NPARTS = (gsize-1)/SORT_MEMORY + 1;

    // Recheck with the data-dependent partition count before scheme allocation.
    ensure_file_limit(NPARTS);

    MOD_LEN = 1;
    while (MOD_LEN < KMER)
      MOD_LEN <<= 1;
    MOD_LEN <<= 1;
    MOD_MSK = MOD_LEN-1;

    MAX_SUPER = Determine_Scheme(block);
    Free_First_Block(block);

    SMER = MAX_SUPER + KMER - 1;
    SLEN_BITS = 0;
    for (val = MAX_SUPER; val > 0; val >>= 1)
      SLEN_BITS += 1;
    SLEN_BIT_MASK = (0x1u << SLEN_BITS)-1;
    SLEN_BYTES = (SLEN_BITS+7) >> 3;
    SMER_BYTES = (SMER*2+7) >> 3;
    SMER_WORD  = SMER_BYTES+SLEN_BYTES;
    KMER_WORD  = KMER_BYTES+2;
    PLEN_BYTES = (SLEN_BITS+8) >> 3;
    TMER_WORD  = KMER_BYTES+2;
  }

  Split_Kmers(IOPACK,SELF_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  clone_invalid_streams(SELF_ROOT,DUAL_ROOT);
  Split_Table(SELF_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  Free_Input_Partition(IOPACK);
  IOPACK = NULL;

  Sorting_Dual(SELF_PATH,SELF_ROOT,DUAL_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  SELF_TABLE_TOUCHED = 1;
  Merge_Tables(SELF_PATH,SELF_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  SELF_PROFILE_TOUCHED = 1;
  Merge_Profiles(SELF_PATH,SELF_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  RELATIVE_PROFILE_TOUCHED = 1;
  Merge_Profiles(REL_PATH,DUAL_ROOT);
  if (VERBOSE)
    time_to(stderr,0);

  free(NUM_RID);
  Free_Kmer_Stream(PRO_TABLE);
  free(REL_PATH);
  free(DUAL_ROOT);
  free(SELF_PATH);
  free(SELF_ROOT);
  free(SORT_PATH);
  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  if (VERBOSE)
    time_to(stderr,1);

  return 0;
}
