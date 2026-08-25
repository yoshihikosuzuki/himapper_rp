/*  dump_kpos  –  Print a block-indexed .kpos container as ASCII.
 *
 *  Usage:  dump_kpos <file.kpos> [block]
 *
 *  Header layout (24 bytes):
 *    magic[4]     "KPOS"
 *    ver          uint8
 *    flags        uint8
 *    idx_bytes    uint8
 *    reserved     uint8
 *    K            int32
 *    entry_size   int32
 *    n_entries    uint64
 *
 *  Each entry (16 bytes):
 *    code   uint64   – 2-bit encoded k-mer (lower 2*K bits)
 *    rpos   int32    – position within the read/contig
 *    rid    uint32   – block-local read/contig id (0-based) with the genome
 *                      copy-tier packed into the top 3 bits (see copynum.h)
 *
 *  The sidecar <file.kpos>.idx is an array of uint64 entry offsets, one
 *  for each block boundary.
 *
 *  Output columns:
 *    block(1-based)  local_read_id(1-based)  rpos  kmer  tier
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "copynum.h"

int main(int argc, char *argv[])
{
  uint64_t selected_block = 0;
  FILE    *fp, *idx;
  char    *idx_path;

  if (argc < 2 || argc > 3)
    { fprintf(stderr, "Usage: %s <file.kpos> [block]\n", argv[0]);
      return 1;
    }

  if (argc == 3)
    { char *end;
      long long value = strtoll(argv[2], &end, 10);
      if (*argv[2] == '\0' || *end != '\0' || value < 1)
        { fprintf(stderr, "Error: block must be a positive integer\n");
          return 1;
        }
      selected_block = (uint64_t)value;
    }

  fp = fopen(argv[1], "rb");
  if (fp == NULL)
    { perror(argv[1]); return 1; }

  char     magic[4];
  uint8_t  ver, flags, idx_bytes, reserved;
  int32_t  K, entry_size;
  uint64_t n_entries;

  if (fread(magic, 4, 1, fp) != 1 ||
      fread(&ver, 1, 1, fp) != 1 ||
      fread(&flags, 1, 1, fp) != 1 ||
      fread(&idx_bytes, 1, 1, fp) != 1 ||
      fread(&reserved, 1, 1, fp) != 1 ||
      fread(&K, sizeof(int32_t), 1, fp) != 1 ||
      fread(&entry_size, sizeof(int32_t), 1, fp) != 1 ||
      fread(&n_entries, sizeof(uint64_t), 1, fp) != 1)
    { fprintf(stderr, "Error: truncated header\n"); fclose(fp); return 1; }

  if (memcmp(magic, "KPOS", 4) != 0)
    { fprintf(stderr, "Error: bad magic (expected KPOS)\n"); fclose(fp); return 1; }
  if (ver != 2 || idx_bytes != sizeof(uint64_t) ||
      entry_size != 16 || K < 1 || K > 32)
    { fprintf(stderr, "Error: unsupported KPOS format\n"); fclose(fp); return 1; }

  off_t data_start = ftello(fp);
  size_t idx_path_len = strlen(argv[1]) + strlen(".idx") + 1;
  idx_path = (char *)malloc(idx_path_len);
  if (idx_path == NULL)
    { perror("malloc"); fclose(fp); return 1; }
  snprintf(idx_path, idx_path_len, "%s.idx", argv[1]);

  idx = fopen(idx_path, "rb");
  if (idx == NULL)
    { perror(idx_path); free(idx_path); fclose(fp); return 1; }
  if (fseeko(idx, 0, SEEK_END) != 0)
    { perror(idx_path); free(idx_path); fclose(idx); fclose(fp); return 1; }

  off_t idx_size = ftello(idx);
  if (idx_size < (off_t)(2 * sizeof(uint64_t)) ||
      idx_size % (off_t)sizeof(uint64_t) != 0)
    { fprintf(stderr, "Error: malformed KPOS index %s\n", idx_path);
      free(idx_path); fclose(idx); fclose(fp); return 1;
    }

  uint64_t nblocks =
      (uint64_t)(idx_size / (off_t)sizeof(uint64_t)) - 1;
  uint64_t *offsets =
      (uint64_t *)malloc((size_t)(nblocks + 1) * sizeof(uint64_t));
  if (offsets == NULL)
    { perror("malloc"); free(idx_path); fclose(idx); fclose(fp); return 1; }

  rewind(idx);
  if (fread(offsets, sizeof(uint64_t), (size_t)nblocks + 1, idx) !=
      (size_t)nblocks + 1)
    { fprintf(stderr, "Error: truncated KPOS index %s\n", idx_path);
      free(offsets); free(idx_path); fclose(idx); fclose(fp); return 1;
    }
  fclose(idx);

  if (offsets[0] != 0 || offsets[nblocks] != n_entries)
    { fprintf(stderr, "Error: KPOS data/index totals disagree\n");
      free(offsets); free(idx_path); fclose(fp); return 1;
    }
  for (uint64_t block = 0; block < nblocks; block++)
    if (offsets[block] > offsets[block + 1])
      { fprintf(stderr, "Error: KPOS block offsets are not monotone\n");
        free(offsets); free(idx_path); fclose(fp); return 1;
      }
  if (selected_block > nblocks)
    { fprintf(stderr, "Error: block %llu is outside 1-%llu\n",
              (unsigned long long)selected_block,
              (unsigned long long)nblocks);
      free(offsets); free(idx_path); fclose(fp); return 1;
    }

  printf("# magic=KPOS  ver=%u  flags=%u  idx_bytes=%u  K=%d  "
         "entry_size=%d  n_entries=%llu  nblocks=%llu\n",
         ver, flags, idx_bytes, K, entry_size,
         (unsigned long long)n_entries, (unsigned long long)nblocks);
  printf("# block\tlocal_read_id\trpos\tkmer\ttier\n");

  uint64_t first_block = selected_block ? selected_block : 1;
  uint64_t last_block  = selected_block ? selected_block : nblocks;
  for (uint64_t block = first_block; block <= last_block; block++)
    { uint64_t first = offsets[block - 1];
      uint64_t last  = offsets[block];
      off_t position =
          data_start + (off_t)(first * (uint64_t)entry_size);

      if (fseeko(fp, position, SEEK_SET) != 0)
        { perror(argv[1]); free(offsets); free(idx_path); fclose(fp); return 1; }

      for (uint64_t e = first; e < last; e++)
        { uint64_t code;
          int32_t  rpos;
          uint32_t rid;
          if (fread(&code, sizeof(uint64_t), 1, fp) != 1 ||
              fread(&rpos, sizeof(int32_t), 1, fp) != 1 ||
              fread(&rid, sizeof(uint32_t), 1, fp) != 1)
            { fprintf(stderr, "Error: truncated entry at %llu\n",
                      (unsigned long long)e);
              free(offsets); free(idx_path); fclose(fp); return 1;
            }

          char kmer[33];
          uint64_t tmp = code;
          for (int b = K - 1; b >= 0; b--)
            { kmer[b] = "ACGT"[tmp & 3];
              tmp >>= 2;
            }
          kmer[K] = '\0';

          printf("%llu\t%d\t%d\t%s\t%u\n",
                 (unsigned long long)block, CN_VAL(rid) + 1,
                 rpos, kmer, CN_TIER(rid));
        }
    }

  free(offsets);
  free(idx_path);
  fclose(fp);
  return 0;
}
