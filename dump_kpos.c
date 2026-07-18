/*  dump_kpos  –  Print .kpos file contents in human-readable ASCII.
 *
 *  Usage:  dump_kpos <file.kpos>
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
 *    rid    uint32   – read/contig id (0-based) with the genome copy-tier packed
 *                      into the top 3 bits (see copynum.h)
 *
 *  Output columns:  contig_id(1-based)  rpos  kmer  tier
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "copynum.h"   // CN_VAL / CN_TIER: split the rid into id + copy-tier

int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <file.kpos>\n", argv[0]);
    return 1;
  }

  FILE *fp = fopen(argv[1], "rb");
  if (!fp) { perror(argv[1]); return 1; }

  /* --- Read header (24 bytes) --- */
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

  printf("# magic=KPOS  ver=%u  flags=%u  idx_bytes=%u  K=%d  entry_size=%d  n_entries=%llu\n",
         ver, flags, idx_bytes, K, entry_size, (unsigned long long)n_entries);

  /* --- Decode entries --- */
  int Kd = (K > 64) ? 32 : K;   /* max bases we can decode from 64-bit code */

  uint64_t count = 0;
  for (uint64_t e = 0; e < n_entries; e++) {
    uint64_t code;
    int32_t  rpos;
    uint32_t rid;
    if (fread(&code, sizeof(uint64_t), 1, fp) != 1 ||
        fread(&rpos, sizeof(int32_t),  1, fp) != 1 ||
        fread(&rid,  sizeof(uint32_t), 1, fp) != 1)
      { fprintf(stderr, "Error: truncated entry at %llu\n", (unsigned long long)e);
        break;
      }

    /* Decode 2-bit code → ACGT string */
    char kmer[65];
    uint64_t tmp = code;
    for (int b = Kd - 1; b >= 0; b--) {
      kmer[b] = "ACGT"[tmp & 3];
      tmp >>= 2;
    }
    kmer[Kd] = '\0';

    printf("%d\t%d\t%s\t%u\n",
           CN_VAL(rid) + 1, rpos, kmer, CN_TIER(rid));
    count++;
  }

  if (count != n_entries)
    fprintf(stderr, "Warning: expected %llu entries, read %llu\n",
            (unsigned long long)n_entries, (unsigned long long)count);

  fclose(fp);
  return 0;
}
