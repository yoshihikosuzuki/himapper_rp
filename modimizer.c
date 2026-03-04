#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "DB.h"
#include "nthash.h"

// #define DEBUG_HASH

const int K = 32;
const int Km1 = K - 1;
const int MOD = 1001;
const int MOD_THRES = 1;

// ASCII base → 2-bit value  (A=0, C=1, G=2, T=3), matching the DAZZ_DB
// numeric encoding used by himapper's tuple_thread.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

// Bitmask covering 2*K bits.  For K >= 32 the full k-mer string does not fit
// in a uint64_t; we keep only the lower 64 bits (last 32 bases), consistent
// with how himapper handles K > 32.
static const uint64_t KMASK =
    (K >= 32) ? UINT64_MAX : ((UINT64_C(1) << (2 * K)) - 1);

// Compute canonical ntHash for every k-mer in the read, select seeds whose
// (hVal % MOD) < MOD_THRES, and write each selected entry as
//   { uint64_t code,  int32_t rpos,  int32_t rid }
// to `out`.  The `code` field holds the 2-bit encoding of the k-mer string
// (like himapper's tuple_thread), masked to 64 bits.
//
//   seq  – pointer to position (K-1) in the read; seq[-(K-1)] is the first
//           base and seq[0] is the last base of the first k-mer.
//   plen – number of k-mers (= rlen - K + 1).
//
// Returns the number of entries written.
int find_seeds(const char *seq, const int plen, const int K,
               const int read_id, FILE *out) {
  uint64_t hVal, fhVal = 0, rhVal = 0;
  int written = 0;

  // First k-mer: starts at seq − K + 1.
  hVal = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);

  if (hVal % MOD < MOD_THRES) {
    const char *kp = seq - K + 1;
    uint64_t c = 0;
    for (int x = 0; x < K; x++)
      c = (c << 2) | base2bit[(unsigned char)kp[x]];
    c &= KMASK;
    int32_t rpos = 0, rid = read_id;
    fwrite(&c,    sizeof(uint64_t), 1, out);
    fwrite(&rpos, sizeof(int32_t),  1, out);
    fwrite(&rid,  sizeof(int32_t),  1, out);
    written++;
#ifdef DEBUG_HASH
    { char cstr[65]; uint64_t tmp = c;
      for (int b = K - 1; b >= 0; b--) { cstr[b] = "ACGT"[tmp & 3]; tmp >>= 2; }
      cstr[K] = '\0';
      fprintf(stderr,
              "  seed: read=%d rpos=%d kmer=%.*s code=%s hVal=%lu f=%lu r=%lu\n",
              rid, rpos, K, kp, cstr, hVal, fhVal, rhVal);
    }
#endif
  }

  for (int i = 1; i < plen; i++) {
    // outgoing base: seq[i − K],  incoming base: seq[i]
    hVal = NTC64_c((unsigned char)seq[i - K], (unsigned char)seq[i],
                   K, &fhVal, &rhVal);

    if (hVal % MOD < MOD_THRES) {
      const char *kp = seq - K + 1 + i;
      uint64_t c = 0;
      for (int x = 0; x < K; x++)
        c = (c << 2) | base2bit[(unsigned char)kp[x]];
      c &= KMASK;
      int32_t rpos = (int32_t)i, rid = read_id;
      fwrite(&c,    sizeof(uint64_t), 1, out);
      fwrite(&rpos, sizeof(int32_t),  1, out);
      fwrite(&rid,  sizeof(int32_t),  1, out);
      written++;
#ifdef DEBUG_HASH
      { char cstr[65]; uint64_t tmp = c;
        for (int b = K - 1; b >= 0; b--) { cstr[b] = "ACGT"[tmp & 3]; tmp >>= 2; }
        cstr[K] = '\0';
        fprintf(stderr,
                "  seed: read=%d rpos=%d kmer=%.*s code=%s hVal=%lu f=%lu r=%lu\n",
                rid, rpos, K, kp, cstr, hVal, fhVal, rhVal);
      }
#endif
    }
  }

  return written;
}

int main(int argc, char *argv[]) {
  DAZZ_DB _db, *db = &_db;
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <db> <out.kpos>\n", argv[0]);
    return 1;
  }
  Open_DB(argv[1], db);

  FILE *out = fopen(argv[2], "wb");
  if (!out) {
    perror("fopen");
    return 1;
  }

  // Write simple header placeholder (24 bytes):
  // magic(4) ver(1) flags(1) idx_bytes(1) reserved(1) K(int32) entry_size(int32) n_entries(uint64)
  char magic[4] = {'K', 'P', 'O', 'S'};
  fwrite(magic, 4, 1, out);
  uint8_t ver = 1, flags = 0, idx_bytes = 0, reserved = 0;
  fwrite(&ver, 1, 1, out);
  fwrite(&flags, 1, 1, out);
  fwrite(&idx_bytes, 1, 1, out);
  fwrite(&reserved, 1, 1, out);
  int32_t Kval = K;
  fwrite(&Kval, sizeof(int32_t), 1, out);
  int32_t entry_size = (int32_t)(sizeof(uint64_t) + sizeof(int32_t) + sizeof(int32_t));
  fwrite(&entry_size, sizeof(int32_t), 1, out);
  uint64_t n_entries = 0;
  fwrite(&n_entries, sizeof(uint64_t), 1, out);

  {
    char *seq = New_Read_Buffer(db);

    uint64_t total_written = 0;
    for (int id = 0; id < db->nreads; id++) {
      DAZZ_READ *r = db->reads + id;
      int rlen = r->rlen;
      int plen = rlen - Km1;

      Load_Read(db, id, seq, 2);
    //   printf("%s\n", seq);
      int written = find_seeds(seq + Km1, plen, K, id, out);
      total_written += (uint64_t)written;
    }

    // Seek back and write the total entry count into header (offset 16)
    if (fseek(out, 16, SEEK_SET) == 0) {
      fwrite(&total_written, sizeof(uint64_t), 1, out);
    }
    fprintf(stderr, "Seeds written: %llu\n", (unsigned long long)total_written);

    free(seq-1);
  }

  fclose(out);

  return 0;
}
