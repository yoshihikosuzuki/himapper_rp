#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "DB.h"
#include "nthash.h"
#include "libfastk.h"  // FASTK profile API (compile with -IFASTK)

// #define DEBUG_HASH

const int K = 32;
const int Km1 = K - 1;
#define WSIZE 2000   // sliding window width for (count, hash) minimizer
#define MOD   1001   // hash modulus for tiebreaking (prime > WSIZE recommended)

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

// Compare (ca, ha) < (cb, hb) lexicographically: count first, hash to break ties.
static inline int tuple_lt(uint16_t ca, uint64_t ha, uint16_t cb, uint64_t hb)
{ return (ca < cb) || (ca == cb && ha < hb); }

// Select seeds using a sliding-window (count, hash) minimizer.
//
// For every window of WSIZE consecutive k-mers the one with the lexicographically
// smallest (profile[i], hVal[i]) tuple is selected.  The result is the union of
// all per-window selections.
//
//   seq     – pointer to position (K-1) in the read (same convention as before)
//   plen    – number of k-mers (rlen - K + 1)
//   K       – k-mer length
//   profile – FastK uint16 count for positions 0..plen-1
//   read_id – 0-based read index
//   out     – binary output file
//
// Returns the number of entries written.
int find_seeds(const char    *seq,
               const int      plen,
               const int      K,
               const uint16  *profile,
               const int      read_id,
               FILE          *out)
{
  if (plen <= 0) return 0;

  // Step 1: compute all canonical hash values.
  uint64_t *hvals = (uint64_t *) malloc((size_t)plen * sizeof(uint64_t));
  if (!hvals) { perror("malloc hvals"); exit(1); }
  { uint64_t fhVal = 0, rhVal = 0;
    hvals[0] = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);
    for (int i = 1; i < plen; i++)
      hvals[i] = NTC64_c((unsigned char)seq[i - K], (unsigned char)seq[i],
                         K, &fhVal, &rhVal);
    // Reduce to [0, 1000] so that hash ties are common — enables robust winnowing
    // to reliably favour the leftmost (earliest) occurrence within a window.
    for (int i = 0; i < plen; i++)
      hvals[i] %= MOD;
  }

  // Step 2: sliding-window minimum with a monotone deque.
  //   The deque stores indices in increasing position order.
  //   Front = index with the minimum (count, hash) in the current window.
  //   We pop from the back any index j whose tuple >= the incoming tuple
  //   (it can never win a future window while the new element is around).
  int  *dq      = (int  *) malloc((size_t)plen * sizeof(int));
  bool *is_seed = (bool *) calloc((size_t)plen,  sizeof(bool));
  if (!dq || !is_seed) { perror("malloc deque/is_seed"); exit(1); }

  int dq_head = 0, dq_tail = 0;
  for (int i = 0; i < plen; i++)
    { // Robust winnowing: evict back elements whose tuple is STRICTLY GREATER
      // than the current tuple.  On a tie we keep the existing (leftmost) entry,
      // so ties are broken in favour of the element already inside the window.
      while (dq_head < dq_tail &&
             tuple_lt(profile[i],             hvals[i],
                      profile[dq[dq_tail-1]], hvals[dq[dq_tail-1]]))
        dq_tail--;
      dq[dq_tail++] = i;

      // Evict front elements that have slid out of the window.
      while (dq[dq_head] < i - WSIZE + 1)
        dq_head++;

      // Only consider complete windows of exactly WSIZE k-mers (standard minimizer definition).
      if (i >= WSIZE - 1)
        is_seed[dq[dq_head]] = true;
    }

  // Step 3: write selected seeds.
  int written = 0;
  for (int i = 0; i < plen; i++)
    { if (!is_seed[i]) continue;
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
                "  seed: read=%d rpos=%d kmer=%.*s code=%s count=%u hash=%lu\n",
                rid, rpos, K, kp, cstr, (unsigned)profile[i], hvals[i]);
      }
#endif
    }

  free(is_seed);
  free(dq);
  free(hvals);
  return written;
}

int main(int argc, char *argv[]) {
  DAZZ_DB _db, *db = &_db;
  if (argc < 4) {
    fprintf(stderr, "Usage: %s <db> <profiles>[.prof] <out.kpos>\n", argv[0]);
    fprintf(stderr, "  <db>       : DAZZ reference dam/db to build k-mer index from\n");
    fprintf(stderr, "  <profiles> : .prof file from FastK/Profex\n");
    fprintf(stderr, "  <out.kpos> : output k-mer position file\n");
    return 1;
  }
  Open_DB(argv[1], db);

  // Open FastK profile (argv[2]). Must be opened before the main loop.
  Profile_Index *P = Open_Profiles(argv[2]);
  if (P == NULL)
    { fprintf(stderr, "%s: Cannot open profile file '%s'\n", argv[0], argv[2]);
      return 1;
    }
  if (P->kmer != K)
    { fprintf(stderr, "%s: profile K=%d does not match K=%d\n", argv[0], P->kmer, K);
      Free_Profiles(P);
      return 1;
    }
  fprintf(stderr, "Profile: K=%d  nreads=%lld\n",
          P->kmer, (long long)P->nbase[P->nparts-1]);

  int     pmax    = 20000;
  uint16 *profile = (uint16 *) malloc((size_t)pmax * sizeof(uint16));
  if (!profile) { perror("malloc profile"); return 1; }

  FILE *out = fopen(argv[3], "wb");
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

      // Fetch k-mer count profile for this read (0-based id).
      int fetched = Fetch_Profile(P, (int64)id, pmax, profile);
      if (fetched > pmax)
        { pmax    = (int)(1.2 * fetched) + 1000;
          profile = (uint16 *) realloc(profile, (size_t)pmax * sizeof(uint16));
          if (!profile) { perror("realloc profile"); exit(1); }
          Fetch_Profile(P, (int64)id, pmax, profile);
        }

      int written = find_seeds(seq + Km1, plen, K, profile, id, out);
      total_written += (uint64_t)written;
    }

    // Seek back and write the total entry count into header (offset 16)
    if (fseek(out, 16, SEEK_SET) == 0) {
      fwrite(&total_written, sizeof(uint64_t), 1, out);
    }
    fprintf(stderr, "Seeds written: %llu\n", (unsigned long long)total_written);

    free(seq-1);
  }

  free(profile);
  Free_Profiles(P);

  fclose(out);

  return 0;
}
