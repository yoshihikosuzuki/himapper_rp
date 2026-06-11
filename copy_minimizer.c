#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DB.h"
#include "libfastk.h"  // FASTK profile API (compile with -IFASTK)
#include "nthash.h"

// #define DEBUG_HASH
// #define DEBUG_NW

const int K = 32;
const int Km1 = K - 1;
#define WSIZE 100   // sliding window width for (count, hash) minimizer
#define MOD 2003    // hash modulus for tiebreaking (prime > WSIZE recommended)
#define RSPACE 100  // minimum spacing between uniformly sparse seeds

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

// Compare (ca, ha) < (cb, hb) lexicographically: count first, hash to break
// ties.
static inline int tuple_lt(uint16_t ca, uint64_t ha, uint16_t cb, uint64_t hb) {
  return (ca < cb) || (ca == cb && ha < hb);
}

// Candidate for uniformly sparse thinning.
typedef struct {
  int pos;
  int nw_val;
  uint16_t g_count;
  uint64_t hash;
} Candidate;

// Compare candidates by descending h_nw priority:
//   1. nw      descending (more windows → higher priority)
//   2. g_count ascending  (rarer in genome → better seed)
//   3. hash    ascending  (stable tiebreaker)
static int cand_cmp(const void* a, const void* b) {
  const Candidate* ca = (const Candidate*)a;
  const Candidate* cb = (const Candidate*)b;
  if (ca->nw_val != cb->nw_val) return (cb->nw_val > ca->nw_val) ? 1 : -1;
  if (ca->g_count != cb->g_count) return (ca->g_count > cb->g_count) ? 1 : -1;
  if (ca->hash != cb->hash) return (ca->hash > cb->hash) ? 1 : -1;
  return 0;
}

// Select seeds using a sliding-window (count, hash) minimizer.
//
// For every window of WSIZE consecutive k-mers the one with the
// lexicographically smallest (profile[i], hVal[i]) tuple is selected.  The
// result is the union of all per-window selections.
//
//   seq     – pointer to position (K-1) in the read (same convention as before)
//   plen    – number of k-mers (rlen - K + 1)
//   K       – k-mer length
//   profile – FastK uint16 count for positions 0..plen-1
//   read_id – 0-based read index
//   out     – binary output file
//
// Returns the number of entries written.
int find_seeds(const char* seq, const int plen, const int K,
               const uint16* profile, const int read_id, FILE* out) {
  if (plen <= 0) return 0;

  // Step 1: compute all canonical hash values.
  uint64_t* hvals = (uint64_t*)malloc((size_t)plen * sizeof(uint64_t));
  if (!hvals) {
    perror("malloc hvals");
    exit(1);
  }
  {
    uint64_t fhVal = 0, rhVal = 0;
    hvals[0] = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);
    for (int i = 1; i < plen; i++)
      hvals[i] = NTC64_c((unsigned char)seq[i - K], (unsigned char)seq[i], K,
                         &fhVal, &rhVal);
    // Reduce to [0, 1000] so that hash ties are common — enables robust
    // winnowing to reliably favour the leftmost (earliest) occurrence within a
    // window.
    for (int i = 0; i < plen; i++) hvals[i] %= MOD;
  }

  // Step 2: sliding-window minimum with a monotone deque.
  //   The deque stores indices in increasing position order.
  //   Front = index with the minimum (count, hash) in the current window.
  //   We pop from the back any index j whose tuple >= the incoming tuple
  //   (it can never win a future window while the new element is around).
  int* dq = (int*)malloc((size_t)plen * sizeof(int));
  bool* is_seed = (bool*)calloc((size_t)plen, sizeof(bool));
  int* nw = (int*)calloc(
      (size_t)plen,
      sizeof(int));  // # of windows for which position i is minimizer
  if (!dq || !is_seed || !nw) {
    perror("malloc deque/is_seed/nw");
    exit(1);
  }

  int dq_head = 0, dq_tail = 0;
  for (int i = 0; i < plen; i++) {  // Robust winnowing: evict back elements
                                    // whose tuple is STRICTLY GREATER
    // than the current tuple.  On a tie we keep the existing (leftmost) entry,
    // so ties are broken in favour of the element already inside the window.
    while (dq_head < dq_tail &&
           tuple_lt(profile[i], hvals[i], profile[dq[dq_tail - 1]],
                    hvals[dq[dq_tail - 1]]))
      dq_tail--;
    dq[dq_tail++] = i;

    // Evict front elements that have slid out of the window.
    while (dq[dq_head] < i - WSIZE + 1) dq_head++;

    // Only consider complete windows of exactly WSIZE k-mers (standard
    // minimizer definition).
    if (i >= WSIZE - 1) {
      is_seed[dq[dq_head]] = true;
      nw[dq[dq_head]]++;
    }
#ifdef DEBUG_HASH
    if (i >= WSIZE - 1) {
      int win_start = i - WSIZE + 1;
      int min_pos = dq[dq_head];
      fprintf(stderr, "WINDOW read=%d [%d..%d]  minimizer=@%d g=%u h=%lu\n",
              read_id, win_start, i, min_pos, (unsigned)profile[min_pos],
              hvals[min_pos]);
      for (int j = win_start; j <= i; j++) {
        const char* wkp = seq - K + 1 + j;
        char wkstr[65];
        for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
        wkstr[K] = '\0';
        fprintf(stderr, "  j=%-6d g=%-5u h=%-4lu  %s%s\n", j,
                (unsigned)profile[j], hvals[j], wkstr,
                (j == min_pos) ? "  <-- MIN" : "");
      }
    }
#endif
  }

  // Step 2.5: Uniformly sparse thinning.
  // Sort candidates by descending h_nw = (NW, -h_minimax) priority, then
  // greedily accept only those whose position is ≥ RSPACE from every
  // already-accepted seed.
  {
    int ncand = 0;
    for (int i = 0; i < plen; i++)
      if (is_seed[i]) ncand++;

    if (ncand > 1) {
      Candidate* cands = (Candidate*)malloc((size_t)ncand * sizeof(Candidate));
      if (!cands) {
        perror("malloc cands");
        exit(1);
      }
      {
        int ci = 0;
        for (int i = 0; i < plen; i++)
          if (is_seed[i]) {
            cands[ci].pos = i;
            cands[ci].nw_val = nw[i];
            cands[ci].g_count = profile[i];
            cands[ci].hash = hvals[i];
            ci++;
          }
      }
      qsort(cands, (size_t)ncand, sizeof(Candidate), cand_cmp);

      // Reset is_seed; re-mark only accepted positions.
      memset(is_seed, 0, (size_t)plen * sizeof(bool));

      // Maintain accepted positions in sorted order (by position) for
      // O(log n) nearest-neighbour distance checks.
      int* acc = (int*)malloc((size_t)ncand * sizeof(int));
      if (!acc) {
        perror("malloc acc");
        exit(1);
      }
      int nacc = 0;

      for (int ci = 0; ci < ncand; ci++) {
        int p = cands[ci].pos;
        // Binary search for insertion point in acc[].
        int lo = 0, hi = nacc;
        while (lo < hi) {
          int mid = (lo + hi) / 2;
          if (acc[mid] < p)
            lo = mid + 1;
          else
            hi = mid;
        }
        // Check left and right neighbours in accepted set.
        int ok = 1;
        if (lo > 0 && p - acc[lo - 1] < RSPACE) ok = 0;
        if (lo < nacc && acc[lo] - p < RSPACE) ok = 0;
        if (ok) {
          memmove(acc + lo + 1, acc + lo, (size_t)(nacc - lo) * sizeof(int));
          acc[lo] = p;
          nacc++;
          is_seed[p] = true;
        }
      }
      free(acc);
      free(cands);
    }
  }

#ifdef DEBUG_NW
  // Print per-position window-count; mark final selected seeds with *****.
  fprintf(stderr, "NW read=%d plen=%d\n", read_id, plen);
  for (int i = 0; i < plen; i++) {
    const char* wkp = seq - K + 1 + i;
    char wkstr[65];
    for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
    wkstr[K] = '\0';
    fprintf(stderr, "  i=%-6d nw=%-5d g=%-5u h=%-4lu  %s  %s%s\n", i, nw[i],
            (unsigned)profile[i], hvals[i], wkstr, nw[i] ? "*" : "",
            is_seed[i] ? "****" : "");
  }
#endif

  // Step 3: write selected seeds.
  int written = 0;
  for (int i = 0; i < plen; i++) {
    if (!is_seed[i]) continue;
    const char* kp = seq - K + 1 + i;
    uint64_t c = 0;
    for (int x = 0; x < K; x++) c = (c << 2) | base2bit[(unsigned char)kp[x]];
    c &= KMASK;
    int32_t rpos = (int32_t)i, rid = read_id;
    fwrite(&c, sizeof(uint64_t), 1, out);
    fwrite(&rpos, sizeof(int32_t), 1, out);
    fwrite(&rid, sizeof(int32_t), 1, out);
    written++;
  }

  free(is_seed);
  free(nw);
  free(dq);
  free(hvals);
  return written;
}

int main(int argc, char* argv[]) {
  DAZZ_DB _db, *db = &_db;
  if (argc < 4) {
    fprintf(stderr, "Usage: %s <db> <profiles>[.prof] <out.kpos>\n", argv[0]);
    fprintf(stderr,
            "  <db>       : DAZZ reference dam/db to build k-mer index from\n");
    fprintf(stderr, "  <profiles> : .prof file from FastK/Profex\n");
    fprintf(stderr, "  <out.kpos> : output k-mer position file\n");
    return 1;
  }
  Open_DB(argv[1], db);

  // Open FastK profile (argv[2]). Must be opened before the main loop.
  Profile_Index* P = Open_Profiles(argv[2]);
  if (P == NULL) {
    fprintf(stderr, "%s: Cannot open profile file '%s'\n", argv[0], argv[2]);
    return 1;
  }
  if (P->kmer != K) {
    fprintf(stderr, "%s: profile K=%d does not match K=%d\n", argv[0], P->kmer,
            K);
    Free_Profiles(P);
    return 1;
  }
  fprintf(stderr, "Profile: K=%d  nreads=%lld\n", P->kmer,
          (long long)P->nbase[P->nparts - 1]);

  int pmax = 20000;
  uint16* profile = (uint16*)malloc((size_t)pmax * sizeof(uint16));
  if (!profile) {
    perror("malloc profile");
    return 1;
  }

  FILE* out = fopen(argv[3], "wb");
  if (!out) {
    perror("fopen");
    return 1;
  }

  // Write simple header placeholder (24 bytes):
  // magic(4) ver(1) flags(1) idx_bytes(1) reserved(1) K(int32)
  // entry_size(int32) n_entries(uint64)
  char magic[4] = {'K', 'P', 'O', 'S'};
  fwrite(magic, 4, 1, out);
  uint8_t ver = 1, flags = 0, idx_bytes = 0, reserved = 0;
  fwrite(&ver, 1, 1, out);
  fwrite(&flags, 1, 1, out);
  fwrite(&idx_bytes, 1, 1, out);
  fwrite(&reserved, 1, 1, out);
  int32_t Kval = K;
  fwrite(&Kval, sizeof(int32_t), 1, out);
  int32_t entry_size =
      (int32_t)(sizeof(uint64_t) + sizeof(int32_t) + sizeof(int32_t));
  fwrite(&entry_size, sizeof(int32_t), 1, out);
  uint64_t n_entries = 0;
  fwrite(&n_entries, sizeof(uint64_t), 1, out);

  {
    char* seq = New_Read_Buffer(db);

    uint64_t total_written = 0;
    for (int id = 0; id < db->nreads; id++) {
      DAZZ_READ* r = db->reads + id;
      int rlen = r->rlen;
      int plen = rlen - Km1;

      Load_Read(db, id, seq, 2);

      // Fetch k-mer count profile for this read (0-based id).
      int fetched = Fetch_Profile(P, (int64)id, pmax, profile);
      if (fetched > pmax) {
        pmax = (int)(1.2 * fetched) + 1000;
        profile = (uint16*)realloc(profile, (size_t)pmax * sizeof(uint16));
        if (!profile) {
          perror("realloc profile");
          exit(1);
        }
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

    free(seq - 1);
  }

  free(profile);
  Free_Profiles(P);

  fclose(out);

  return 0;
}
