#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "DB.h"
#include "nthash.h"
#include "libfastk.h"
#include "copynum.h"

// #define DEBUG_HASH
// #define DEBUG_NW

int K = 32;          // k-mer length; maximum 32   [-k]
int WSIZE  = 100;    // sliding window width for (count, hash) minimizer   [-w]
int RSPACE = 100;    // minimum spacing between uniformly sparse seeds     [-s]
static const int MOD = 2003;   // hash modulus for tiebreaking (prime > WSIZE recommended)

// ASCII base → 2-bit value  (A=0, C=1, G=2, T=3), matching the DAZZ_DB
// numeric encoding used by himapper's tuple_thread.
static const uint8_t base2bit[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
};

// Bitmask covering the low 2*K bits of a k-mer's 2-bit code; set in main() once K is
// known (K is capped at 32, so K==32 gives the full uint64 mask).
static uint64_t KMASK;

// Compare (ca_g, ca_r, ha) < (cb_g, cb_r, hb):
//   primary:   lower genome count  → more unique in genome   → better seed
//   secondary: higher reads count  → better covered in reads → better seed
//   tertiary:  lower hash % MOD    → stable tiebreaker
static inline int tuple_lt(uint16_t ca_g, uint16_t ca_r, uint64_t ha,
                            uint16_t cb_g, uint16_t cb_r, uint64_t hb)
{ if (ca_g != cb_g) return ca_g < cb_g;
  if (ca_r != cb_r) return ca_r > cb_r;   // higher reads count wins
  return ha < hb;
}

// Candidate for uniformly sparse thinning — stores sort-key fields so the
// qsort comparator does not need external array access.
typedef struct {
  int      pos;
  int      nw_val;
  uint16_t g_count;
  uint16_t r_count;
  uint64_t hash;
} Candidate;

// Compare candidates by descending h_nw priority:
//   1. nw      descending (more windows → higher priority)
//   2. g_count ascending  (rarer in genome → better seed)
//   3. r_count descending (more covered by reads → better seed)
//   4. hash    ascending  (stable tiebreaker)
static int cand_cmp(const void *a, const void *b)
{ const Candidate *ca = (const Candidate *)a;
  const Candidate *cb = (const Candidate *)b;
  if (ca->nw_val != cb->nw_val)
    return (cb->nw_val > ca->nw_val) ? 1 : -1;
  if (ca->g_count != cb->g_count)
    return (ca->g_count > cb->g_count) ? 1 : -1;
  if (ca->r_count != cb->r_count)
    return (cb->r_count > ca->r_count) ? 1 : -1;
  if (ca->hash != cb->hash)
    return (ca->hash > cb->hash) ? 1 : -1;
  return 0;
}

// Select seeds using a sliding-window (genome_count, -reads_count, hash) minimizer.
//
// For every window of WSIZE consecutive k-mers the one with the lexicographically
// smallest (prof_genome[i], -prof_reads[i], hVal[i]) tuple is selected.
// The result is the union of all per-window selections.
//
//   seq         – pointer to position (K-1) in the read
//   plen        – number of k-mers (rlen - K + 1)
//   K           – k-mer length
//   prof_genome – FastK uint16 genome k-mer counts, positions 0..plen-1
//   prof_reads  – FastK uint16 reads  k-mer counts, positions 0..plen-1
//   read_id     – 0-based read index
//   out         – binary output file
//
// Returns the number of entries written.
int find_seeds(const char    *seq,
               const int      plen,
               const int      K,
               const uint16  *prof_genome,
               const uint16  *prof_reads,
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
    // Reduce to [0, MOD) so that hash ties are common — enables robust winnowing
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
  int  *nw      = (int  *) calloc((size_t)plen,  sizeof(int));   // # of windows for which position i is minimizer
  if (!dq || !is_seed || !nw) { perror("malloc deque/is_seed/nw"); exit(1); }

  int dq_head = 0, dq_tail = 0;
  for (int i = 0; i < plen; i++)
    { // Robust winnowing: evict back elements whose tuple is STRICTLY GREATER
      // than the current tuple.  On a tie we keep the existing (leftmost) entry,
      // so ties are broken in favour of the element already inside the window.
      if (prof_genome[i] != 0)
        { while (dq_head < dq_tail &&
                 tuple_lt(prof_genome[i],             prof_reads[i],             hvals[i],
                          prof_genome[dq[dq_tail-1]], prof_reads[dq[dq_tail-1]], hvals[dq[dq_tail-1]]))
            dq_tail--;
          dq[dq_tail++] = i;
        }

      // Evict front elements that have slid out of the window.  With g==0 positions
      // skipped the deque can be empty, so guard against dq_head == dq_tail.
      while (dq_head < dq_tail && dq[dq_head] < i - WSIZE + 1)
        dq_head++;

      // Only consider complete windows of exactly WSIZE k-mers (standard minimizer
      // definition).  Emit only if the window holds a g>=1 candidate (deque non-empty);
      // an all-g==0 window produces no seed.
      if (i >= WSIZE - 1 && dq_head < dq_tail)
        { is_seed[dq[dq_head]] = true;
          nw[dq[dq_head]]++;
        }
#ifdef DEBUG_HASH
      if (i >= WSIZE - 1 && dq_head < dq_tail)
      { int win_start = i - WSIZE + 1;
        int min_pos   = dq[dq_head];
        fprintf(stderr, "WINDOW read=%d [%d..%d]  minimizer=@%d g=%u r=%u h=%lu\n",
                read_id, win_start, i, min_pos,
                (unsigned)prof_genome[min_pos], (unsigned)prof_reads[min_pos], hvals[min_pos]);
        for (int j = win_start; j <= i; j++) {
          const char *wkp = seq - K + 1 + j;
          char wkstr[65];
          for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
          wkstr[K] = '\0';
          fprintf(stderr, "  j=%-6d g=%-5u r=%-5u h=%-4lu  %s%s\n",
                  j,
                  (unsigned)prof_genome[j], (unsigned)prof_reads[j], hvals[j],
                  wkstr,
                  (j == min_pos) ? "  <-- MIN" : "");
        }
      }
#endif
    }

  // Step 3: Uniformly sparse thinning.
  // Sort candidates by descending h_nw = (NW, -h_minimax) priority, then
  // greedily accept only those whose position is ≥ RSPACE from every
  // already-accepted seed.
  { int ncand = 0;
    for (int i = 0; i < plen; i++)
      if (is_seed[i]) ncand++;

    if (ncand > 1)
      { Candidate *cands = (Candidate *) malloc((size_t)ncand * sizeof(Candidate));
        if (!cands) { perror("malloc cands"); exit(1); }
        { int ci = 0;
          for (int i = 0; i < plen; i++)
            if (is_seed[i])
              { cands[ci].pos     = i;
                cands[ci].nw_val  = nw[i];
                cands[ci].g_count = prof_genome[i];
                cands[ci].r_count = prof_reads[i];
                cands[ci].hash    = hvals[i];
                ci++;
              }
        }
        qsort(cands, (size_t)ncand, sizeof(Candidate), cand_cmp);

        // Reset is_seed; re-mark only accepted positions.
        memset(is_seed, 0, (size_t)plen * sizeof(bool));

        // Maintain accepted positions in sorted order (by position) for
        // O(log n) nearest-neighbour distance checks.
        int *acc = (int *) malloc((size_t)ncand * sizeof(int));
        if (!acc) { perror("malloc acc"); exit(1); }
        int nacc = 0;

        for (int ci = 0; ci < ncand; ci++)
          { int p = cands[ci].pos;
            // Binary search for insertion point in acc[].
            int lo = 0, hi = nacc;
            while (lo < hi)
              { int mid = (lo + hi) / 2;
                if (acc[mid] < p) lo = mid + 1;
                else              hi = mid;
              }
            // Check left and right neighbours in accepted set.
            int ok = 1;
            if (lo > 0    && p - acc[lo-1] < RSPACE) ok = 0;
            if (lo < nacc && acc[lo] - p   < RSPACE) ok = 0;
            if (ok)
              { memmove(acc + lo + 1, acc + lo,
                        (size_t)(nacc - lo) * sizeof(int));
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
  for (int i = 0; i < plen; i++)
    { const char *wkp = seq - K + 1 + i;
      char wkstr[65];
      for (int b2 = 0; b2 < K; b2++) wkstr[b2] = wkp[b2];
      wkstr[K] = '\0';
      fprintf(stderr, "  i=%-6d nw=%-5d g=%-5u r=%-5u h=%-4lu  %s  %s%s\n",
              i, nw[i],
              (unsigned)prof_genome[i], (unsigned)prof_reads[i], hvals[i],
              wkstr, nw[i] ? "*" : "", is_seed[i] ? "****" : "");
    }
#endif

  // Step 4: write selected seeds.
  int written = 0;
  for (int i = 0; i < plen; i++)
    { if (!is_seed[i]) continue;
      const char *kp = seq - K + 1 + i;
      uint64_t c = 0;
      for (int x = 0; x < K; x++)
        c = (c << 2) | base2bit[(unsigned char)kp[x]];
      c &= KMASK;
      int32_t rpos = (int32_t)i;
      // Every seed's rid carries the genome copy-tier in its top 3 bits (see copynum.h):
      // the low 29 bits stay ref/read id. Strip the tier with CN_VAL wherever it needs the id.
      // read_id < CN_MAXVAL (2^29) is enforced in main(), so we pack unconditionally.
      uint32_t rid = CN_PACK(read_id, cn_tier_of(prof_genome[i]));
      fwrite(&c,    sizeof(uint64_t), 1, out);
      fwrite(&rpos, sizeof(int32_t),  1, out);
      fwrite(&rid,  sizeof(uint32_t), 1, out);
      written++;
    }

  free(is_seed);
  free(nw);
  free(dq);
  free(hvals);
  return written;
}

int main(int argc, char *argv[]) {
  DAZZ_DB _db, *db = &_db;
  // Optional density flags (-w<WSIZE> -s<RSPACE>) precede the 4 positional args.
  { int i, j = 1;
    for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'k': K      = atoi(argv[i]+2); break;
          case 'w': WSIZE  = atoi(argv[i]+2); break;
          case 's': RSPACE = atoi(argv[i]+2); break;
          default:
            fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[i]);
            return 1;
        }
      } else
        argv[j++] = argv[i];
    }
    argc = j;
  }
  if (argc < 5 || WSIZE < 1 || RSPACE < 0) {
    fprintf(stderr, "Usage: %s [-k<32>] [-w<100>] [-s<100>] <db> <ref.prof> <ref.reads.prof> <out.kpos>\n", argv[0]);
    fprintf(stderr, "  -k<N>           : k-mer length (<= 32)\n");
    fprintf(stderr, "  -w<N>           : minimizer window width  (larger => sparser seeds)\n");
    fprintf(stderr, "  -s<N>           : min spacing of uniform sparse seeds (larger => sparser)\n");
    fprintf(stderr, "  <db>            : DAZZ reference dam/db\n");
    fprintf(stderr, "  <ref.prof>      : FastK genome k-mer count profile\n");
    fprintf(stderr, "  <ref.reads.prof>: FastK reads  k-mer count profile\n");
    fprintf(stderr, "  <out.kpos>      : output k-mer position file\n");
    return 1;
  }
  if (K > 32) {
    fprintf(stderr, "%s: K=%d out of range; K must <= 32\n", argv[0], K);
    return 1;
  }
  int Km1 = K - 1;
  KMASK = (K >= 32) ? UINT64_MAX : ((UINT64_C(1) << (2 * K)) - 1);
  fprintf(stderr, "K=%d WSIZE=%d RSPACE=%d\n", K, WSIZE, RSPACE);
  Open_DB(argv[1], db);
  // The copy-tier occupies the rid's top 3 bits (bits 29..31), so the contig/read id must
  // fit in the low 29 bits.  A larger DB would overflow the id into the tier bits and make
  // himapper mis-identify reads — fail loudly rather than emit a corrupt .kpos.
  if (db->nreads > (int) CN_MAXVAL) {
    fprintf(stderr, "%s: %d reads/contigs exceeds the copy-tier rid limit %d (2^29); "
                    "the id would collide with the tier bits\n",
            argv[0], db->nreads, (int) CN_MAXVAL);
    return 1;
  }

  // Open genome profile (argv[2]).
  Profile_Index *Pg = Open_Profiles(argv[2]);
  if (Pg == NULL)
    { fprintf(stderr, "%s: Cannot open genome profile '%s'\n", argv[0], argv[2]);
      return 1;
    }
  if (Pg->kmer != K)
    { fprintf(stderr, "%s: genome profile K=%d does not match K=%d\n", argv[0], Pg->kmer, K);
      Free_Profiles(Pg);
      return 1;
    }
  fprintf(stderr, "Genome profile:  K=%d  nreads=%lld\n",
          Pg->kmer, (long long)Pg->nbase[Pg->nparts-1]);

  // Open reads profile (argv[3]).
  Profile_Index *Pr = Open_Profiles(argv[3]);
  if (Pr == NULL)
    { fprintf(stderr, "%s: Cannot open reads profile '%s'\n", argv[0], argv[3]);
      Free_Profiles(Pg);
      return 1;
    }
  if (Pr->kmer != K)
    { fprintf(stderr, "%s: reads profile K=%d does not match K=%d\n", argv[0], Pr->kmer, K);
      Free_Profiles(Pg); Free_Profiles(Pr);
      return 1;
    }
  fprintf(stderr, "Reads  profile:  K=%d  nreads=%lld\n",
          Pr->kmer, (long long)Pr->nbase[Pr->nparts-1]);

  int     pmax       = 20000;
  uint16 *prof_genome = (uint16 *) malloc((size_t)pmax * sizeof(uint16));
  uint16 *prof_reads  = (uint16 *) malloc((size_t)pmax * sizeof(uint16));
  if (!prof_genome || !prof_reads) { perror("malloc profile"); return 1; }

  FILE *out = fopen(argv[4], "wb");
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
  int32_t entry_size = (int32_t)(sizeof(uint64_t) + sizeof(int32_t) + sizeof(uint32_t));
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

      // Fetch genome k-mer counts for this contig/read (0-based id).
      int fetched = Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
      if (fetched > pmax)
        { pmax       = (int)(1.2 * fetched) + 1000;
          prof_genome = (uint16 *) realloc(prof_genome, (size_t)pmax * sizeof(uint16));
          prof_reads  = (uint16 *) realloc(prof_reads,  (size_t)pmax * sizeof(uint16));
          if (!prof_genome || !prof_reads) { perror("realloc profile"); exit(1); }
          Fetch_Profile(Pg, (int64)id, pmax, prof_genome);
        }

      // Fetch reads k-mer counts for the same position range.
      Fetch_Profile(Pr, (int64)id, pmax, prof_reads);

      int written = find_seeds(seq + Km1, plen, K, prof_genome, prof_reads, id, out);
      total_written += (uint64_t)written;
    }

    // Seek back and write the total entry count into header (offset 16)
    if (fseek(out, 16, SEEK_SET) == 0) {
      fwrite(&total_written, sizeof(uint64_t), 1, out);
    }
    fprintf(stderr, "Seeds written: %llu\n", (unsigned long long)total_written);

    free(seq-1);
  }

  free(prof_genome);
  free(prof_reads);
  Free_Profiles(Pg);
  Free_Profiles(Pr);

  fclose(out);

  return 0;
}
