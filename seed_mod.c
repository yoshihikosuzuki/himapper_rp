#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "DB.h"
#include "nthash.h"

#define DEBUG_HASH

const int K = 40;
const int Km1 = K - 1;
const int MOD = 1001;
const int MOD_THRES = 1;

static void kmer_hash(const char *seq, int *hash, int plen, int K) {
  // canonical, forward, and reverse-strand hash values
  uint64_t hVal, fhVal, rhVal;

  fhVal = rhVal = 0;
  hVal = NTC64_b(seq - K + 1, K, &fhVal, &rhVal);  // initial hash value
  hash[0] = (int)(hVal % MOD);

#ifdef DEBUG_HASH
  fprintf(stderr, "Canonical hash values:\n");
  if (hash[0] < MOD_THRES)
    fprintf(stderr,
            "  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n", 0,
            K, seq - K + 1, hash[0], hVal, fhVal, rhVal);
#endif

  for (int i = 1; i < plen; i++) {
    hVal = NTC64_c(seq[i - K], seq[i], K, &fhVal,
                   &rhVal);  // consecutive hash values
    hash[i] = (int)(hVal % MOD);

#ifdef DEBUG_HASH
    if (hash[i] < MOD_THRES)
      fprintf(stderr,
              "  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n",
              i, K, seq - K + 1 + i, hash[i], hVal, fhVal, rhVal);
#endif
  }

  return;
}

void find_seeds(const char *seq, const char *class, int *hash, const int plen,
                const int K) {
  // Compute canonical hash for every k-mer
  kmer_hash(seq + K - 1, hash, plen, K);

  return;
}

int main(int argc, char *argv[]) {
  DAZZ_DB _db, *db = &_db;
  Open_DB(argv[1], db);

  {
    int rlen_max = db->maxlen;
    char *seq = New_Read_Buffer(db);
    int *hash = Malloc((rlen_max + 1) * sizeof(int), "Hash array");

    for (int id = 0; id < db->nreads; id++) {
      DAZZ_READ *r = db->reads + id;
      int rlen = r->rlen;
      int plen = rlen-Km1;

      Load_Read(db, id, seq, 2);
      printf("%s\n", seq);
      kmer_hash(seq+Km1,hash,plen,K);

      // for (int i = 0; i < plen; i ++) {
      //   printf("")
      // }
    }
  }

  return 0;
}
