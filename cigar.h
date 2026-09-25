#ifndef HIMAPPER_CIGAR_H
#define HIMAPPER_CIGAR_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

/* BAM operation codes, with lengths stored in the high 28 bits. */

enum
  { CG_I  = 1,
    CG_D  = 2,
    CG_EQ = 7,
    CG_X  = 8
  };

typedef struct
  { uint32_t *ops;
    size_t    n, cap;
  } Cigar;

void Cigar_Push(Cigar *c, int op, int len);
void Cigar_Append(Cigar *dst, const Cigar *src);
void Cigar_Reverse(Cigar *c);
void Cigar_Free(Cigar *c);

/* Split at an exact (query,reference) coordinate; failure leaves dst empty. */

int Cigar_Slice(Cigar *dst, const Cigar *src, int q, int r, int suffix);

void Cigar_Count(const Cigar *c, int64_t *q, int64_t *r, int64_t *matches,
                 int64_t *columns, int64_t *edits);
int Cigar_Print(FILE *out, const Cigar *c, int reverse);

#endif
