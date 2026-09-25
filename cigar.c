#include "cigar.h"

#include <stdlib.h>
#include <limits.h>

void Cigar_Push(Cigar *c, int op, int len)
{ const unsigned maxrun = 0x0fffffff;

  if (len <= 0)
    return;
  while (len)
    { unsigned part = (unsigned) len > maxrun ? maxrun : (unsigned) len;

      if (c->n && (c->ops[c->n-1] & 15) == (unsigned) op)
        { unsigned old = c->ops[c->n-1] >> 4;
          unsigned add = part < maxrun-old ? part : maxrun-old;

          c->ops[c->n-1] += add << 4;
          len  -= add;
          part -= add;
          if (!part)
            continue;
        }
      if (c->n == c->cap)
        { size_t    cap = c->cap ? c->cap + c->cap/2 + 16 : 64;
          uint32_t *p = realloc(c->ops,cap*sizeof(*p));

          if (!p)
            { perror("CIGAR allocation");
              exit(1);
            }
          c->ops = p;
          c->cap = cap;
        }
      c->ops[c->n++] = (part << 4) | op;
      len -= part;
    }
}

void Cigar_Append(Cigar *dst, const Cigar *src)
{ for (size_t i = 0; i < src->n; ++i)
    Cigar_Push(dst,src->ops[i] & 15,src->ops[i] >> 4);
}

void Cigar_Reverse(Cigar *c)
{ for (size_t i = 0, j = c->n; i < j && i < --j; ++i)
    { uint32_t v = c->ops[i];

      c->ops[i] = c->ops[j];
      c->ops[j] = v;
    }
}

void Cigar_Free(Cigar *c)
{ free(c->ops);
  *c = (Cigar) {0};
}

int Cigar_Slice(Cigar *dst, const Cigar *src, int q, int r, int suffix)
{ int x = 0, y = 0, found = 0;

  dst->n = 0;
  for (size_t i = 0; i < src->n; ++i)
    { int op = src->ops[i] & 15, n = src->ops[i] >> 4;
      int dx = op != CG_D, dy = op != CG_I;

      if (!found)
        { int k = dx ? q-x : r-y;

          if (k >= 0 && k <= n && x+dx*k == q && y+dy*k == r)
            { Cigar_Push(dst,op,suffix ? n-k : k);
              found = 1;
            }
          else if (!suffix)
            Cigar_Push(dst,op,n);
        }
      else if (suffix)
        Cigar_Push(dst,op,n);
      x += dx*n;
      y += dy*n;
    }
  if (!found && x == q && y == r)
    found = 1;
  if (!found)
    dst->n = 0;
  return found;
}

void Cigar_Count(const Cigar *c, int64_t *q, int64_t *r, int64_t *matches,
                 int64_t *columns, int64_t *edits)
{ *q = *r = *matches = *columns = *edits = 0;

  for (size_t i = 0; i < c->n; ++i)
    { int op = c->ops[i] & 15, n = c->ops[i] >> 4;

      if (op != CG_D)
        *q += n;
      if (op != CG_I)
        *r += n;
      *columns += n;
      if (op == CG_EQ)
        *matches += n;
      else
        *edits += n;
    }
}

int Cigar_Print(FILE *out, const Cigar *c, int reverse)
{ static const char ops[] = "MIDNSHP=X";

  for (size_t i = 0; i < c->n; ++i)
    { uint32_t v = c->ops[reverse ? c->n-1-i : i];

      if (fprintf(out,"%u%c",v >> 4,ops[v & 15]) < 0)
        return -1;
    }
  return 0;
}
