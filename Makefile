CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = himapper minimaxmer dump_kpos

all: $(ALL)

minimaxmer: minimaxmer.c DB.c DB.h QV.c QV.h libfastk.c libfastk.h gene_core.h copynum.h
	gcc $(CFLAGS) -DSKIP_GENE_CORE -o minimaxmer minimaxmer.c DB.c QV.c libfastk.c -lpthread -lm

himapper: himapper.c map.c map.h align.c align.h DB.c DB.h QV.c QV.h copynum.h
	gcc $(CFLAGS) -o himapper himapper.c map.c align.c DB.c QV.c -lpthread -lm

dump_kpos: dump_kpos.c copynum.h
	gcc $(CFLAGS) -o dump_kpos dump_kpos.c

clean:
	rm -f $(ALL) HPC.himapper modimizer copy_minimizer
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile

.PHONY: all clean package
