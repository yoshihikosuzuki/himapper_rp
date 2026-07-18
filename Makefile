CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = modimizer himapper HPC.himapper copy_minimizer minimaxmer dump_kpos

all: $(ALL)

modimizer: modimizer.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o modimizer modimizer.c DB.c QV.c -lm

copy_minimizer: copy_minimizer.c DB.c DB.h QV.c QV.h libfastk.c libfastk.h gene_core.h
	gcc $(CFLAGS) -DSKIP_GENE_CORE -o copy_minimizer copy_minimizer.c DB.c QV.c libfastk.c -lpthread -lm

minimaxmer: minimaxmer.c DB.c DB.h QV.c QV.h libfastk.c libfastk.h gene_core.h copynum.h
	gcc $(CFLAGS) -DSKIP_GENE_CORE -o minimaxmer minimaxmer.c DB.c QV.c libfastk.c -lpthread -lm

himapper: himapper.c map.c map.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o himapper himapper.c map.c align.c DB.c QV.c -lpthread -lm

HPC.himapper: HPC.himapper.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.himapper HPC.himapper.c DB.c QV.c -lm

dump_kpos: dump_kpos.c copynum.h
	gcc $(CFLAGS) -o dump_kpos dump_kpos.c

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile
