CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = modimizer himapper HPC.himapper

all: $(ALL)

modimizer: seed_mod.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o modimizer seed_mod.c DB.c QV.c -lm

himapper: himapper.c map.c map.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o himapper himapper.c map.c align.c DB.c QV.c -lpthread -lm

HPC.himapper: HPC.himapper.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.himapper HPC.himapper.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile
