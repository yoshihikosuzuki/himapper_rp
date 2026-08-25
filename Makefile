CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
DAZZ_DB_DIR = DAZZ_DB
CPPFLAGS += -I$(DAZZ_DB_DIR)

DAZZ_DB_SRC = $(DAZZ_DB_DIR)/DB.c $(DAZZ_DB_DIR)/QV.c
DAZZ_DB_HDR = $(DAZZ_DB_DIR)/DB.h $(DAZZ_DB_DIR)/QV.h

ALL = himapper minimaxmer dump_kpos

all: $(ALL)

minimaxmer: minimaxmer.c minimaxmer.h nthash.h $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) libfastk.c libfastk.h gene_core.h copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSKIP_GENE_CORE -o minimaxmer minimaxmer.c $(DAZZ_DB_SRC) libfastk.c -lpthread -lm

minimaxmer-lib.o: minimaxmer.c minimaxmer.h nthash.h $(DAZZ_DB_HDR) libfastk.h gene_core.h copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	  -DMINIMAXMER_NO_MAIN -DSKIP_GENE_CORE -c -o $@ minimaxmer.c

himapper: himapper.c map.c map.h align.c align.h minimaxmer.h minimaxmer-lib.o \
          $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) libfastk.c libfastk.h copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSKIP_GENE_CORE -o himapper himapper.c map.c align.c \
	  minimaxmer-lib.o $(DAZZ_DB_SRC) libfastk.c -lpthread -lm

dump_kpos: dump_kpos.c copynum.h
	$(CC) $(CFLAGS) -o dump_kpos dump_kpos.c

clean:
	rm -f $(ALL) minimaxmer-lib.o
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile
