CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

DAZZ_DB_DIR ?= DAZZ_DB
FASTK_DIR ?= FASTK
CPPFLAGS += -I$(DAZZ_DB_DIR)

DAZZ_DB_SRC = $(DAZZ_DB_DIR)/DB.c $(DAZZ_DB_DIR)/QV.c
DAZZ_DB_HDR = $(DAZZ_DB_DIR)/DB.h $(DAZZ_DB_DIR)/QV.h

ALL = himapper minimaxmer dump_kpos

all: $(ALL)

minimaxmer: minimaxmer.c $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) libfastk.c libfastk.h gene_core.h copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSKIP_GENE_CORE -o minimaxmer minimaxmer.c $(DAZZ_DB_SRC) libfastk.c -lpthread -lm

himapper: himapper.c map.c map.h align.c align.h $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -o himapper himapper.c map.c align.c $(DAZZ_DB_SRC) -lpthread -lm

dump_kpos: dump_kpos.c copynum.h
	$(CC) $(CFLAGS) -o dump_kpos dump_kpos.c

local-tools:
	$(MAKE) -C $(DAZZ_DB_DIR) fasta2DB fasta2DAM DBshow DBsplit DBstats DBtrim
	$(MAKE) -C $(FASTK_DIR) deflate.lib
	$(MAKE) -C $(FASTK_DIR) libhts.a
	$(MAKE) -C $(FASTK_DIR) FastK Fastrm

clean:
	rm -f $(ALL) HPC.himapper modimizer copy_minimizer
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile

.PHONY: all local-tools clean package
