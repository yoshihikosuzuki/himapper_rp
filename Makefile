CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
DAZZ_DB_DIR = DAZZ_DB
FASTK_DIR = FASTK
CPPFLAGS += -I$(DAZZ_DB_DIR) -I$(FASTK_DIR)

DAZZ_DB_SRC = $(DAZZ_DB_DIR)/DB.c $(DAZZ_DB_DIR)/QV.c
DAZZ_DB_HDR = $(DAZZ_DB_DIR)/DB.h $(DAZZ_DB_DIR)/QV.h
FASTK_SRC = $(FASTK_DIR)/libfastk.c
FASTK_HDR = $(FASTK_DIR)/libfastk.h $(FASTK_DIR)/gene_core.h

ALL = himapper minimaxmer

all: $(ALL)

local-tools:
	$(MAKE) -C $(DAZZ_DB_DIR)
	$(MAKE) -C $(FASTK_DIR)
	$(MAKE) -C fastkdual FastKdual

minimaxmer: minimaxmer.c minimaxmer.h nthash.h $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) $(FASTK_SRC) $(FASTK_HDR) copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSKIP_GENE_CORE -o minimaxmer minimaxmer.c $(DAZZ_DB_SRC) $(FASTK_SRC) -lpthread -lm

minimaxmer-lib.o: minimaxmer.c minimaxmer.h nthash.h $(DAZZ_DB_HDR) $(FASTK_HDR) copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	  -DMINIMAXMER_NO_MAIN -DSKIP_GENE_CORE -c -o $@ minimaxmer.c

himapper: himapper.c map.c map.h mapq.c mapq.h align.c align.h cigar.c cigar.h paf.c paf.h minimaxmer.h minimaxmer-lib.o \
	          $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) $(FASTK_SRC) $(FASTK_HDR) copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSKIP_GENE_CORE -o himapper himapper.c map.c mapq.c align.c cigar.c paf.c \
	  minimaxmer-lib.o $(DAZZ_DB_SRC) $(FASTK_SRC) -lpthread -lm

clean:
	rm -f $(ALL) minimaxmer-lib.o
	rm -fr *.dSYM
	rm -f himapper.tar.gz

package:
	make clean
	tar -zcf himapper.tar.gz README.md LICENSE *.h *.c Makefile

.PHONY: all fastkdual local-tools clean package
