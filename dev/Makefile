CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
MINIMAX_SELECTOR_FLAGS ?= -DMINIMAX_HASH_MOD=2003 \
                          -DMINIMAX_HASH_REFINE=0 \
                          -DMINIMAX_ROBUST_WINNOW=1 \
                          -DMINIMAX_CANONICAL_WINNOW=1

ROOT := ..
FASTK_DIR := $(ROOT)/FASTK
CPPFLAGS += -I$(ROOT)/DAZZ_DB -I$(FASTK_DIR) -I$(ROOT)

DAZZ_DB_SRC = $(ROOT)/DAZZ_DB/DB.c $(ROOT)/DAZZ_DB/QV.c
DAZZ_DB_HDR = $(ROOT)/DAZZ_DB/DB.h $(ROOT)/DAZZ_DB/QV.h
FASTK_SRC = $(FASTK_DIR)/libfastk.c
FASTK_HDR = $(FASTK_DIR)/libfastk.h $(FASTK_DIR)/gene_core.h
COMMON_SRC = $(DAZZ_DB_SRC) $(FASTK_SRC)
COMMON_DEP = $(DAZZ_DB_SRC) $(DAZZ_DB_HDR) \
             $(ROOT)/nthash.h $(FASTK_SRC) $(FASTK_HDR) \
             $(ROOT)/copynum.h $(ROOT)/minimaxmer.h

PROGRAMS = copy_minimizer modimizer dump_kpos

all: $(PROGRAMS)

copy_minimizer: copy_minimizer.c $(COMMON_DEP)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(MINIMAX_SELECTOR_FLAGS) \
	  -DSKIP_GENE_CORE -o $@ copy_minimizer.c $(COMMON_SRC) -lpthread -lm

modimizer: modimizer.c $(COMMON_DEP)
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	  -DSKIP_GENE_CORE -o $@ modimizer.c $(COMMON_SRC) -lpthread -lm

dump_kpos: dump_kpos.c $(ROOT)/copynum.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ dump_kpos.c

clean:
	rm -f $(PROGRAMS)

.PHONY: all clean
