FASTK_DIR ?= ../FASTK

CC       ?= gcc
CPPFLAGS ?=
CFLAGS   ?= -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDFLAGS  ?=

FASTK_SOURCES = $(FASTK_DIR)/io.c \
                $(FASTK_DIR)/split.c \
                $(FASTK_DIR)/table.c \
                $(FASTK_DIR)/merge.c \
                $(FASTK_DIR)/MSDsort.c \
                $(FASTK_DIR)/LSDsort.c \
                $(FASTK_DIR)/libfastk.c

FASTK_HEADERS = $(FASTK_DIR)/FastK.h \
                $(FASTK_DIR)/libfastk.h \
                $(FASTK_DIR)/gene_core.h \
                $(FASTK_DIR)/gene_core.c

FASTK_LIBS = $(FASTK_DIR)/LIBDEFLATE/libdeflate.a \
             $(FASTK_DIR)/HTSLIB/libhts.a

HTS_LIBS ?= -lz -lm -lbz2 -llzma -lcurl

.PHONY: all clean

all: FastKdual compare_profiles

FastKdual: FastKdual.c FastKdual.h count_dual.c $(FASTK_DIR)/count.c \
           $(FASTK_SOURCES) $(FASTK_HEADERS) $(FASTK_LIBS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	  -I$(FASTK_DIR) -I$(FASTK_DIR)/HTSLIB \
	  -o $@ FastKdual.c count_dual.c $(FASTK_SOURCES) \
	  $(FASTK_LIBS) -lpthread $(HTS_LIBS)

compare_profiles: compare_profiles.c $(FASTK_DIR)/libfastk.c \
                  $(FASTK_DIR)/libfastk.h $(FASTK_DIR)/gene_core.h \
                  $(FASTK_DIR)/gene_core.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) \
	  -I$(FASTK_DIR) -o $@ compare_profiles.c \
	  $(FASTK_DIR)/libfastk.c -lpthread -lm

clean:
	rm -f FastKdual compare_profiles
