#!/bin/bash

ml fastk dazz_db

mkdir -p tmp

K=32

fasta2DAM ref ref.fasta
fasta2DAM reads reads.fasta

FastK -k${K} -v -t1 -p -Ptmp ref.fasta
FastK -k${K} -v -t1 -p -Ptmp reads.fasta

FastK -k${K} -v -t1 -p:ref -Nreads.ref -Ptmp reads.dam
FastK -k${K} -v -t1 -p:reads -Nref.reads -Ptmp ref.dam
