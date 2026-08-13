# himapper_rp

Fast, haplotype-aware HiFi read mapper using relative k-mer profiles.

## Local DAZZ_DB and FASTK

himapper is built against the `DAZZ_DB/` clone in this repository, and the
preparation workflows call `DAZZ_DB/fasta2DAM`, `DAZZ_DB/DBshow`,
`DAZZ_DB/DBsplit`, and `FASTK/FastK` by absolute path.  Build the local
dependencies and himapper with:

```sh
make local-tools
make
```

On compute nodes without the system curl development files, load the dependency
documented by `~/code/install-scripts/pg/app/fastk.sh` while building:

```sh
source /bio/lmod/lmod/init/bash
module purge
module use /bio/package/.modulefiles
module load curl/7.65.3
make local-tools
module purge
make
```

The local DAZZ_DB parser accepts block sizes wider than the historical
11-character field while retaining the fixed-width stub layout when writing.

## Current workflow

### One-pass reads profiles

`fastkdual/FastKdual` is a separate FastK-derived executable.  It leaves the
vendored `FASTK/FastK` source and behavior unchanged, but reuses one input
partition and weighted k-mer sort to produce both the reads self
table/profile and the reads-vs-reference profile:

```sh
fastkdual/FastKdual -k32 -T16 -t1 \
  -p:ref -Nreads -Rreads.ref -P"$TMPDIR/fastk" reads.dam
```

The command creates the same `reads.hist`, `reads.ktab`, `reads.prof`, and
`reads.ref.prof` artifacts that previously required separate self and relative
FastK runs.  Build it together with the local dependencies using
`make local-tools`.

### Seed selection and mapping

`minimaxmer` takes a whole, partitioned DAZZLER DB and creates one KPOS
container plus its block-offset sidecar. DB blocks are processed in parallel
with `-T`, while entries remain block-local and code-sorted within each block.

```sh
minimaxmer -T16 -l10 ref.dam ref.prof ref.reads.prof .ref.minimax.kpos
minimaxmer -T16 reads.dam reads.ref.prof reads.prof .reads.minimax.kpos
```

These commands also create `.ref.minimax.kpos.idx` and
`.reads.minimax.kpos.idx`.

Long homopolymer, dinucleotide-satellite, and trinucleotide-satellite tracts
are excluded from the reference command above.  The option is off by default
(`-l-1`), as used by the read command.  `minimaxmer` does not infer the DB
role; the caller decides whether to pass `-l`.  A positive value `N` removes
selected reference K-mers that overlap any exact period-1, -2, or -3 tract
whose complete-motif span is strictly greater than `N` bases.  Thus the
recommended chr6 cutoff of 10 first excludes 11-base homopolymers, 6 complete
dinucleotide units (12 bases), and 4 complete trinucleotide units (12 bases).
This cutoff improved the measured chr6 result at `w=s=100`; its whole-genome
generality has not yet been tested.  Overlap means any of the K bases, not only
the stored K-mer start.
Filtering occurs after sparse thinning, so enabled output is an exact subset
of the corresponding default output and no replacement minimizer is added.
The `PERIODIC_MASK_V1` counters named `*_mask_bases` report the full maximal
tract span actually masked, including a possible terminal partial motif unit;
the cutoff comparison itself always uses complete motif units.

The default seed selector keeps the historical
`(genome_count, -reads_count, canonical_ntHash % 2003)` priority.  When a
window has several occurrences of the same minimum tuple, it uses robust
winnowing: the selected occurrence is retained while it remains a minimum.
The traversal direction of each tied subrun is chosen with
`D(W) = #G + #T - #A - #C` over the `K + W - 1` spelled bases, so the selected
positions mirror between a sequence and its reverse complement.  `D(W)` is
computed lazily only for tied windows.

Map the whole reads DB in one process. himapper loads the reference KPOS once,
then processes every reads DB block internally.  A fixed worker team is reused
throughout each reads block instead of creating a new pthread gang for every
sort, match, and report phase:

```sh
himapper -1 -T16 \
  ref.dam .ref.minimax.kpos \
  reads.dam .reads.minimax.kpos
```

Read seed generation can instead be fused with mapping.  Serial mode gives
seed generation and mapping the full `-T` thread budget in turn and passes each
sorted KPOS block directly in memory.  It is the no-overlap control for the
fused path.  Compared with separate standalone seed and map processes, its
difference includes removal of intermediate KPOS write/read, duplicate read-DB
loading, process-boundary overhead, and a small scheduling difference; it must
not be attributed to file I/O alone:

```sh
himapper -1 -T16 \
  --read-seed-profiles reads.ref.prof reads.prof \
  --seed-mode serial --seed-threads 16 \
  --seed-window 100 --seed-spacing 100 \
  ref.dam .ref.minimax.kpos reads.dam generate
```

Fused/generated read seeds are not periodic-filtered.  To use low-complexity
masking, create the reference KPOS with standalone `minimaxmer -l10` as shown
above and pass it with ordinary unfiltered read KPOS or generated read seeds.

Pipeline mode overlaps generation of the next block with mapping of the
current block while keeping the combined budget at `-T` threads.  If `S` is
given by `--seed-threads`, seed generation uses `S` workers and mapping uses
`T-S`.  The measured 16-thread default is 6 seed workers plus 10 map workers:

```sh
himapper -1 -T16 \
  --read-seed-profiles reads.ref.prof reads.prof \
  --seed-mode pipeline --seed-threads 6 \
  --seed-window 100 --seed-spacing 100 \
  ref.dam .ref.minimax.kpos reads.dam generate
```

The generated read KPOS is not written by default.  Add
`--save-read-kpos .reads.minimax.kpos` only when the same seeds must be reused;
omit it from timing runs that are intended to measure the no-intermediate-file
path.

The fused path functionally supports serial budgets from `T=1` through 128
and pipeline budgets from `T=2` through 128.  Serial always uses `T/T` in
sequence; pipeline accepts any split satisfying `1 <= S < T`.  Only the
16-thread `S=6`, `T-S=10` pipeline has been measured on the WGS benchmark, so
performance at other totals or splits is unverified.  Exact map worker counts
are selected automatically for the fused path.  Standalone mapping from an
existing read KPOS retains the historical power-of-two thread rounding.
Serial seed generation keeps up to `T` generated DB/KPOS blocks resident so
its RSS grows with the thread budget; the T64/T128 WGS memory footprint has not
been measured.  Pipeline keeps up to `S` generated blocks in flight and is the
safer default at large thread counts.

The WGS stage-07 script defaults to pipeline mode and chooses the nearest
integer to `3*T/8` seed workers (6 at `T=16`), clamped to `1..T-1`.  Override
it with `HIMAPPER_SEED_MODE=serial|pipeline` and, for pipeline mode,
`HIMAPPER_SEED_THREADS=S`.  Serial rejects a conflicting seed-thread override
because it uses `T/T`.  Set `HIMAPPER_SAVE_READ_KPOS` to a nonempty output path
only when the generated read KPOS must be retained; the default timed run
writes timing to `seed-map.time` and does not save it.

A multi-block reads DB produces `reads.1.ref.las`,
`reads.2.ref.las`, and so on. A single-block DB produces `reads.ref.las`.

When converting LAS output back to FASTA scaffold coordinates, first retain
the N-delimited DAM layout and pass it to the placement converter:

```sh
DBshow -n ref.dam > ref.dam.layout
LAshow reads ref.dam reads.ref.las |
  work_real/scripts/lashow_to_placement.py \
    --ref-layout ref.dam.layout reads.names ref.names > mapping.tsv
```

The layout is required because a FASTA scaffold containing N runs becomes
multiple DAM contigs.  Forward LAshow coordinates already contain `fpulse`;
the converter adjusts only reverse hits for himapper's per-contig reverse
complement convention.  By default, the converter reconstructs DALIGNER
START/NEXT chains and emits the maximum-query-span chain for each read.
Pass `--all` only for secondary diagnostics; it emits every reconstructed
candidate chain in input order, one chain per row.  Because this produces
multiple rows per read, do not use it as input to primary evaluation or
coverage plotting.





---

Below are the original description of DAMAPPER

# The Damapper Library

## _Author:  Gene Myers_
## _First:   July 11, 2016_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/damapper-commands).

The commands below permit one to find the best, or several good alternative,
locations that a read may have been sampled from in a reference genome.  The
code is designed for PACBIO reads with an assumption of a 15% average error
rate and occational large drop out regions where the quality is exceedingly
low.

```
1. damapper [-vpzCN] [-k<int(20)>] [-t<int>] [-M<int>] [-T<int(4)>]
                     [-e<double(.85)] [-s<int(100)>] [-n<double(1.)>]
                     [-m<track>]+ <ref:db|dam> <reads:db|dam> ...
```

Search the reference data base \<ref\> for the best matches to each read in the list of
databases or database blocks \<reads\>.  The parameters -v, -k, -t, -M, -T -e, -s,
and -m are exactly as for the daligner (see here).  The parameter -z asks that LAs are
sorted in pile order as opposed to map order (see the -a option of daligner for which this
is the negation).  Matches at the expected correlation
given by -e are sought and unlike daligner they can be of any length (i.e. short).
Alignments at much lower rates are not forced in order to report a single matching
interval, but rather when a read has bad stretches, the result is a chain of local
alignments over the good segments of the read.  Damapper reports the best chain
covering a segment of a read, and may report several chains covering disjoint segments
(e.g. the case of a chimeric read).  If the -n option is given then all chains that
are within given fraction of optimal are also reported, e.g. -n.95 reports all matches
within 95% of the top match.

For a given reference data base X and read block Y, damapper produces the single file
Y.X..las.
Note the order of X and Y is reversed
from that on the command line (which is different than the conventions for
daligner).
The reversal of the order of the blocks emphasizes that in the result, the A-reads
are the mapped reads, and the B-reads are contigs of the reference.  Each file is
sorted in order of the A-reads, and if a match is a chain of local alignments the
la's in the chain occur in increasing order of A-coordinate.  Moreover, if several
near optimal matches are reported, matches appear in order of score starting
with the best.

Itf the -C option is set, then damapper also outputs a file X.Y.las for a given
block pair that contains all
the same matches as in X.Y.las but where the A-read is a contig of the
reference and the B-read is a mapped read.  And if the -N options is set, then the
files Y.X.las is not produced.

The -p option requests that damapper produce a repeat profile track for each read.
For each trace-point sized interval (established by the -s parameter) the number of
times, c, this segment is involved in a distinct alignment to the reference is
estimated.  0 is recorded for segments that don't match anything in the reference,
1 for segments that match uniquely, and otherwise floor(log<sub>10</sub>c/10) up to a cap of
40 corresponding to 10,000 copies.  Observe carefully that the track has the same
form as intrinsic quality values.  They can be output by DBdump and DaViewer is able
to graphically display them.  These profiles should make it obvious when a read does
not have a unique location in a reference sequence due to its being entirely or almost
entirely repetitive.

```
2. HPC.damapper [-vbpzCN]
                [-k<int(20)>] [-t<int>] [-M<int>] [-e<double(.85)] [-s<int(100)]
                [-n<double(1.)>] [-m<track>]+ [-B<int(4)>] [-T<int(4)>] [-f<name>]
                <ref:db|dam> <reads:db|dam> [<first:int>[-<last:int>]]
```

HPC.damapper writes a UNIX shell script to the standard output that maps every read in
blocks \<first\> to \<last\> of database \<reads\> to a reference sequence \<ref\>.  If \<last\>
is missing then only the single block \<first\> is mapped, and if \<first\> is also missing
then all blocks of the database are mapped.  Except for the -B and -f options, all
other options are passed through to damapper, save of -v which is passed to every
program in the script.  The -B option determines the maximum number of blocks mapped
per call to damapper in the script.

The command script output by HPC.damapper and other HPC.\<x\> programs consists of
command blocks each of which begins with a comment line (begins with #) followed by a
potentially long list of lines each containing a shell command.  Command blocks whose
comment mentions "jobs" and gives the number of said in parenthesis, we call parallel
blocks because each command line in the block can be sent to a node in a cluster for
independent execution, i.e. none of the commands in a block depend on another in the
block.  The remaining command blocks we call house-keeping blocks because they can be
executed by the shell on the launch/server node and the commands are either checking
the integrity of .las files with LAcheck, or removing intermediate files with rm. Each
block should be performed in the order given and should complete before the next
block is performed.

If the -f option is set, then each command block is written to a file with a name of
the form \<name\>.#.\<description\> where \<name\> is specified by the user in the -f option
argument, # gives the order in which the command block in the given file is to be
performed in relation to other command block files, and \<description\> is a (very)
short symbolic reminder of what the block is doing.  For example,
"HPC.damapper -fJOBS REF DB" would produce the files:

```
    JOBS.01.OVL
    JOBS.02.CHECK.OPT
```

The files with the suffix .OPT are optional and need not be executed albeit we highly
recommend that one run the CHECK block.
