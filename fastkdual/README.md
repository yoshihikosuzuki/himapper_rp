# FastKdual

`FastKdual` produces these two FastK products from one input partition and one
super-mer/weighted-k-mer sort:

- a self histogram, k-mer table, and profile;
- a profile relative to an existing FastK table.

The original `FASTK/FastK` source and executable are not modified. `FastKdual`
is separate because its phase 2 keeps the sorted weighted k-mers resident long
enough to build one dual-count inverse list. Each inverse record stores one
position index shared by the self and relative counts. The two profile channels
are then compressed and written sequentially from that immutable list.

## Build

The vendored FastK tree must already have its static HTSlib and libdeflate
libraries built:

```sh
make -C FASTK libhts.a deflate.lib
make -C fastkdual
```

To use a different FastK source tree from the same repository revision, or a
verified source-compatible revision:

```sh
make -C fastkdual FASTK_DIR=/path/to/FASTK
```

`count_dual.c` deliberately reuses private `count.c` helpers, so arbitrary
upstream FastK revisions are not ABI/source-compatible build targets.

## Usage

```text
FastKdual [-k<int(40)>] [-t<int(1)>] -p:<table>[.ktab]
          -R<relative_output_root>
          [-c] [-bc<int>] [-v] [-N<self_output_root>]
          [-P<dir($TMPDIR)>] [-M<int(12)>] [-T<int(4)>]
          <source> ...
```

Both `-p:<table>` and `-R<relative_output_root>` are required. The self and
relative output basenames must differ because their temporary profile fragments
share `-P`.

The vendored FastK input/table/profile modules contain legacy shell calls.
Therefore FastKdual validates `-P`, `-N`, `-R`, `-p`, every positional source,
and the current working directory before `Partition_Input`. Paths may contain
ASCII letters and digits plus:

```text
/ . _ - + : = , @ %
```

Whitespace, quotes, glob characters, command separators/substitutions, and
other shell-sensitive characters are rejected explicitly. This restriction
prevents a user-controlled path from reaching a legacy shell command.

On an error after processing starts, FastKdual performs best-effort cleanup
without a shell or glob. It unlinks only the known self/relative output names
and the enumerated `T`, `L`, `U`, `NS`, and `P` temporary files for the two
selected roots. An output kind is eligible for removal only after this run has
opened it or entered the merge stage that replaces it, so an earlier validation
or resource-limit failure does not delete an existing completed output. It
never recursively removes a directory.

Replacing an existing same-name output is not transactional: after its merge
stage starts, a later allocation, input-fragment, or write failure may remove
the previous output along with the partial replacement. Use a fresh output root
or preserve the old artifacts separately when they must survive a failed rerun.

For the HiMapper reads profiles:

```sh
fastkdual/FastKdual \
  -k32 -T16 -t1 \
  -p:ref \
  -Nreads \
  -Rreads.ref \
  -P"$TMPDIR" \
  reads.dam
```

This produces:

```text
reads.hist
reads.ktab and .reads.ktab.*
reads.prof and .reads.{pidx,prof}.*
reads.ref.prof and .reads.ref.{pidx,prof}.*
```

With `-v`, phase 2 reports cumulative CPU and wall time separately for the
shared sort, self channel, and relative channel.

## Phase-2 memory layout

Let `B=KMAX_BYTES`, `K=KMER_BYTES`, and let `S` be the number of weighted
k-mer records in the largest partition. The relevant record widths are:

```text
weighted k-mer       W = K + 2 + B bytes
ordinary inverse     C = B + 1 bytes
dual-count inverse   D = B + 3 bytes
```

Keeping two ordinary inverse arrays would require `(W+2C)S` bytes until profile
fragment construction. FastKdual instead requires `(W+D)S` while joining.
After the dual inverse has been radix-sorted, both backing allocations are
shrunk: one to `DS` for the immutable dual records and one to `CS` for the
reused profile workspace. During each profile-link sort the corresponding
terms are therefore `(D+C)S`, plus the unchanged sorted-super-mer and
`2*PROF_BYTES*nmers` link arrays.

A one-channel-at-a-time implementation without the dual record does not achieve
this bound with the existing FastK radix sort: sorting the first inverse needs a
second `CS` target, while the weighted k-mer list must remain intact for the
relative join. It consequently reaches `W+2C` again unless a new in-place sort
or a disk spill is introduced.

The `W+2C` comparison above is against the initial two-inverse FastKdual
prototype, not against two sequential invocations of ordinary FastK. Ordinary
FastK only retains one channel's inverse state at a time. FastKdual targets
lower total CPU, wall time, repeated input I/O, and temporary sorting work; its
peak RSS can remain higher than either individual FastK invocation.

In independent single-run work_sim checks, peak RSS was approximately:

```text
T16, one sort partition: ordinary FastK max 865,436 KB; FastKdual 1,026,004 KB (+18.6%)
T16, -M1, 3 sort parts: ordinary FastK max 435,648 KB; FastKdual   436,716 KB (+0.25%)
```

For the same runs, replacing the two sequential FastK invocations reduced
wall/CPU from 3.98/39.72 seconds to 2.05/22.92 seconds in the one-partition
case, and from 2.99/27.91 seconds to 1.56/18.48 seconds with `-M1`.
These measurements are workload-, allocator-, and cache-dependent. Whole-genome
FastKdual peak RSS has not yet been measured, so the work_sim percentages must
not be projected to HG002 WGS as a confirmed value.

## Exactness check

`compare_profiles` fully decodes every read and compares every profile count:

```sh
fastkdual/compare_profiles baseline/reads     dual/reads
fastkdual/compare_profiles baseline/reads.ref dual/reads.ref
```

The work_sim validation additionally compares the histogram, table stub and
all hidden table/profile shards byte for byte against fresh `FastK` self and
relative runs made from the same `reads.dam`, `k=32`, and thread count.
`validation/worksim-t16-m1.sbatch` is the exact script used for Slurm job
422426 (`T=16`, `-M1`); its persistent log is
`worksim-t16-m1-422426.log`.
