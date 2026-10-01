# Development seed selectors

This directory contains experimental seed-selector variants that are not part
of the production build:

- `copy_minimizer`: copy-number-oriented minimizer variant
- `modimizer`: modulo-hash selector variant
- `dump_kpos`: block-indexed KPOS container inspector and ASCII dumper

Build or clean both programs from the repository root with:

```sh
make -C dev
make -C dev clean
```

The standalone command lines are:

```text
copy_minimizer [-T<N>] [-k<N>] [-w<N>] [-s<N>] [-l<N>] <db> <genome.prof> <out.kpos>
modimizer      [-T<N>] [-k<N>]                   [-l<N>] <db> <genome.prof> <out.kpos>
dump_kpos <file.kpos> [block]
```

`MINIMAX_SELECTOR_FLAGS` applies only to `copy_minimizer`.  It can be
overridden in the same way as in the root Makefile, for example:

```sh
make -C dev copy_minimizer MINIMAX_SELECTOR_FLAGS='-DMINIMAX_HASH_MOD=0 -DMINIMAX_HASH_REFINE=0 -DMINIMAX_ROBUST_WINNOW=1 -DMINIMAX_CANONICAL_WINNOW=1'
```

## Mapper algorithm experiments

[mapper_experiments/](mapper_experiments/README.md) contains the C adapter for applying
minimap2 chaining/alignment functions to himapper seeds. The experiment-specific
integration, frozen dependencies, build products and results are in
[work_eval](../work_eval/README.md). The regular mapper build does not depend on this adapter.

## MAPQ tests

[mapq_tests/](mapq_tests/README.md) contains standalone score/catalog tests and
read-block bookkeeping tests. Real-read compatibility and resource measurements
are kept under `work_eval/experiments/20260925_mapq/`.
