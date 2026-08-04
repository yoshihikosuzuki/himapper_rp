# Development seed selectors

This directory contains experimental seed-selector variants that are not part
of the production build:

- `copy_minimizer`: copy-number-oriented minimizer variant
- `modimizer`: modulo-hash selector variant

Build or clean both programs from the repository root with:

```sh
make -C dev
make -C dev clean
```

The selector build flags can be overridden in the same way as for the root
Makefile, for example:

```sh
make -C dev MINIMAX_SELECTOR_FLAGS='-DMINIMAX_HASH_MOD=0 -DMINIMAX_HASH_REFINE=0 -DMINIMAX_ROBUST_WINNOW=1 -DMINIMAX_CANONICAL_WINNOW=1'
```
