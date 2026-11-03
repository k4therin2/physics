#!/bin/bash
# runs the C_l^(22) evaluations (single-threaded BLAS; 4-core box)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
run() { python3 cl22.py "$@" > cl22_$7.log 2>&1; }
( run 1 2 8 24 40 24 eo8_a; run 1 2 4 24 40 24 eo4_a ) &
( run 1 2 8 40 48 24 eo8_b; run 1 2 4 40 48 24 eo4_b ) &
( run 1 2 8 56 60 24 eo8_c ) &
wait
