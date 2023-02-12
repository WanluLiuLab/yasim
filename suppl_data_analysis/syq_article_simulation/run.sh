#!/usr/bin/env bash
set -ue
python -m yasim llr


mkdir -p diff_read_depth
cd diff_read_depth || exit 1

# Read depth
for mean_depth in 10 25 40 55 70 85 100; do
    pass
done

for complexity_level in 1 3 5 7 9; do
    pass
done

for accuracy in 0.6 0.7 0.8 0.9 1.0; do
    pass
done

for five_prime_break in ; do
    pass
done

# NIpG
# RLEN
# ERROR RATE
# RC

