set -e

gcc quadratic_sieve.c -o quadratic_sieve -O3 -lgmp -lm
time ./quadratic_sieve