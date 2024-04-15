set -e

gcc quadratic_sieve.c -o quadratic_sieve -I /opt/homebrew/opt/gmp/include -L /opt/homebrew/opt/gmp/lib -l gmp -O3
time ./quadratic_sieve