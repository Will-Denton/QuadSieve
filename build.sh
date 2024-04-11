set -e

gcc quadratic_sieve.c -o quadratic_sieve -I$(brew --prefix gmp)/include -L$(brew --prefix gmp)/lib -lgmp
./quadratic_sieve