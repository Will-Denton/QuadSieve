set -e

gcc quadratic_sieve.c -o quadratic_sieve -O3 -lgmp -lm `pkg-config --cflags --libs glib-2.0`
valgrind ./quadratic_sieve