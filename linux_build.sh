set -e

gcc quadratic_sieve.c -o quadratic_sieve \
-Ofast \
-lgmp \
-lm \
`pkg-config --cflags --libs glib-2.0` \
`pkg-config --cflags --libs mpfr`

valgrind --leak-check=full ./quadratic_sieve