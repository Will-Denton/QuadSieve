set -e

gcc quadratic_sieve.c -o quadratic_sieve -I /opt/homebrew/opt/gmp/include -L /opt/homebrew/opt/gmp/lib -lgmp -I /opt/homebrew/opt/glib/include/glib-2.0 -I /opt/homebrew/opt/glib/lib/glib-2.0/include -L /opt/homebrew/opt/glib/lib -lglib-2.0 -Ofast
time ./quadratic_sieve