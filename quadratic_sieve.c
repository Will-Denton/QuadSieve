#include <stdio.h>
#include <gmp.h>

void sieve(mpz_t n, int B, int S, mpz_t* factor1, mpz_t* factor2) {
    // Trivial factors
    mpz_set_ui(*factor1, 1);
    mpz_set(*factor2, n);
}

int main() {
    mpz_t n;
    mpz_init(n);
    mpz_set_str(n, "16921456439215439701", 10); // base 10

    int B = 15000;
    int S = 15000000;

    // Factors of n
    mpz_t factor1;
    mpz_init(factor1);
    mpz_t factor2;
    mpz_init(factor2);

    sieve(n, B, S, &factor1, &factor2);
    gmp_printf("n: %Zd, factors: (%Zd, %Zd)\n", n, factor1, factor2);

    // Clear memory
    mpz_clear(n);
    mpz_clear(factor1);
    mpz_clear(factor2);

    return 0;
}
