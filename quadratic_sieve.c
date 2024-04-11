#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>

int* sieve_of_eratosthenes(int B, int* factor_base_length) {
    // Dynamic array to store primes
    int capacity = 1;
    int* factor_base = (int*)malloc(capacity * sizeof(int));
    *factor_base_length = 0;

    // Allocate memory for prime check array
    bool* is_primes = (bool*)malloc((B + 1) * sizeof(bool));
    for(int i = 0; i <= B; i++) {
        is_primes[i] = true;
    }
    is_primes[0] = false;
    is_primes[1] = false;

    // Sieve of Eratosthenes
    for(int i = 2; i < ceil(sqrt(B+1)); i++) {
        if(is_primes[i]) {
            // Dynamically resize factor_base if necessary
            if (*factor_base_length >= capacity) {
                capacity *= 2;
                factor_base = (int*)realloc(factor_base, capacity * sizeof(int));
            }

            factor_base[*factor_base_length] = i;
            (*factor_base_length)++;

            for(int j = i * i; j <= B; j += i) {
                is_primes[j] = false;
            }
        }
    }

    for (int i = ceil(sqrt(B+1)); i <= B; i++) {
        if(is_primes[i]) {
            // Dynamically resize factor_base if necessary
            if (*factor_base_length >= capacity) {
                capacity *= 2;
                factor_base = (int*)realloc(factor_base, capacity * sizeof(int));
            }

            factor_base[*factor_base_length] = i;
            (*factor_base_length)++;
        }
    }

    free(is_primes);

    return factor_base;
}

void sieve(mpz_t n, int B, int S, mpz_t* factor1, mpz_t* factor2) {
    /*
        ceil(sqrt(n))
    */
    mpz_t root_n;
    mpz_init(root_n);

    mpf_t n_f;
    mpf_init(n_f);
    mpf_set_z(n_f, n);

    mpf_sqrt(n_f, n_f);
    mpf_ceil(n_f, n_f);
    mpz_set_f(root_n, n_f);
    mpf_clear(n_f);

    /*
        sieve_of_eratosthenes
    */
    puts("Creating sieve...");
    int factor_base_length = 0;
    int* primes_under_B = sieve_of_eratosthenes(B, &factor_base_length);

    printf("Primes up to %d:\n", B);
    for(int i = 0; i < factor_base_length; i++) {
        printf("%d ", primes_under_B[i]);
    }
    printf("\n");

    free(primes_under_B);

    // Trivial factors (delete later)
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
