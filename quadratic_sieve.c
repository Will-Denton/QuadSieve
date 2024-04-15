#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>

int* sieve_of_eratosthenes(int B, int* num_primes_under_B) {
    // Initialize array of booleans
    bool* is_prime = malloc((B + 1) * sizeof(bool));
    for (int i = 0; i < B + 1; i++) {
        is_prime[i] = true;
    }
    is_prime[0] = is_prime[1] = false;

    // Sieve of Eratosthenes
    for (int i = 2; i < ((int)sqrt(B))+1; i++) {
        if (is_prime[i]) {
            for (int j = i * i; j < B+1; j += i) {
                is_prime[j] = false;
            }
        }
    }

    // Count number of primes under B
    *num_primes_under_B = 0;
    for (int i = 2; i < B+1; i++) {
        if (is_prime[i]) {
            (*num_primes_under_B)++;
        }
    }

    // Get primes
    int* primes_under_B = malloc(*num_primes_under_B * sizeof(int));
    int primes_index = 0;
    for (int i = 2; i < B + 1; i++) {
        if (is_prime[i]) {
            primes_under_B[primes_index++] = i;
        }
    }

    // Free memory
    free(is_prime);

    return primes_under_B;
}

int quadratic_residue(mpz_t n, int p) {
    mpz_t mpz_p;
    mpz_init_set_ui(mpz_p, p);  // Convert int to mpz_t

    if (mpz_cmp_ui(mpz_p, 2) == 0) {
        return 1;
    }

    mpz_t exp, result;
    mpz_inits(exp, result, NULL);

    // Calculate exp = (p - 1) // 2
    mpz_sub_ui(exp, mpz_p, 1);
    mpz_fdiv_q_ui(exp, exp, 2);

    // Calculate n^exp (mod p)
    mpz_powm(result, n, exp, mpz_p);
    int res = (int)mpz_get_si(result);  // NOTE: This cast is only safe since p is an int

    mpz_clears(exp, result, NULL);
    return res;
}

mpz_t* get_factor_base(int* primes, int num_primes, mpz_t n, int* factor_base_size) {
    // NOTE: We are temporarily allocating more memory than we need here
    mpz_t* factor_base = malloc(num_primes * sizeof(mpz_t));

    *factor_base_size = 0;
    for (int i = 0; i < num_primes; i++) {
        if (quadratic_residue(n, primes[i]) == 1) {
            mpz_init_set_ui(factor_base[(*factor_base_size)++], primes[i]);
        }
    }

    // Resize the array to the actual size needed
    if (*factor_base_size < num_primes) {
        mpz_t* resized_factor_base = realloc(factor_base, (*factor_base_size) * sizeof(mpz_t));
        factor_base = resized_factor_base;
    }

    return factor_base;
}

double* get_sieve_log(int S, mpz_t n) {
    puts("Creating Sieve Log...");

    // root_n = np.float64(np.ceil(math.sqrt(n)))
    mpz_t root_n_mpz;
    mpz_init(root_n_mpz);
    mpz_sqrt(root_n_mpz, n);
    if (mpz_perfect_square_p(n) == 0) { // if n is not a perfect square
        mpz_add_ui(root_n_mpz, root_n_mpz, 1);
    }
    double root_n = mpz_get_d(root_n_mpz);

    // Create sieve
    double* sieve = malloc(S * sizeof(double));
    for (int i = 0; i < S; i++) {
        double current = i + root_n;
        double value = current * current - mpz_get_d(n);
        sieve[i] = log(value);
    }

    // Clean up
    free(sieve);
    mpz_clear(root_n_mpz);

    return sieve;
}


void sieve(mpz_t n, int B, int S, mpz_t* factor1, mpz_t* factor2) {
    /*
    root_n = math.ceil(math.sqrt(n))
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
    primes_under_B = sieve_of_eratosthenes(B)
    */
    int num_primes_under_B;
    int* primes_under_B = sieve_of_eratosthenes(B, &num_primes_under_B);

    /*
    factor_base = get_factor_base(primes_under_B, n)
    */
    int factor_base_size;
    mpz_t* factor_base = get_factor_base(primes_under_B, num_primes_under_B, n, &factor_base_size);

    /*
    sieve = get_sieve_log(S, n)
    */
    // NOTE: We get slightly different results here than the Python code
    // This might be due to floating point precision issues
    double* sieve = get_sieve_log(S, n);

    /*
    sieve_primes_log(n, factor_base, S, sieve)
    */

    // Free memory
    free(primes_under_B);

    // Trivial factors (delete later)
    mpz_set_ui(*factor1, 1);
    mpz_set(*factor2, n);
}

int main() {
    mpz_t n;
    mpz_init(n);
    mpz_set_str(n, "46839566299936919234246726809", 10); // base 10

    int B = 15000;
    int S = 15000000;

    // Nontrivial factors of n
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
