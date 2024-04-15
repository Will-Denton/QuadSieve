#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>

int* sieve_of_eratosthenes(int B, int* num_primes_under_B) {
    // Initialize array of booleans
    bool* is_prime = malloc((B + 1) * sizeof(bool));
    if (is_prime == NULL) {
        puts("ERROR: Unable to allocate memory for is_prime.");
        exit(1);
    }
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
    if (primes_under_B == NULL) {
        puts("ERROR: Unable to allocate memory for primes_under_B.");
        exit(1);
    }
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
    if (p == 0) {
        return 1;
    }

    mpz_t mpz_p;
    mpz_init_set_ui(mpz_p, p);  // Convert int to mpz_t

    mpz_t exp, result;
    mpz_inits(exp, result, NULL);

    // Calculate exp = (p - 1) // 2
    mpz_sub_ui(exp, mpz_p, 1);
    mpz_fdiv_q_ui(exp, exp, 2);

    // Calculate n^exp (mod p)
    mpz_powm(result, n, exp, mpz_p);
    int res = (int)mpz_get_si(result);  // NOTE: This cast is only safe since p is an int

    mpz_clears(mpz_p, exp, result, NULL);
    return res;
}

int* get_factor_base(int* primes, int num_primes, mpz_t n, int* factor_base_size) {
    // NOTE: We are temporarily allocating more memory than we need here
    int* factor_base = malloc(num_primes * sizeof(int));
    if (factor_base == NULL) {
        puts("ERROR: Unable to allocate memory for factor_base.");
        exit(1);
    }

    *factor_base_size = 0;
    for (int i = 0; i < num_primes; i++) {
        if (quadratic_residue(n, primes[i]) == 1) {
            factor_base[(*factor_base_size)++] = primes[i];
        }
    }

    // Resize the array to the actual size needed
    if (*factor_base_size < num_primes) {
        int* resized_factor_base = realloc(factor_base, (*factor_base_size) * sizeof(int));
        factor_base = resized_factor_base;
    }

    return factor_base;
}

// Computes the ceiling of the square root of n
void ceil_sqrt(mpz_t result, const mpz_t n) {
    mpz_init(result);
    mpz_sqrt(result, n);
    if (mpz_perfect_square_p(n) == 0) {  // If n is not a perfect square
        mpz_add_ui(result, result, 1);
    }
}

double* get_sieve_log(int S, mpz_t n) {
    puts("Starting get_sieve_log...");

    // root_n = np.float64(np.ceil(math.sqrt(n)))
    mpz_t root_n_mpz;
    ceil_sqrt(root_n_mpz, n);
    double root_n = mpz_get_d(root_n_mpz);

    // Create sieve
    double* sieve = malloc(S * sizeof(double));
    if (sieve == NULL) {
        puts("ERROR: Unable to allocate memory for sieve.");
        exit(1);
    }

    for (int i = 0; i < S; i++) {
        double current = i + root_n;
        double value = current * current - mpz_get_d(n);
        sieve[i] = log(value);
    }

    // Clean up
    mpz_clear(root_n_mpz);

    return sieve;
}

void shanks_tonelli(mpz_t n, int p, int *root_mod_p_1, int *root_mod_p_2) {
    // calculate the square roots of n mod p
    // DID THIS https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    if (p == 2) {
        *root_mod_p_1 = 1;
        return;
    }

    mpz_t Q, z, c, t, R, tmp, p_mpz, b;
    mpz_inits(Q, z, c, t, R, tmp, p_mpz, b, NULL);
    mpz_set_ui(p_mpz, p);

    // 1. find S and Q such that p - 1 = Q * 2^S
    mpz_sub_ui(Q, p_mpz, 1);
    unsigned long S = 0;
    while (mpz_even_p(Q)) {
        mpz_fdiv_q_2exp(Q, Q, 1);
        S++;
    }

    // 2. find a quadratic non-residue z
    mpz_set_ui(z, 2);
    while (quadratic_residue(z, p) != p - 1) {
        mpz_add_ui(z, z, 1);
    }

    // 3. initialize M, c, t, R
    unsigned long M = S;
    mpz_powm_ui(c, z, mpz_get_ui(Q), p_mpz);
    mpz_powm_ui(t, n, mpz_get_ui(Q), p_mpz);
    mpz_add_ui(tmp, Q, 1);
    mpz_fdiv_q_2exp(tmp, tmp, 1);
    mpz_powm_ui(R, n, mpz_get_ui(tmp), p_mpz);

    // 4. loop
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            *root_mod_p_1 = mpz_get_si(R);
            mpz_sub(tmp, p_mpz, R);
            *root_mod_p_2 = mpz_get_si(tmp);
            mpz_clears(Q, z, c, t, R, tmp, p_mpz, b, NULL);
            return;
        }

        unsigned long i = 1;
        mpz_powm_ui(tmp, t, 1 << i, p_mpz); // 2**i = 1 << i
        while (mpz_cmp_ui(tmp, 1) != 0) {
            i++;
            mpz_powm_ui(tmp, t, 1 << i, p_mpz);
        }

        mpz_powm_ui(b, c, 1 << (M - i - 1), p_mpz);
        M = i;
        mpz_powm_ui(c, b, 2, p_mpz);

        // t = (t * pow(b, 2, p)) % p
        mpz_powm_ui(tmp, b, 2, p_mpz);
        mpz_mul(t, t, tmp);
        mpz_mod(t, t, p_mpz);

        mpz_mul(R, R, b);
        mpz_mod(R, R, p_mpz);
    }
}

void sieve_primes_log(mpz_t n, int* factor_base, int factor_base_size, int S, double* log_sieve) {
    puts("Starting sieve_primes_log...");

    mpz_t root_n, tmp, p_mpz, root_mod_p_1_mpz, root_mod_p_2_mpz;
    ceil_sqrt(root_n, n);
    mpz_inits(tmp, p_mpz, root_mod_p_1_mpz, root_mod_p_2_mpz, NULL);

    for (int i=0; i<factor_base_size; i++) {
        int p = factor_base[i];
        mpz_set_ui(p_mpz, p);

        int root_mod_p_1 = -1;
        int root_mod_p_2 = -1;
        shanks_tonelli(n, p, &root_mod_p_1, &root_mod_p_2);
        if (root_mod_p_1 != -1) {
            // Calculate x = (root_mod_p_1 - root_n) % p
            mpz_set_ui(root_mod_p_1_mpz, root_mod_p_1);
            mpz_sub(tmp, root_mod_p_1_mpz, root_n);
            mpz_mod(tmp, tmp, p_mpz);
            int x = mpz_get_si(tmp);

            for (int j=x; j<S; j+=p) {
                log_sieve[j] -= log(p);
            }
        }
        if (root_mod_p_2 != -1) {
            // Calculate x = (root_mod_p_2 - root_n) % p
            mpz_set_ui(root_mod_p_2_mpz, root_mod_p_2);
            mpz_sub(tmp, root_mod_p_2_mpz, root_n);
            mpz_mod(tmp, tmp, p_mpz);
            int x = mpz_get_si(tmp);

            for (int j=x; j<S; j+=p) {
                log_sieve[j] -= log(p);
            }
        }
    }

    mpz_clears(root_n, tmp, p_mpz, root_mod_p_1_mpz, root_mod_p_2_mpz, NULL);
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
    int* factor_base = get_factor_base(primes_under_B, num_primes_under_B, n, &factor_base_size);

    /*
    sieve = get_sieve_log(S, n)
    */
    // NOTE: We get slightly different results here than the Python code
    // This might be due to floating point precision issues
    double* sieve = get_sieve_log(S, n); // size S

    /*
    sieve_primes_log(n, factor_base, S, sieve)
    */
    sieve_primes_log(n, factor_base, factor_base_size, S, sieve);

    // Free memory
    free(primes_under_B);
    free(factor_base);
    free(sieve);
    mpz_clear(root_n);

    // Trivial factors (delete later)
    mpz_set_ui(*factor1, 1);
    mpz_set(*factor2, n);
}

int main() {
    mpz_t n;
    mpz_init(n);

    // mpz_set_str(n, "46839566299936919234246726809", 10); // base 10
    // int B = 15000;
    // int S = 15000000;

    mpz_set_str(n, "6172835808641975203638304919691358469663", 10); // base 10
    int B = 30000;
    int S = 1000000000;

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
