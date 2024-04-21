#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <glib.h>

int* sieve_of_eratosthenes(int B, int* num_primes_under_B) {
    // Initialize array of booleans
    bool* is_prime = malloc((B + 1) * sizeof(bool));
    if (is_prime == NULL) {
        puts("ERROR: Unable to allocate memory for is_prime.");
        exit(1);
    }
    // TODO: Calloc is faster than malloc + loop
    // TODO: We can use a bit array instead of a boolean array
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

    // TODO: Some of these variables can just be ints, which may improve performance
    // If we switch over, we would need to implement our own modular exponentiation function using Pingala's algorithm
    mpz_t z, c, t, R, tmp, p_mpz, b;
    mpz_inits(z, c, t, R, tmp, p_mpz, b, NULL);
    mpz_set_ui(p_mpz, p);

    // 1. find S and Q such that p - 1 = Q * 2^S
    int Q = p - 1;
    unsigned long S = 0;
    while (Q % 2 == 0) {
        Q = Q / 2;
        S++;
    }

    // 2. find a quadratic non-residue z
    mpz_set_ui(z, 2);
    while (quadratic_residue(z, p) != p - 1) {
        mpz_add_ui(z, z, 1);
    }

    // 3. initialize M, c, t, R
    unsigned long M = S;
    mpz_powm_ui(c, z, Q, p_mpz);
    mpz_powm_ui(t, n, Q, p_mpz);
    mpz_set_ui(tmp, Q+1);
    mpz_fdiv_q_2exp(tmp, tmp, 1);
    mpz_powm_ui(R, n, mpz_get_ui(tmp), p_mpz);

    // 4. loop
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            *root_mod_p_1 = mpz_get_si(R);
            mpz_sub(tmp, p_mpz, R);
            *root_mod_p_2 = mpz_get_si(tmp);
            mpz_clears(z, c, t, R, tmp, p_mpz, b, NULL);
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

int* get_B_smooth_factors(mpz_t b, int* factor_base, int factor_base_size, int* factors_size) {
    // trial division to find the factors of b
    int* factors = malloc(factor_base_size * sizeof(int));
    if (factors == NULL) {
        puts("ERROR: Unable to allocate memory for factors.");
        exit(1);
    }

    *factors_size = 0;
    for (int i = 0; i < factor_base_size; i++) {
        int p = factor_base[i];
        while (mpz_divisible_ui_p(b, p) != 0) {
            mpz_divexact_ui(b, b, p); // fast since we know it's divisible
            factors[(*factors_size)++] = p;
        }
    }

    // Resize the array to the actual size needed
    if (*factors_size < factor_base_size) {
        int* resized_factors = realloc(factors, (*factors_size) * sizeof(int));
        factors = resized_factors;
    }

    return factors;
}

void compute_b(mpz_t b, int i, mpz_t root_n, mpz_t n) {
    mpz_set_si(b, i);
    mpz_add(b, b, root_n);

    // Calculate (i + root_n)^2
    mpz_mul(b, b, b);

    // Compute (i + root_n)^2 - n
    mpz_sub(b, b, n);
}

void get_factor_vector(int* factors, int factors_size, int* factor_base, int factor_base_size, int* exponent_vector) {
    // create a vector of the exponents of the factors in the factor base
    GHashTable* lookup_factor_index = g_hash_table_new(NULL, NULL);
    for (int i=0; i<factor_base_size; i++) {
        g_hash_table_insert(lookup_factor_index, GINT_TO_POINTER(factor_base[i]), GINT_TO_POINTER(i));
    }

    for (int i = 0; i < factors_size; i++) {
        int key = GPOINTER_TO_INT(g_hash_table_lookup(lookup_factor_index, GINT_TO_POINTER(factors[i])));
        exponent_vector[key] += 1;
    }

    g_hash_table_destroy(lookup_factor_index);
}

void create_matrix(double* sieve, int sieve_size, mpz_t root_n, int* factor_base, int factor_base_size, mpz_t n, GArray* matrix, GArray* as_vector, GHashTable* factor_exponent_dict) {
    double epsilon = 0.01;

    mpz_t b;
    mpz_init(b);

    for (int i=0; i<sieve_size; i++) {
        // TODO: Can do this in an earlier step
        if (sieve[i] < epsilon) {
            compute_b(b, i, root_n, n);
            int factors_size;
            int* factors = get_B_smooth_factors(b, factor_base, factor_base_size, &factors_size);

            int* exponent_vector = calloc(factor_base_size, sizeof(int));
            if (exponent_vector == NULL) {
                puts("ERROR: Unable to allocate memory for exponent_vector.");
                exit(1);
            }
            get_factor_vector(factors, factors_size, factor_base, factor_base_size, exponent_vector);

            bool* exponent_vector_mod_2 = malloc(factor_base_size * sizeof(bool));
            if (exponent_vector_mod_2 == NULL) {
                puts("ERROR: Unable to allocate memory for exponent_vector_mod_2.");
                exit(1);
            }
            int sum = 0;
            for (int j=0; j<factor_base_size; j++) {
                int mod_2 = exponent_vector[j] % 2;
                sum += mod_2;
                exponent_vector_mod_2[j] = mod_2;
            }
            if (sum == 0) {
                continue;
            }

            // matrix.append(exponent_vector_mod_2)
            g_array_append_val(matrix, exponent_vector_mod_2);

            // as_vector.append(i + root_n)
            mpz_t* i_plus_root_n = g_new(mpz_t, 1);
            mpz_init(*i_plus_root_n);
            mpz_add_ui(*i_plus_root_n, root_n, i);
            g_array_append_val(as_vector, i_plus_root_n);

            // factor_exponent_dict[i + root_n] = exponent_vector
            char* key_str = mpz_get_str(NULL, 10, *i_plus_root_n);
            g_hash_table_insert(factor_exponent_dict, key_str, exponent_vector);

            if (matrix->len >= 2*factor_base_size) {
                break;
            }

            free(factors);
        }
    }

    mpz_clear(b);
}

void find_linear_dependencies(GArray* dependencies, GArray* matrix, int factor_base_size) {
    // find the linear dependencies in the matrix
    // follows this paper https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf

    int n = matrix->len;
    int m = factor_base_size;

    // marks = np.zeros(n, dtype=bool)
    bool* marks = calloc(n, sizeof(bool));
    if (marks == NULL) {
        puts("ERROR: Unable to allocate memory for marks.");
        exit(1);
    }

    for (int i = 0; i < m; i++) {
        int piv = -1;
        for (int j = 0; j < n; j++) {
            bool* row = g_array_index(matrix, bool*, j);
            if (row[i] == 1) { // matrix[j][i] == 1
                marks[j] = true;
                piv = j;
                break;
            }
        }

        if (piv != -1) {
            bool* pivotRow = g_array_index(matrix, bool*, piv);
            for (int k=0; k<m; k++) {
                if (k != i && pivotRow[k] == 1) { // matrix[piv][k] == 1
                    for (int j=0; j<n; j++) {
                        bool* row = g_array_index(matrix, bool*, j);
                        if (row[i] == 1) { // If matrix[j][i] == 1,
                            row[k] = row[k] ^ 1; // Flip matrix[j][k]
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (!marks[i]) {
            GArray* dependent_list = g_array_new(FALSE, FALSE, sizeof(int));
            g_array_append_val(dependent_list, i);
            bool* row = g_array_index(matrix, bool*, i);
            
            for (int j = 0; j < m; j++) {
                if (row[j] == 1) {  // matrix[i][j] == 1
                    for (int k = 0; k < n; k++) {
                        if (g_array_index(matrix, bool*, k)[j] == 1) {  // matrix[k][j] == 1
                            g_array_append_val(dependent_list, k);
                            break;
                        }
                    }
                }
            }
            g_array_append_val(dependencies, dependent_list);
        }
    }

    free(marks);
}

void calculate_as_product(GArray* dependencies, GArray* exponent_as, mpz_t as_product) {
    mpz_set_ui(as_product, 1);
    for (int i = 0; i < dependencies->len; i++) {
        int row = g_array_index(dependencies, int, i);
        mpz_t* val = g_array_index(exponent_as, mpz_t*, row);
        mpz_mul(as_product, as_product, *val);
    }
}

void calculate_primes_product(GArray* dependencies, GHashTable* factor_exponent_dict, GArray* as_vector, int* factor_base, int factor_base_size, mpz_t primes_product) {
    mpz_set_ui(primes_product, 1);
    int* prime_power_vector = calloc(factor_base_size, sizeof(int));
    
    for (int i = 0; i < dependencies->len; i++) {
        int row = g_array_index(dependencies, int, i);

        mpz_t* key_mpz = g_array_index(as_vector, mpz_t*, row);
        char* key_str = mpz_get_str(NULL, 10, *key_mpz);
        int* factor_exponents = g_hash_table_lookup(factor_exponent_dict, key_str);

        for (int j = 0; j < factor_base_size; j++) {
            prime_power_vector[j] = prime_power_vector[j] + factor_exponents[j];
        }
    }

    // prime_power_vector = prime_power_vector // 2
    for (int i = 0; i < factor_base_size; i++) {
        prime_power_vector[i] = prime_power_vector[i] / 2;
    }

    mpz_t power;
    mpz_init(power);
    for (int i=0; i<factor_base_size; i++) {
        int p = factor_base[i];

        mpz_ui_pow_ui(power, p, prime_power_vector[i]);

        mpz_mul(primes_product, primes_product, power);
    }
    mpz_clear(power);
}

void euclidian_algorithm(mpz_t result, mpz_t a, mpz_t b) {
    mpz_t current_a, current_b;
    mpz_init_set(current_a, a);
    mpz_init_set(current_b, b);

    while (mpz_sgn(current_b) != 0) { // Continue while current_b is not zero
        mpz_t temp;
        mpz_init(temp);
        
        // temp = current_a % current_b
        mpz_mod(temp, current_a, current_b);
        
        // current_a = current_b
        mpz_set(current_a, current_b);
        
        // current_b = temp
        mpz_set(current_b, temp);
        
        mpz_clear(temp);
    }

    mpz_set(result, current_a);
    mpz_clear(current_a);
    mpz_clear(current_b);
}

bool return_factors(mpz_t factor1, mpz_t factor2, GArray* dependencies, GArray* as_vector, GHashTable* factor_exponent_dict, int* factor_base, int factor_base_size, mpz_t n) {
    if (dependencies->len == 0) {
        return 1;
    }
    
    for(int i = 0; i < dependencies->len; i++) { 
        GArray* dependency = g_array_index(dependencies, GArray*, i);

        mpz_t as_product;
        mpz_init(as_product);
        mpz_t primes_product;
        mpz_init(primes_product);
        
        calculate_as_product(dependency, as_vector, as_product);
        calculate_primes_product(dependency, factor_exponent_dict, as_vector, factor_base, factor_base_size, primes_product);

        mpz_t f;
        mpz_init(f);
        mpz_t a;
        mpz_init(a);
        mpz_sub(a, primes_product, as_product);
        euclidian_algorithm(f, a, n);
        if (mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, n) != 0) {
            mpz_set(factor1, f);
            mpz_divexact(factor2, n, f);
            return 0;
        }
    }

    return 1;
}

void sieve(mpz_t n, int B, int S, mpz_t factor1, mpz_t factor2) {
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

    /*
    matrix, as_vector, factor_exponent_dict = create_matrix(sieve, root_n, factor_base, n)
    */
    GArray* matrix = g_array_new(FALSE, FALSE, sizeof(bool*));
    GArray* as_vector = g_array_new(FALSE, FALSE, sizeof(mpz_t*));
    GHashTable* factor_exponent_dict = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    create_matrix(sieve, S, root_n, factor_base, factor_base_size, n, matrix, as_vector, factor_exponent_dict);

    // dependencies = find_linear_dependencies(matrix_rr)
    GArray* dependencies = g_array_new(FALSE, FALSE, sizeof(GArray*));
    find_linear_dependencies(dependencies, matrix, factor_base_size);

    // return return_factors(dependencies, as_vector, factor_exponent_dict, factor_base, n)
    bool fail = return_factors(factor1, factor2, dependencies, as_vector, factor_exponent_dict, factor_base, factor_base_size, n);
    if (fail) {
        puts("ERROR: No nontrivial factors found.");
    }

    /*
    Free memory
    */
    // matrix
    for (int i = 0; i < matrix->len; i++) {
        free(g_array_index(matrix, bool*, i));
    }
    g_array_free(matrix, TRUE);

    // as_vector
    for (int i = 0; i < as_vector->len; i++) {
        mpz_t* value = g_array_index(as_vector, mpz_t*, i);
        mpz_clear(*value);
        g_free(value);
    }
    g_array_free(as_vector, TRUE);

    // dependencies
    for (int i = 0; i < dependencies->len; i++) {
        GArray* inner_array = g_array_index(dependencies, GArray*, i);
        g_array_free(inner_array, TRUE);
    }
    g_array_free(dependencies, TRUE);

    g_hash_table_destroy(factor_exponent_dict);

    free(primes_under_B);
    free(factor_base);
    free(sieve);
    mpz_clear(root_n);
}

int main() {
    mpz_t n;
    mpz_init(n);

    mpz_set_str(n, "16921456439215439701", 10); // base 10
    int B = 2000;
    int S = 4000000;

    // mpz_set_str(n, "46839566299936919234246726809", 10); // base 10
    // int B = 15000;
    // int S = 15000000;

    // mpz_set_str(n, "6172835808641975203638304919691358469663", 10); // base 10
    // int B = 30000;
    // int S = 1000000000;

    // Nontrivial factors of n
    mpz_t factor1;
    mpz_init(factor1);
    mpz_t factor2;
    mpz_init(factor2);

    sieve(n, B, S, factor1, factor2);
    gmp_printf("n: %Zd, factors: (%Zd, %Zd)\n", n, factor1, factor2);

    // Clear memory
    mpz_clear(n);
    mpz_clear(factor1);
    mpz_clear(factor2);

    return 0;
}
