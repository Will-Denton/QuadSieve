#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <glib.h>
#include <assert.h>
#include <mpfr.h>
#include <pthread.h>

int* sieve_of_eratosthenes(int B, int* num_primes_under_B) {
    /*
    Runs the sieve of eratosthenes algorithm to find all primes < B.
    */

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
    /*
    Calculates the legendre symbol (n/p) = n^((p-1)/2) mod p.
    */
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
    /*
    Creates the factor base for the quadratic sieve.
    In order for shanks-tonelli to work, the factor base must contain p such that n is a quadratic residue mod p.
    Limit the primes found from the sieve of eratosthenes to only those that are quadratic residues.
    */

   // allocate memory for factor_base
    int* factor_base = malloc(num_primes * sizeof(int));
    if (factor_base == NULL) {
        puts("ERROR: Unable to allocate memory for factor_base.");
        exit(1);
    }

    // Find the size of the factor base
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


void ceil_sqrt(mpz_t result, const mpz_t n) {
    /*
    Computes the ceiling of the square root of n.
    */
    mpz_init(result);
    mpz_sqrt(result, n);
    if (mpz_perfect_square_p(n) == 0) {  // If n is not a perfect square
        mpz_add_ui(result, result, 1);
    }
}


// struct to store data local to the threads running get_sieve_log
typedef struct {
    int start;
    int end;
    __float128 start_val;
    double n_double;
    double* sieve;
} ThreadData;



void* compute_sieve(void* arg) {
    /*
    calculates ln((i + root_n)^2 - n) for each thrads given range.
    */
    ThreadData* data = (ThreadData*)arg;
    __float128 current = data->start_val; // use __float128 to avoid floating point errors
    for (int i = data->start; i < data->end; i++, current += 1.0) {
        __float128 value = current * current - data->n_double;
        data->sieve[i] = (double) log(value);
    }
    return NULL;
}

double* get_sieve_log(int S, mpz_t n) {
    /*
    Given a bound S, create a list of numbers to sieve.
    The sieve list is given by ln((i + root_n)^2 - n) for i in range(S).
    All floats are 64 bit since numpy list compression cant use anything larger.
    This causes a lot of error to accumulate, so epsilon is used to offset this later.
    */
    puts("Starting get_sieve_log...");
    mpz_t root_n_mpz;
    ceil_sqrt(root_n_mpz, n);
    double root_n = mpz_get_d(root_n_mpz);
    mpz_clear(root_n_mpz);
    double n_double = mpz_get_d(n);

    // Create sieve
    double* sieve = malloc(S * sizeof(double));
    if (sieve == NULL) {
        puts("ERROR: Unable to allocate memory for sieve.");
        exit(1);
    }

    // Create threads
    int num_threads = 10;
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];
    int length_per_thread = S / num_threads;
    
    __float128 current = root_n;
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].start = i * length_per_thread;
        thread_data[i].end = (i == num_threads - 1) ? S : (i + 1) * length_per_thread;
        thread_data[i].start_val = current + i * length_per_thread;
        thread_data[i].n_double = n_double;
        thread_data[i].sieve = sieve;
        pthread_create(&threads[i], NULL, compute_sieve, &thread_data[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    return sieve;
}


void shanks_tonelli(mpz_t n, int p, int *root_mod_p_1, int *root_mod_p_2) {
    /*
    Standard shanks-tonelli algorithm to find the square roots of n mod p.
    */
    if (p == 2) {
        *root_mod_p_1 = 1;
        return;
    }

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
    /*
    Sieve the numbers in the sieve that are divisible by primes in the factor base.
    Shanks-tonelli is used to find the square roots of n mod p and solutions to this
    are used to find the numbers in the sieve that are divisible by p in the sieve.
    All values in the sieve that are about 0 after this process are B-smooth.
    */
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
            mpz_set_ui(root_mod_p_1_mpz, root_mod_p_1);
            mpz_sub(tmp, root_mod_p_1_mpz, root_n);
            mpz_mod(tmp, tmp, p_mpz);
            int x = mpz_get_si(tmp);
            for (int j=x; j<S; j+=p) {
                log_sieve[j] -= log(p);
            }
        }
        if (root_mod_p_2 != -1) {
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


void get_B_smooth_factors(mpz_t b, int* factor_base, int factor_base_size, GArray* factors) {
    /*
    The sieve process does not find the factors of b, only that b is B-smooth.
    This function finds the factors of b that are in the factor base using trial division.
    */
    for (int i = 0; i < factor_base_size; i++) {
        int p = factor_base[i];
        while (mpz_divisible_ui_p(b, p) != 0 && mpz_cmp_ui(b, 0) != 0) {
            mpz_divexact_ui(b, b, p); // fast since we know it's divisible
            g_array_append_val(factors, p);
        }
    }
}

void compute_b(mpz_t b, int i, mpz_t root_n, mpz_t n) {
    /*
    Calculates b
    */
    mpz_set_si(b, i);
    mpz_add(b, b, root_n);

    // Calculate (i + root_n)^2
    mpz_mul(b, b, b);

    // Compute (i + root_n)^2 - n
    mpz_sub(b, b, n);
}

void get_factor_vector(GArray* factors, int* factor_base, int factor_base_size, int* exponent_vector) {
    /*
    Creates the exponent vector for the factors of b in the factor base.
    The exponent vector is given as the number of times each factor in the factor base divides b.
    */
    GHashTable* lookup_factor_index = g_hash_table_new(NULL, NULL);
    for (int i=0; i<factor_base_size; i++) {
        g_hash_table_insert(lookup_factor_index, GINT_TO_POINTER(factor_base[i]), GINT_TO_POINTER(i));
    }

    for (int i = 0; i < factors->len; i++) {
        int factor = g_array_index(factors, int, i);
        int key = GPOINTER_TO_INT(g_hash_table_lookup(lookup_factor_index, GINT_TO_POINTER(factor)));
        exponent_vector[key] += 1;
    }

    g_hash_table_destroy(lookup_factor_index);
}

void create_matrix_log(double* sieve, int sieve_size, mpz_t root_n, int* factor_base, int factor_base_size, mpz_t n, GArray* matrix, GArray* as_vector, GHashTable* factor_exponent_dict) {
    /*
    Creates the matrix mod 2 for finding dependencies in the sieve.
    The matrix is created by finding the B-smooth numbers in the sieve and creating the exponent vectors for each a value.
    The goal is to have a matrix with more rows than columns to find linear dependencies.
    */
    puts("Starting create_matrix_log...");

    double epsilon = 0.0001;
    mpz_t b;
    mpz_init(b);

    for (int i=0; i<sieve_size; i++) {
        if (sieve[i] < epsilon) {
            compute_b(b, i, root_n, n);
            GArray* factors = g_array_new(FALSE, FALSE, sizeof(int));
            if (factors == NULL) {
                puts("ERROR: Unable to allocate GArray for factors.");
                exit(1);
            }
            get_B_smooth_factors(b, factor_base, factor_base_size, factors);

            int* exponent_vector = calloc(factor_base_size, sizeof(int));
            if (exponent_vector == NULL) {
                puts("ERROR: Unable to allocate memory for exponent_vector.");
                exit(1);
            }
            get_factor_vector(factors, factor_base, factor_base_size, exponent_vector);

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

            g_array_append_val(matrix, exponent_vector_mod_2);

            mpz_t* i_plus_root_n = g_new(mpz_t, 1);
            mpz_init(*i_plus_root_n);
            mpz_add_ui(*i_plus_root_n, root_n, i);
            g_array_append_val(as_vector, i_plus_root_n);

            char* key_str = mpz_get_str(NULL, 10, *i_plus_root_n);
            g_hash_table_insert(factor_exponent_dict, key_str, exponent_vector);

            if (matrix->len >= 2*factor_base_size) {
                break;
            }

            g_array_free(factors, TRUE);
        }
    }

    mpz_clear(b);
}


void find_linear_dependencies(GArray* dependencies, GArray* matrix, int factor_base_size) {
    /*
    Finds the linear dependencies in the matrix using gaussian elimination mod 2.
    Follows the algorithm presented in this paper: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    */
    puts("Starting find_linear_dependencies...");

    int n = matrix->len;
    int m = factor_base_size;

    bool* marks = calloc(n, sizeof(bool));
    if (marks == NULL) {
        puts("ERROR: Unable to allocate memory for marks.");
        exit(1);
    }

    /*
    algorithm presented in the paper:
    1. Search for Aij = 1 in column j
    2. If found then mark row i as an independent row
    3. For all other columns k != j, if Akj = 1, then xor column k with column j
    4. repeat for all columns
    */
    for (int i = 0; i < m; i++) {
        int piv = -1;
        for (int j = 0; j < n; j++) {
            bool* row = g_array_index(matrix, bool*, j);
            if (row[i] == 1) { 
                marks[j] = true;
                piv = j;
                break;
            }
        }
        if (piv != -1) {
            bool* pivotRow = g_array_index(matrix, bool*, piv);
            for (int k=0; k<m; k++) {
                if (k != i && pivotRow[k] == 1) { 
                    for (int j=0; j<n; j++) {
                        bool* row = g_array_index(matrix, bool*, j);
                        if (row[i] == 1) { 
                            row[k] = row[k] ^ 1;
                        }
                    }
                }
            }
        }
    }

    /*
    All unmarked rows are dependent on the marked rows
    loop through all unmarked rows and create a list of rows that the unmarked row is dependent on
    */
    for (int i = 0; i < n; i++) {
        if (!marks[i]) {
            GArray* dependent_list = g_array_new(FALSE, FALSE, sizeof(int));
            g_array_append_val(dependent_list, i);
            bool* row = g_array_index(matrix, bool*, i);
            
            for (int j = 0; j < m; j++) {
                if (row[j] == 1) {  
                    for (int k = 0; k < n; k++) {
                        if (g_array_index(matrix, bool*, k)[j] == 1) {  
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
    /*
    Calculates the product of the as for (a1...an)^2 = (p1^e1...pn^en)^2 mod n.
    */
    mpz_set_ui(as_product, 1);
    for (int i = 0; i < dependencies->len; i++) {
        int row = g_array_index(dependencies, int, i);
        mpz_t* val = g_array_index(exponent_as, mpz_t*, row);
        mpz_mul(as_product, as_product, *val);
    }
}


void calculate_primes_product(GArray* dependencies, GHashTable* factor_exponent_dict, GArray* as_vector, int* factor_base, int factor_base_size, mpz_t primes_product) {
    /*
    Calculates the product of the primes for (a1...an)^2 = (p1^e1...pn^en)^2 mod n.
    */
    mpz_set_ui(primes_product, 1);
    int* prime_power_vector = calloc(factor_base_size, sizeof(int));
    
    for (int i = 0; i < dependencies->len; i++) {
        int row = g_array_index(dependencies, int, i);

        mpz_t* key_mpz = g_array_index(as_vector, mpz_t*, row);
        char* key_str = mpz_get_str(NULL, 10, *key_mpz);
        int* factor_exponents = g_hash_table_lookup(factor_exponent_dict, key_str);
        free(key_str);

        for (int j = 0; j < factor_base_size; j++) {
            prime_power_vector[j] = prime_power_vector[j] + factor_exponents[j];
        }
    }

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
    free(prime_power_vector);
}


void euclidian_algorithm(mpz_t result, mpz_t a, mpz_t b) {
    /*
    Finds the greatest common divisor of a and b using the euclidian algorithm.
    */
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


void return_factors(mpz_t factor1, mpz_t factor2, GArray* dependencies, GArray* as_vector, GHashTable* factor_exponent_dict, int* factor_base, int factor_base_size, mpz_t n) {
    /*
    Returns the factors of n by calculating the product of the as and the product of the primes.
    Uses the basic principle that given (a1...an)^2 = (p1^e1...pn^en)^2 mod n, then a factor is = gcd(a1...an - p1^e1...pn^en, n).
    */
    if (dependencies->len == 0) {
        return;
    }

    mpz_t as_product, primes_product, f, a;
    mpz_inits(as_product, primes_product, f, a, NULL);

    for(int i = 0; i < dependencies->len; i++) { 
        GArray* dependency = g_array_index(dependencies, GArray*, i);

        calculate_as_product(dependency, as_vector, as_product);
        calculate_primes_product(dependency, factor_exponent_dict, as_vector, factor_base, factor_base_size, primes_product);

        mpz_sub(a, primes_product, as_product);
        euclidian_algorithm(f, a, n);
        if (mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, n) != 0) {
            mpz_set(factor1, f);
            mpz_divexact(factor2, n, f);

            break;
        }
    }

    mpz_clears(as_product, primes_product, f, a, NULL);
}


void sieve_log(mpz_t n, int B, int S, mpz_t factor1, mpz_t factor2) {
    /*
    Main sieve function for the quadratic sieve.
    */

    // Check to see if n is a perfect square
    if (mpz_perfect_square_p(n) != 0) {
        gmp_printf("n: %Zd is a perfect square.\n", n);
        mpz_sqrt(factor1, n);
        mpz_set(factor2, factor1);
        return;
    }

    mpz_t root_n;
    mpz_init(root_n);
    ceil_sqrt(root_n, n);

    // Initialize the factor base
    int num_primes_under_B;
    int* primes_under_B = sieve_of_eratosthenes(B, &num_primes_under_B);
    
    // Get the true factor base
    int factor_base_size;
    int* factor_base = get_factor_base(primes_under_B, num_primes_under_B, n, &factor_base_size);

    // Get the sieve
    double* sieve = get_sieve_log(S, n); // size S

    // Sieve the primes
    sieve_primes_log(n, factor_base, factor_base_size, S, sieve);

    // Create the matrix
    GArray* matrix = g_array_new(FALSE, FALSE, sizeof(bool*));
    GArray* as_vector = g_array_new(FALSE, FALSE, sizeof(mpz_t*));
    GHashTable* factor_exponent_dict = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    create_matrix_log(sieve, S, root_n, factor_base, factor_base_size, n, matrix, as_vector, factor_exponent_dict);

    // Find linear dependencies
    GArray* dependencies = g_array_new(FALSE, FALSE, sizeof(GArray*));
    find_linear_dependencies(dependencies, matrix, factor_base_size);

    printf("Number of dependencies: %d with matrix shape: %d x %d\n", dependencies->len, matrix->len, factor_base_size);

    // return factors
    return_factors(factor1, factor2, dependencies, as_vector, factor_exponent_dict, factor_base, factor_base_size, n);

    // Free memory
    for (int i = 0; i < matrix->len; i++) {
        free(g_array_index(matrix, bool*, i));
    }
    g_array_free(matrix, TRUE);

    for (int i = 0; i < as_vector->len; i++) {
        mpz_t* value = g_array_index(as_vector, mpz_t*, i);
        mpz_clear(*value);
        g_free(value);
    }
    g_array_free(as_vector, TRUE);

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


int main(int argc, char* argv[]) {
    // Initialization
    mpz_t n, factor1, factor2;
    int B, S;
    mpz_init(n);
    mpz_init(factor1);
    mpz_init(factor2);

    // Commandline argument parsing: ./quadratic_sieve <number to factor> <B> <S>
    if(argc != 4) {
        puts("Usage: ./quadratic_sieve <number to factor> <B> <S>");
        exit(1);
    } else {
        if(mpz_set_str(n, argv[1], 10) != 0) {
            puts("ERROR: Unable to parse number.");
            exit(1);
        }
        B = atoi(argv[2]);
        S = atoi(argv[3]);
        if (B <= 0 || S <= 0) {
            puts("ERROR: B and S must be positive integers.");
            exit(1);
        }
    }

    gmp_printf("Starting sieve with n: %Zd, B: %d, S: %d\n", n, B, S);

    // Run and time the quadratic sieve
    clock_t start, end;
    start = clock();

    sieve_log(n, B, S, factor1, factor2);

    end = clock();
    double time_taken = ((double)end - start) / CLOCKS_PER_SEC * 1000;

    printf("Time taken: %f ms\n", time_taken);
    gmp_printf("%Zd has factors: (%Zd, %Zd)\n", n, factor1, factor2);

    // Cleanup
    mpz_clear(n);
    mpz_clear(factor1);
    mpz_clear(factor2);

    return 0;
}
