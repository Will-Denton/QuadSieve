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