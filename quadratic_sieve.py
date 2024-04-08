import math
import numpy as np

# Steps for the quadratic sieve:
# 1. Choose a smoothness bound B. The number π(B), denoting the number of prime numbers less than B, will control both the length of the vectors and the number of vectors needed.
# 2. Use sieving to locate π(B) + 1 numbers ai such that bi = (ai2 mod n) is B-smooth.
# 3. Factor the bi and generate exponent vectors mod 2 for each one.
# 4. Use linear algebra to find a subset of these vectors which add to the zero vector. Multiply the corresponding ai together and give the result mod n the name a; similarly, multiply the bi together which yields a B-smooth square b2.
# 5. We are now left with the equality a2 = b2 mod n from which we get two square roots of (a2 mod n), one by taking the square root in the integers of b2 namely b, and the other the a computed in step 4.
# 6. We now have the desired identity. Compute the GCD of n with the difference (or sum) of a and b. This produces a factor, although it may be a trivial factor (n or 1). If the factor is trivial, try again with a different linear dependency or different a.


DO_ASSERTS = True


def sieve_of_eratosthenes(B):
    # calculate primes up to B
    # https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
    factor_base = np.array([], dtype=int)
    is_primes = np.ones(B + 1, dtype=bool)
    is_primes[0] = False
    is_primes[1] = False
    for i in range(2, math.ceil(math.sqrt(B + 1))):
        if is_primes[i]:
            factor_base = np.append(factor_base, i)
            for j in range(i * i, B + 1, i):
                is_primes[j] = False
    for i in range(math.ceil(math.sqrt(B + 1)), B + 1):
        if is_primes[i]:
            factor_base = np.append(factor_base, i)
    return factor_base


def quadratic_residue(n, p):
    # check if n is a quadratic residue mod p (n has a square root mod p)
    # https://en.wikipedia.org/wiki/Quadratic_residue
    assert isinstance(n, int) and isinstance(p, int)
    if p == 2:
        return 1

    n = n % p
    r = (p-1) // 2

    a = 1
    while r != 0:
        if r % 2 == 0:
            r = r // 2
            n = n*n % p
        else:
            r = r - 1
            a = a*n % p
    return a


def get_factor_base(primes, n):
    # in order for shanks-tonelli to work, the factor base must contain p such that n is a quadratic residue mod p
    # limit the primes found from the sieve of eratosthenes to only those that are quadratic residues
    factor_base = np.array([], dtype=int)
    for p in primes:
        if quadratic_residue(n, int(p)) == 1:
            factor_base = np.append(factor_base, p)
    return factor_base


def get_sieve(S, n):
    # create sieve of size S
    sieve = []
    root_n = math.ceil(math.sqrt(n))
    for i in range(S):
        a = i + root_n
        b = a*a - n
        sieve.append(b)
    return sieve


def shanks_tonelli(n, p):
    # calculate the square roots of n mod p
    #DID THIS https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    if p == 2:
        return [1]
    
    # 1. find S and Q such that p - 1 = Q * 2^S
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q = Q // 2
        S += 1

    if S == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    
    # 2. find a quadratic non-residue z
    z = 2
    while quadratic_residue(z, p) != p-1:
        z += 1

    # 3. initialize M, c, t, R
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)

    # 4. loop
    while True:
        if t == 1:
            return (R, p - R)
        i = 1
        while pow(t, 2**i, p) != 1:
            i += 1
        b = pow(c, 2**(M-i-1), p)
        M = i
        c = pow(b, 2, p)
        t = (t * pow(b, 2, p)) % p
        R = (R * b) % p


def sieve_primes(n, factor_base, S, sieve):
    # sieve the numbers in the sieve that are divisible by primes in the factor base
    # shanks-tonelli is used to find the square roots of n mod p and solutions to this
    # are used to find the numbers in the sieve that are divisible by p in the sieve
    root_n = math.ceil(math.sqrt(n))
    for p in factor_base:
        roots_mod_p = shanks_tonelli(n, int(p))
        for root in roots_mod_p:
            x = (root - root_n) % p
            for i in range(x, S, p):
                sieve[i] = sieve[i] // p


def trial_division(b, factor_base):
    # trial division to find the factors of b
    factors = np.array([], dtype=int)
    for p in factor_base:
        p = int(p)
        while b % p == 0:
            b = b // p
            factors = np.append(factors, p)
    return factors


def get_B_smooth_factors(b, factor_base):
    return trial_division(b, factor_base)


def get_factor_vector(factors, factor_base):
    # create a vector of the exponents of the factors in the factor base
    exponent_vector = np.zeros(len(factor_base), dtype=int)
    lookup_factor_index = {factor: i for i, factor in enumerate(factor_base)}
    for factor in factors:
        exponent_vector[lookup_factor_index[factor]] += 1
    return exponent_vector


def create_matrix(sieve, original_sieve, root_n, factor_base):
    # create the matrix mod 2 for finding dependencies
    # as vector is vector of a values corresponding to the b values in the sieve
    # factor_exponent_dict is a dictionary of the exponent vectors for each a value
    matrix = []
    as_vector = []
    bs_vector = []
    factor_exponent_dict = {}
    for i in range(len(sieve)):
        if sieve[i] == 1:
            factors = get_B_smooth_factors(original_sieve[i], factor_base)
            exponent_vector = get_factor_vector(factors, factor_base)
            if np.sum(exponent_vector % 2) == 0:
                continue
            matrix.append((exponent_vector % 2).tolist())
            as_vector.append(i + root_n)
            bs_vector.append(original_sieve[i])
            factor_exponent_dict[i + root_n] = exponent_vector
            if len(matrix) >= 2 * len(factor_base): # if rows is double the columns, break
                break
    matrix = np.array(matrix, dtype=bool)
    return matrix, as_vector, bs_vector, factor_exponent_dict


def find_linear_dependencies(matrix):
    # find the linear dependencies in the matrix
    # follows this paper https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    n,m = matrix.shape
    marks = np.zeros(n, dtype=bool)
    for i in range(m):
        col = matrix[:, i]
        piv = -1
        for j in range(n):
            if col[j] == 1:
                marks[j] = True
                piv = j
                break
        if piv != -1:
            for k in range(m):
                if k != i:
                    col2 = matrix[:, k]
                    if col2[piv] == 1:
                        matrix[:, k] = (matrix[:, k] ^ col) 

    dependencies = []
    for i in range(n):
        if marks[i] == False:
            dependent_list = []
            dependent_list.append(i)
            dependent_row = matrix[i]
            for j in range(m):
                if dependent_row[j] == 1:
                    dependent_column = matrix[:, j]
                    for k in range(n):
                        if dependent_column[k] == 1:
                            dependent_list.append(k)
                            break
            dependencies.append(dependent_list)
    return dependencies


def calculate_as_product(depenencies, exponent_as):
    as_product = 1
    for row in depenencies:
        as_product *= exponent_as[row]
        as_product = as_product
    return as_product


def calculate_primes_product(depenencies, factor_exponent_dict, as_vector, factor_base):
    primes_product = 1
    prime_power_vector = np.zeros(len(factor_base), dtype=int)
    for row in depenencies:
        prime_power_vector = prime_power_vector + factor_exponent_dict[as_vector[row]]
    prime_power_vector = prime_power_vector // 2

    for i, p in enumerate(factor_base):
        primes_product *= pow(int(p), int(prime_power_vector[i]))
        primes_product = primes_product

    return primes_product


def euclidian_algorithm(a, b):
    while b != 0:
        a, b = b, a % b
    return a

def sieve(n, B, S):
    root_n = math.ceil(math.sqrt(n))

    print("Creating sieve...")
    primes_under_B = sieve_of_eratosthenes(B)
    factor_base = get_factor_base(primes_under_B, n)
    sieve = get_sieve(S, n)
    original_sieve = np.copy(sieve)

    print("Created sieve, starting shanks tonelli...")
    sieve_primes(n, factor_base, S, sieve)

    print("Finished shanks tonelli, starting matrix creation...")
    matrix, as_vector, bs_vector, factor_exponent_dict = create_matrix(sieve, original_sieve, root_n, factor_base)
    print(f"matrix: {matrix.shape}")

    print("Finished matrix creation, finding linear dependencies...")
    matrix_rr = np.copy(matrix)
    dependencies = find_linear_dependencies(matrix_rr)
        
    # assert that all dependencies are actually dependencies
    if DO_ASSERTS:
        for d in dependencies:
            total = np.zeros(len(factor_base), dtype=bool)  
            for i in d:
                row = matrix[i]
                total = total ^ row
            assert np.array_equal(total, np.zeros(len(factor_base), dtype=bool))

        # assert prime factors actually multiply to a
        for i in range(matrix.shape[0]):
            a = as_vector[i]
            calc_b = a*a - n
            b = bs_vector[i]
            assert calc_b == b

            f = factor_exponent_dict[a]
            primes_product = 1
            for j, p in enumerate(factor_base):
                primes_product *= pow(int(p), int(f[j]))

            factors = []
            for i in range(len(f)):
                if f[i] != 0:
                    factors.append(factor_base[i] ** f[i])
            if primes_product != b:
                print(get_B_smooth_factors(b, factor_base))
                print(f"factors: {factors}")
                print(f"primes_product: {primes_product}")
                print(f"b: {b}")
            assert primes_product == b

    print("Finished finding linear dependencies, looking for factors...")
    print(f"number of dependencies: {len(dependencies)}")
    for depedency in dependencies:
        as_product = calculate_as_product(depedency, as_vector)
        primes_product = calculate_primes_product(depedency, factor_exponent_dict, as_vector, factor_base)
        f = euclidian_algorithm(primes_product - as_product, n)
        
        if f != 1 and f != n:
            return f, n // f
        
    return None

if __name__ == "__main__":
    n, B, S = 46839566299936919234246726809, 15000, 15000000
    # n, B, S = 16921456439215439701, 2000, 4000000
    # n, B, S = 6172835808641975203638304919691358469663, 15000, 25000000
    print(f"n: {n}, factors: {sieve(n, B, S)}")
    # primes_under_B = sieve_of_eratosthenes(B)
    # factor_base = get_factor_base(primes_under_B, n)
    # print(get_B_smooth_factors(34539161229042973095,factor_base))
