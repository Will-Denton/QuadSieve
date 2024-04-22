import math
import numpy as np
from tqdm import tqdm
import sys

# Steps for the quadratic sieve:
# 1. Choose a smoothness bound B. The number π(B), denoting the number of prime numbers less than B, will control both the length of the vectors and the number of vectors needed.
# 2. Use sieving to locate π(B) + 1 numbers ai such that bi = (ai2 mod n) is B-smooth.
# 3. Factor the bi and generate exponent vectors mod 2 for each one.
# 4. Use linear algebra to find a subset of these vectors which add to the zero vector. Multiply the corresponding ai together and give the result mod n the name a; similarly, multiply the bi together which yields a B-smooth square b2.
# 5. We are now left with the equality a2 = b2 mod n from which we get two square roots of (a2 mod n), one by taking the square root in the integers of b2 namely b, and the other the a computed in step 4.
# 6. We now have the desired identity. Compute the GCD of n with the difference (or sum) of a and b. This produces a factor, although it may be a trivial factor (n or 1). If the factor is trivial, try again with a different linear dependency or different a.


def sieve_of_eratosthenes(B):
    is_prime = np.ones(B + 1, dtype=bool)
    is_prime[:2] = False
    for i in tqdm(range(2, int(B**0.5) + 1), desc="Sieve of Eratosthenes: "):
        if is_prime[i]:
            is_prime[i * i : B + 1 : i] = False
    return np.nonzero(is_prime)[0]


def quadratic_residue(n, p):
    assert isinstance(n, int) and isinstance(p, int)
    if p == 2:
        return 1

    return pow(n, (p - 1) // 2, p)


def get_factor_base(primes, n):
    # in order for shanks-tonelli to work, the factor base must contain p such that n is a quadratic residue mod p
    # limit the primes found from the sieve of eratosthenes to only those that are quadratic residues
    factor_base = []
    for p in tqdm(primes, desc="Creating Factor Base: "):
        if quadratic_residue(n, int(p)) == 1:
            factor_base.append(p)
    return np.array(factor_base, dtype=int)


def get_sieve(S, n):
    # create sieve of size S
    sieve = []
    root_n = math.ceil(math.sqrt(n))
    for i in tqdm(range(S), desc="Creating Sieve: "):
        a = i + root_n
        b = a * a - n
        sieve.append(b)
    return sieve


def get_sieve_log(S, n):
    print("Creating Sieve Log ...")
    # everythig is float64 since numpy list compression cant do python bigInts with size > 64 bits
    # this causes a lot of error to accumulate, so epsilon is used to offset this later
    root_n = np.float64(np.ceil(math.sqrt(n)))
    i_values = np.arange(0, S, dtype=np.float64)
    sieve = np.log((i_values + root_n) ** 2 - np.float64(n))
    return sieve


def shanks_tonelli(n, p):
    # calculate the square roots of n mod p
    # DID THIS https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    if p == 2:
        return [1]

    # 1. find S and Q such that p - 1 = Q * 2^S
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q = Q // 2
        S += 1

    # 2. find a quadratic non-residue z
    z = 2
    while quadratic_residue(z, p) != p - 1:
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
        b = pow(c, 2 ** (M - i - 1), p)
        M = i
        c = pow(b, 2, p)
        t = (t * pow(b, 2, p)) % p
        R = (R * b) % p


def sieve_primes(n, factor_base, S, sieve):
    # sieve the numbers in the sieve that are divisible by primes in the factor base
    # shanks-tonelli is used to find the square roots of n mod p and solutions to this
    # are used to find the numbers in the sieve that are divisible by p in the sieve
    # super ineffiecient for small values of p
    root_n = math.ceil(math.sqrt(n))
    for p in tqdm(factor_base, desc="Preforming Shanks-Tonelli: "):
        roots_mod_p = shanks_tonelli(n, int(p))
        for root in roots_mod_p:
            x = (root - root_n) % p
            for i in range(x, S, p):
                sieve[i] = sieve[i] // p


def sieve_primes_log(n, factor_base, S, log_sieve):
    root_n = math.ceil(math.sqrt(n))
    for p in tqdm(factor_base, desc="Preforming Sieve-Primes Numpy: "):
        roots_mod_p = shanks_tonelli(n, int(p))
        for root in roots_mod_p:
            x = (root - root_n) % p
            log_sieve[x:S:p] -= np.log(p)


def get_B_smooth_factors(b, factor_base):
    # trial division to find the factors of b
    factors = []
    for p in factor_base:
        while b % int(p) == 0:
            b = b // int(p)
            factors.append(int(p))
    return np.array(factors, dtype=int)


def get_factor_vector(factors, factor_base):
    # create a vector of the exponents of the factors in the factor base
    exponent_vector = np.zeros(len(factor_base), dtype=int)
    lookup_factor_index = {factor: i for i, factor in enumerate(factor_base)}
    for factor in factors:
        exponent_vector[lookup_factor_index[factor]] += 1
    return exponent_vector


def create_matrix(sieve, root_n, factor_base, n):
    # create the matrix mod 2 for finding dependencies
    # as vector is vector of a values corresponding to the b values in the sieve
    # factor_exponent_dict is a dictionary of the exponent vectors for each a value
    epsilon = 1e-2  # this offsets the ammout of error accumulated by switching to float64 earlier
    matrix = []
    as_vector = []
    factor_exponent_dict = {}
    for i in tqdm(range(len(sieve)), desc="Creating Matrix: "):
        if sieve[i] < epsilon:
            factors = get_B_smooth_factors((i + root_n) ** 2 - n, factor_base)
            exponent_vector = get_factor_vector(factors, factor_base)
            if np.sum(exponent_vector % 2) == 0:
                continue
            matrix.append((exponent_vector % 2).tolist())
            as_vector.append(i + root_n)
            factor_exponent_dict[i + root_n] = exponent_vector
            if len(matrix) >= 2 * len(
                factor_base
            ):  # if rows is double the columns, break
                break
    matrix = np.array(matrix, dtype=bool)
    return matrix, as_vector, factor_exponent_dict


def find_linear_dependencies(matrix):
    # find the linear dependencies in the matrix
    # follows this paper https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    n, m = matrix.shape
    marks = np.zeros(n, dtype=bool)
    for i in range(m):
        col = matrix[:, i]
        piv = -1
        for j in range(n):
            if col[j] == 1:  # matrix[j][i] == 1
                marks[j] = True
                piv = j
                break
        if piv != -1:
            for k in range(m):
                if k != i:
                    col2 = matrix[:, k]
                    if col2[piv] == 1:  # matrix[piv][k] == 1
                        matrix[:, k] = matrix[:, k] ^ col  # matrix[:, k] = matrix[:, k] ^ matrix[:, i]

    dependencies = []
    for i in range(n):
        if marks[i] == False:
            dependent_list = []
            dependent_list.append(i)
            dependent_row = matrix[i]
            for j in range(m):
                if dependent_row[j] == 1:  # matrix[i][j] == 1
                    dependent_column = matrix[:, j]
                    for k in range(n):
                        if dependent_column[k] == 1:  # matrix[k][j] == 1
                            dependent_list.append(k)
                            break
            dependencies.append(dependent_list)
    return dependencies


def calculate_as_product(depenencies, exponent_as):
    as_product = 1
    for row in depenencies:
        as_product *= exponent_as[row]
    return as_product


def calculate_primes_product(depenencies, factor_exponent_dict, as_vector, factor_base):
    primes_product = 1
    prime_power_vector = np.zeros(len(factor_base), dtype=int)
    for row in depenencies:
        prime_power_vector = prime_power_vector + factor_exponent_dict[as_vector[row]]
    prime_power_vector = prime_power_vector // 2

    for i, p in enumerate(factor_base):
        primes_product *= pow(int(p), int(prime_power_vector[i]))

    return primes_product


def euclidian_algorithm(a, b):
    while b != 0:
        a, b = b, a % b
    return a


def return_factors(dependencies, as_vector, factor_exponent_dict, factor_base, n):
    if len(dependencies) == 0:
        return None

    for dependency in dependencies:
        as_product = calculate_as_product(dependency, as_vector)
        primes_product = calculate_primes_product(
            dependency, factor_exponent_dict, as_vector, factor_base
        )
        f = euclidian_algorithm(primes_product - as_product, n)
        if f != 1 and f != n:
            return f, n // f

    return None


def sieve(n, B, S):
    root_n = math.ceil(math.sqrt(n))

    primes_under_B = sieve_of_eratosthenes(B)
    factor_base = get_factor_base(primes_under_B, n)
    sieve = get_sieve_log(S, n)

    sieve_primes_log(n, factor_base, S, sieve)

    print("Finished shanks tonelli, starting matrix creation...")
    matrix, as_vector, factor_exponent_dict = create_matrix(
        sieve, root_n, factor_base, n
    )

    print(
        f"Finished matrix creation with size: {matrix.shape}, finding linear dependencies..."
    )
    matrix_rr = np.copy(matrix)
    dependencies = find_linear_dependencies(matrix_rr)

    # save matrix to file
    # np.save("matrix.npy", matrix)
    with open("matrix.txt", "w") as f:
        for row in matrix:
            f.write("[")
            for r in row:
                if r:
                    f.write("1 ")
                else:
                    f.write("0 ")
                # delete the last space
            f.seek(f.tell() - 1)
            f.write("]")
            f.write("\n")

    # save dependencies to file - dependencies have different sizes so cant save as numpy array
    with open("dependencies_easy.txt", "w") as f:
        for d in dependencies:
            f.write(str(d) + "\n")

    #save as_vector to file
    # as_vector_np = np.array(as_vector)
    # np.save("as_vector.npy", as_vector_np)

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

        f = factor_exponent_dict[a]
        primes_product = 1
        for j, p in enumerate(factor_base):
            primes_product *= pow(int(p), int(f[j]))

        factors = []
        for i in range(len(f)):
            if f[i] != 0:
                factors.append(factor_base[i] ** f[i])
        if primes_product != calc_b:
            print(get_B_smooth_factors(b, factor_base))
            print(f"factors: {factors}")
            print(f"primes_product: {primes_product}")
            print(f"b: {calc_b}")
        assert primes_product == calc_b



    print(
        f"Finished finding {len(dependencies)} linear dependencies, looking for factors..."
    )
    return return_factors(dependencies, as_vector, factor_exponent_dict, factor_base, n)


if __name__ == "__main__":
    # n, B, S = 16921456439215439701, 2000, 4000000
    # n, B, S = 46839566299936919234246726809, 15000, 15000000
    n, B, S = 6172835808641975203638304919691358469663, 30000, 1000000000
    print(f"n: {n}, factors: {sieve(n, B, S)}")

    # n2 = 16921456439215439701
    # n, B, S = 46839566299936919234246726809, 15000, 2
    # print(n2.bit_length())
    # get_sieve_log(S, n2)

    # print number of bits in n

    # primes_under_B = sieve_of_eratosthenes(B)
    # factor_base = get_factor_base(primes_under_B, n)
    # sieve_s = get_sieve_numpy(S, n)
    # sieve_primes_numpy(n, factor_base, S, sieve_s)
    # epsilon = 1e-6
    # #print index where sieve_s is less than epsilon
    # sieve_numpy_ans = np.where(sieve_s < epsilon)
    # sieve_numpy_ans = sieve_numpy_ans[0]

    # sieve_n = get_sieve(S, n)
    # sieve_primes(n, factor_base, S, sieve_n)
    # #make empty numpy array of length 0
    # sieve_ans = np.zeros(0, dtype=int)
    # for i in range(S):
    #     if sieve_n[i] == 1:
    #         sieve_ans = np.append(sieve_ans, i)
    # print(sieve_ans)

    # #check to see if sieve_ans and sieve_numpy_ans are the same
    # print(np.array_equal(sieve_ans, sieve_numpy_ans))
    # print(len(sieve_ans), len(sieve_numpy_ans))
