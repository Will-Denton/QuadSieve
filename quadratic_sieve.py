import math
import numpy as np
from tqdm import tqdm


def sieve_of_eratosthenes(B):
    """
    args: (int) B - upper bound for primes
    returns: numpy array of all primes < B

    Runs the sieve of eratosthenes algorithm to find all primes < B.
    
    """
    is_prime = np.ones(B + 1, dtype=bool)                                       # boolean array from 0 to B
    is_prime[:2] = False
    for i in tqdm(range(2, int(B**0.5) + 1), desc="Sieve of Eratosthenes: "):
        if is_prime[i]:                                                         # if i is prime
            is_prime[i * i : B + 1 : i] = False                                 # set all multiples of i to false
    return np.nonzero(is_prime)[0]                                              # return all indices where is_prime is true


def quadratic_residue(n, p):
    """
    args: (int) n, (int) p - n is the number to check if it is a quadratic residue mod p
    returns: 1 if n is a quadratic residue mod p, 0 otherwise

    Calculates the legendre symbol (n/p) = n^((p-1)/2) mod p.
    """
    assert isinstance(n, int) and isinstance(p, int)                            # check if n and p are integers
    if p == 2:
        return 1                                                                # all numbers are quadratic residues mod 2
    return pow(n, (p - 1) // 2, p)                                              # return n^((p-1)/2) mod p


def get_factor_base(primes, n):
    """
    args: (np array) primes - numpy array of prime numbers
          (int) n - the number to factor
    returns: numpy array of primes that n is a quadratic residue mod p

    Creates the factor base for the quadratic sieve.
    In order for shanks-tonelli to work, the factor base must contain p such that n is a quadratic residue mod p.
    Limit the primes found from the sieve of eratosthenes to only those that are quadratic residues.
    """
    factor_base = []
    for p in tqdm(primes, desc="Creating Factor Base: "):
        if quadratic_residue(n, int(p)) == 1:                                   # if n is a quadratic residue mod p
            factor_base.append(p)                                               # add p to the factor base                                   
    return np.array(factor_base, dtype=int)                                     # return the factor base as a numpy array


def get_sieve_log(S, n):
    """
    args: (int) S - the size of the sieve
          (int) n - the number to factor
    returns: log of the list of numbers to sieve
    
    Given a bound S, create a list of numbers to sieve.
    The sieve list is given by ln((i + root_n)^2 - n) for i in range(S).
    All floats are 64 bit since numpy list compression cant use anything larger.
    This causes a lot of error to accumulate, so epsilon is used to offset this later.
    """
    print("Creating Sieve Log ...")
    root_n = np.float64(np.ceil(math.sqrt(n)))
    i_values = np.arange(0, S, dtype=np.float64)                                # create a numpy array from 0 to S
    sieve = np.log((i_values + root_n) ** 2 - np.float64(n))                    # for all i in the array, calculate ln((i + root_n)^2 - n)
    return sieve


def shanks_tonelli(n, p):
    """
    args: (int) n, (int) p - n is the number to find the square root of mod p
    returns: list of square roots of n mod p (there are at most 2 square roots mod p)

    Standard shanks-tonelli algorithm to find the square roots of n mod p.
    """
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


def sieve_primes_log(n, factor_base, S, log_sieve):
    """
    args: (int) n - the number to factor
          (np array) factor_base - the factor base for the quadratic sieve
          (int) S - the size of the sieve
          (np array) sieve - the list of numbers to sieve
    returns: None
    
    Sieve the numbers in the sieve that are divisible by primes in the factor base.
    Shanks-tonelli is used to find the square roots of n mod p and solutions to this
    are used to find the numbers in the sieve that are divisible by p in the sieve.
    All values in the sieve that are about 0 after this process are B-smooth.
    """
    root_n = math.ceil(math.sqrt(n))
    for p in tqdm(factor_base, desc="Preforming Sieve-Primes Numpy: "):
        roots_mod_p = shanks_tonelli(n, int(p))                                 # find the square roots of n mod p
        for root in roots_mod_p:
            x = (root - root_n) % p                                             # solve for x in (x + root_n)^2 - n
            log_sieve[x:S:p] -= np.log(p)                                       # for all x in the sieve, subtract ln(p) where x is divisible by p


def get_B_smooth_factors(b, factor_base):
    """
    args: (int) b - the number to factor
          (np array) factor_base - the factor base for the quadratic sieve
    returns: numpy array of factors of b that are in the factor base

    The sieve process does not find the factors of b, only that b is B-smooth.
    This function finds the factors of b that are in the factor base using trial division.
    """
    factors = []
    for p in factor_base:
        while b % int(p) == 0:                                                  # while b is divisible by p
            b = b // int(p)                                                     # divide b by p                      
            factors.append(int(p))                                              # add p to the list of factors
    return np.array(factors, dtype=int)


def get_factor_vector(factors, factor_base):
    """
    args: (np array) factors - the factors of b that are in the factor base
          (np array) factor_base - the factor base for the quadratic sieve
    returns: numpy array of the exponents of the factors in the factor base

    Creates the exponent vector for the factors of b in the factor base.
    The exponent vector is given as the number of times each factor in the factor base divides b.
    """
    exponent_vector = np.zeros(len(factor_base), dtype=int)
    lookup_factor_index = {factor: i for i, factor in enumerate(factor_base)}   # create a dictionary to lookup the index of a factor
    for factor in factors:
        exponent_vector[lookup_factor_index[factor]] += 1                       # increment the exponent of the factor in the exponent vector  
    return exponent_vector


def create_matrix(sieve, root_n, factor_base, n):
    """
    args: (np array) sieve - the list of numbers to sieve
          (int) root_n - the square root of n
          (np array) factor_base - the factor base for the quadratic sieve
          (int) n - the number to factor
    returns: matrix - the matrix mod 2 for finding dependencies
             as_vector - the list of a values corresponding to the b values in the sieve
             factor_exponent_dict - a dictionary of the exponent vectors for each a value

    Creates the matrix mod 2 for finding dependencies in the sieve.
    The matrix is created by finding the B-smooth numbers in the sieve and creating the exponent vectors for each a value.
    The goal is to have a matrix with more rows than columns to find linear dependencies.
    """
    epsilon = 1e-2                                                              # this offsets the ammout of error accumulated by switching to float64 earlier
    matrix = []
    as_vector = []
    factor_exponent_dict = {}
    for i in tqdm(range(len(sieve)), desc="Creating Matrix: "):
        if sieve[i] < epsilon:                                                  # if the sieve value is less than epsilon (around 0)                       
            factors = get_B_smooth_factors((i + root_n) ** 2 - n, factor_base)  # find the factors of (i + root_n)^2 - n that are in the factor base
            exponent_vector = get_factor_vector(factors, factor_base)           # create the exponent vector for the factors
            if np.sum(exponent_vector % 2) == 0:                                # if all exponents are even, then the B-smooth number is not needed
                continue
            matrix.append((exponent_vector % 2).tolist())
            as_vector.append(i + root_n)
            factor_exponent_dict[i + root_n] = exponent_vector
            if len(matrix) >= 2 * len(factor_base):                             # if rows is double the columns, break
                break
    matrix = np.array(matrix, dtype=bool)
    return matrix, as_vector, factor_exponent_dict


def find_linear_dependencies(matrix):
    """
    args: (np array) matrix - the matrix mod 2 for finding dependencies
    returns: list of linear dependencies in the matrix

    Finds the linear dependencies in the matrix using gaussian elimination mod 2.
    Follows the algorithm presented in this paper: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    
    """
    # algorithm presented in the paper:
    # 1. Search for Aij = 1 in column j
    # 2. If found then mark row i as an independent row
    # 3. For all other columns k != j, if Akj = 1, then xor column k with column j
    # 4. repeat for all columns
    n, m = matrix.shape
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
                        matrix[:, k] = matrix[:, k] ^ col

    # All unmarked rows are dependent on the marked rows
    # loop through all unmarked rows and create a list of rows that the unmarked row is dependent on
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
    """
    args: (list) dependencies - list of dependencies
          (np array) exponent_as - the list of a values corresponding to the b values in the sieve
    returns: the product of the a values in the dependencies

    Calculates the product of the as for (a1...an)^2 = (p1^e1...pn^en)^2 mod n.
    """
    as_product = 1
    for row in depenencies:
        as_product *= exponent_as[row]                                          # multiply all a values together     
    return as_product


def calculate_primes_product(depenencies, factor_exponent_dict, as_vector, factor_base):
    """
    args: (list) dependencies - list of dependencies
          (dict) factor_exponent_dict - a dictionary of the exponent vectors for each a value
          (np array) as_vector - the list of a values corresponding to the b values in the sieve
          (np array) factor_base - the factor base for the quadratic sieve
    returns: the product of the primes in the dependencies

    Calculates the product of the primes for (a1...an)^2 = (p1^e1...pn^en)^2 mod n.
    """
    primes_product = 1
    prime_power_vector = np.zeros(len(factor_base), dtype=int)               # create a vector to store the exponents of the primes
    for row in depenencies:
        prime_power_vector = prime_power_vector + factor_exponent_dict[as_vector[row]] 
    prime_power_vector = prime_power_vector // 2                             # divide the exponents by 2 (taking the square root)           

    for i, p in enumerate(factor_base):
        primes_product *= pow(int(p), int(prime_power_vector[i]))            # multiply all primes to their powers together

    return primes_product


def euclidian_algorithm(a, b):
    """
    args: (int) a, (int) b - the numbers to find the gcd of
    returns: the gcd of a and b

    Finds the greatest common divisor of a and b using the euclidian algorithm.
    """
    while b != 0:
        a, b = b, a % b
    return a


def return_factors(dependencies, as_vector, factor_exponent_dict, factor_base, n):
    """
    args: (list) dependencies - list of dependencies
          (np array) as_vector - the list of a values corresponding to the b values in the sieve
          (dict) factor_exponent_dict - a dictionary of the exponent vectors for each a value
          (np array) factor_base - the factor base for the quadratic sieve
          (int) n - the number to factor
    returns: the factors of n

    Returns the factors of n by calculating the product of the as and the product of the primes.
    Uses the basic principle that given (a1...an)^2 = (p1^e1...pn^en)^2 mod n, then a factor is = gcd(a1...an - p1^e1...pn^en, n).
    """
    if len(dependencies) == 0:                                              # if there are no dependencies return
        return None

    for dependency in dependencies:
        as_product = calculate_as_product(dependency, as_vector)
        primes_product = calculate_primes_product(dependency, factor_exponent_dict, as_vector, factor_base)
        f = euclidian_algorithm(primes_product - as_product, n)             # find the gcd of the product of the as and the product of the primes
        if f != 1 and f != n:                                               # if the gcd is non trivial, return the factors
            return f, n // f

    return None


def sieve(n, B, S):
    """
    args: (int) n - the number to factor
          (int) B - the upper bound for primes
          (int) S - the size of the sieve
    returns: the factors of n

    Main function to factor n using the quadratic sieve.
    """

    # initialize variables
    root_n = math.ceil(math.sqrt(n))
    primes_under_B = sieve_of_eratosthenes(B)
    factor_base = get_factor_base(primes_under_B, n)
    sieve = get_sieve_log(S, n)

    # sieve the numbers in the sieve that are divisible by primes in the factor base
    sieve_primes_log(n, factor_base, S, sieve)

    # find the B-smooth numbers in the sieve and create the matrix mod 2
    print("Finished shanks tonelli, starting matrix creation...")
    matrix, as_vector, factor_exponent_dict = create_matrix(sieve, root_n, factor_base, n)

    # find linear dependencies in the matrix
    print(f"Finished matrix creation with size: {matrix.shape}, finding linear dependencies...")
    matrix_rr = np.copy(matrix)
    dependencies = find_linear_dependencies(matrix_rr)

    # return the factors of n
    return return_factors(dependencies, as_vector, factor_exponent_dict, factor_base, n)


if __name__ == "__main__":
    # uncomment the following lines to test the quadratic sieve

    n, B, S = 16921456439215439701, 2000, 4000000
    # n, B, S = 46839566299936919234246726809, 15000, 15000000
    # n, B, S = 6172835808641975203638304919691358469663, 30000, 1000000000
    print(f"n: {n}, factors: {sieve(n, B, S)}")