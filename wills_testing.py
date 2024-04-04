import math

# import msvcrt
import numpy as np
from sympy import Matrix

# Steps for the quadratic sieve:
# 1. Choose a smoothness bound B. The number π(B), denoting the number of prime numbers less than B, will control both the length of the vectors and the number of vectors needed.
# 2. Use sieving to locate π(B) + 1 numbers ai such that bi = (ai2 mod n) is B-smooth.
# 3. Factor the bi and generate exponent vectors mod 2 for each one.
# 4. Use linear algebra to find a subset of these vectors which add to the zero vector. Multiply the corresponding ai together and give the result mod n the name a; similarly, multiply the bi together which yields a B-smooth square b2.
# 5. We are now left with the equality a2 = b2 mod n from which we get two square roots of (a2 mod n), one by taking the square root in the integers of b2 namely b, and the other the a computed in step 4.
# 6. We now have the desired identity. Compute the GCD of n with the difference (or sum) of a and b. This produces a factor, although it may be a trivial factor (n or 1). If the factor is trivial, try again with a different linear dependency or different a.


NUMBER_OF_VECTORS_IN_SPACE_MOD = 1


def get_smoothness_bound(n):
    # return math.floor(math.exp(math.sqrt(math.log(n) * math.log(math.log(n))))
    return math.ceil(math.sqrt(n))


def incriment_x(x):
    # if x > 0:
    #     return -x
    # else:
    #     return -x + 1
    return x + 1


def sieve_of_eratosthenes(B):
    factor_base = []
    is_primes = [True] * (B + 1)
    for i in range(2, math.ceil(math.sqrt(B + 1))):
        if is_primes[i]:
            factor_base.append(i)
            for j in range(i * i, B + 1, i):
                is_primes[j] = False
    for i in range(math.ceil(math.sqrt(B + 1)), B + 1):
        if is_primes[i]:
            factor_base.append(i)
    return factor_base


def is_B_smooth_trial_division(factor_base, b):
    factors = []
    if b == 0:
        return None
    for p in factor_base:
        # print(f"p: {p}, b1: {b1}")
        while b / p == b // p:
            b = b // p
            factors.append(p)
    if b != 1 or len(factors) == 0:
        return None
    return factors


def calculate_a_product(ns_vector, exponent_as):
    as_product = 1
    for i, a in enumerate(ns_vector):
        if a == 1:
            as_product *= exponent_as[i]
            as_product = as_product % n
    return as_product


# replace with pringalas algorithm
def prime_to_power(p, power, n):
    prime_power = 1
    for i in range(power):
        prime_power *= p
        prime_power = prime_power % n
    return prime_power


def calculate_primes_product(ns_vector, exponent_vectors_actual, factor_base):
    primes_product = 1
    prime_power_vector = np.zeros(len(factor_base), dtype=int)
    for i, a in enumerate(ns_vector):
        if a == 1:
            prime_power_vector += exponent_vectors_actual[i]

    prime_power_vector = prime_power_vector // 2

    print(f"prime_power_vector: {prime_power_vector}")

    for i, p in enumerate(factor_base):
        primes_product *= prime_to_power(p, prime_power_vector[i], n)
        primes_product = primes_product % n
    return primes_product


def sieve(n):
    # Choose smoothness bound B
    # online found that B is approximately equal to exp(sqrt(log(n) * log(log(n))))
    B = get_smoothness_bound(n)

    # Now we have B compute primes up to B using Sieve of Eratosthenes
    factor_base = sieve_of_eratosthenes(B)
    print(factor_base)

    # find b-smooth numbers
    x = 0
    root_n = math.ceil(math.sqrt(n))
    exponent_vectors = None
    exponent_vectors_actual = None
    exponent_as = None
    while True:
        # generate a and b
        a = x + root_n
        b = a**2 % n
        print(f"a: {a}, b: {b}, x: {x}")

        # check if b is B-smooth
        factors = is_B_smooth_trial_division(factor_base, b)
        if factors is None:
            x = incriment_x(x)
            continue
        print(f"factors: {factors}")

        # generate exponent vector
        exponent_vector = np.zeros(len(factor_base), dtype=int)
        exponent_vector_actual = np.zeros(len(factor_base), dtype=int)
        for i, p in enumerate(factor_base):
            exponent_vector[i] = factors.count(p) % 2
            exponent_vector_actual[i] = factors.count(p)

        # add exponent vector to space
        if not np.all(exponent_vector == 0):
            if exponent_vectors is None:
                exponent_vectors = exponent_vector
                exponent_vectors_actual = exponent_vector_actual
                exponent_as = np.array([a])
            elif exponent_vectors.ndim == 1:
                if not np.all(exponent_vectors == exponent_vector):
                    exponent_vectors = np.vstack((exponent_vectors, exponent_vector))
                    exponent_vectors_actual = np.vstack(
                        (exponent_vectors_actual, exponent_vector_actual)
                    )
                    exponent_as = np.append(exponent_as, a)
            elif not np.any(np.all(exponent_vectors == exponent_vector, axis=1)):
                exponent_vectors = np.vstack((exponent_vectors, exponent_vector))
                exponent_vectors_actual = np.vstack(
                    (exponent_vectors_actual, exponent_vector_actual)
                )
                exponent_as = np.append(exponent_as, a)
            # print(f"exponent_vector:\n {exponent_vector}")
        print(f"exponent_vectors:\n {exponent_vectors}")

        # check if we have enough vectors
        if (
            exponent_vectors is not None
            and exponent_vectors.ndim == 2
            and exponent_vectors.shape[0]
            == NUMBER_OF_VECTORS_IN_SPACE_MOD * len(factor_base)
        ):
            break

        x = incriment_x(x)

        # key = msvcrt.getch()
        # if key == b'q':
        #     return 0

    # now we have enough vectors, find the null space
    print(f"exponent_vectors:\n {exponent_vectors}")
    print(f"exponent_vectors_actual:\n {exponent_vectors_actual}")
    print(f"exponent_as:\n {exponent_as}")

    ns_sympy = Matrix(exponent_vectors.T).nullspace()
    ns = np.zeros((len(ns_sympy[0]), len(ns_sympy)))
    for i in range(len(ns_sympy)):
        numpy_array_from_sympy = np.array(ns_sympy[i]).astype(np.float64)
        ns = np.hstack((ns, numpy_array_from_sympy))
    ns = ns.T

    # print(f"null_space:\n {ns}")
    print(f"null_space size: {ns.shape}")

    # for all null space vectors
    for i in range(ns.shape[0]):
        print(f"\n{i} ======================================")
        threshold = 1e-10
        ns_vector = ns[i]
        print(f"null_space_vector:\n {ns_vector}")
        ns_vector = np.where(np.abs(ns_vector) < threshold, 0, 1).flatten()
        print(f"null_space:\n {ns_vector}")

        as_product = calculate_a_product(ns_vector, exponent_as)
        primes_product = calculate_primes_product(
            ns_vector, exponent_vectors_actual, factor_base
        )

        print(f"as_product: {as_product}")
        print(f"primes_product: {primes_product}")

        factor = math.gcd(as_product - primes_product, n)
        print(f"factor: {factor}")

        if factor != 1 and factor != n and factor != -1:
            return factor


if __name__ == "__main__":
    n = 15347
    B = get_smoothness_bound(n)
    f = sieve(n)
    print("\n=========")
    print(f"number: {n}")
    print(f"B: {B}")
    print(f"factor_base: {sieve_of_eratosthenes(B)}")
    print(f"factor: {f}")

    # TODO: doesnt work for big numbers for some reason
