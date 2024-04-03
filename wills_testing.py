import math
import msvcrt
import numpy as np
from scipy.linalg import null_space
# Steps for the quadratic sieve:
# 1. Choose a smoothness bound B. The number π(B), denoting the number of prime numbers less than B, will control both the length of the vectors and the number of vectors needed.
# 2. Use sieving to locate π(B) + 1 numbers ai such that bi = (ai2 mod n) is B-smooth.
# 3. Factor the bi and generate exponent vectors mod 2 for each one.
# 4. Use linear algebra to find a subset of these vectors which add to the zero vector. Multiply the corresponding ai together and give the result mod n the name a; similarly, multiply the bi together which yields a B-smooth square b2.
# 5. We are now left with the equality a2 = b2 mod n from which we get two square roots of (a2 mod n), one by taking the square root in the integers of b2 namely b, and the other the a computed in step 4.
# 6. We now have the desired identity. Compute the GCD of n with the difference (or sum) of a and b. This produces a factor, although it may be a trivial factor (n or 1). If the factor is trivial, try again with a different linear dependency or different a.

def incriment_x(x):
    if x > 0:
        return -x
    else:
        return -x + 1

def sieve(n):
    # Choose smoothness bound B
    # online found that B is approximately equal to exp(sqrt(log(n) * log(log(n))))
    B = math.floor(math.exp(math.sqrt(math.log(n) * math.log(math.log(n)))))
    
    # Now we have B compute primes up to B using Sieve of Eratosthenes
    factor_base = []
    is_primes = [True] * (B + 1)
    for i in range(2, math.ceil(math.sqrt(B + 1))):
        if is_primes[i]:
            factor_base.append(i)
            for j in range(i*i, B+1, i):
                is_primes[j] = False
    for i in range(math.ceil(math.sqrt(B + 1)), B + 1):
        if is_primes[i]:
            factor_base.append(i)
    print(factor_base)

    # find b-smooth numbers
    x = 0
    root_n = math.ceil(math.sqrt(n))
    exponent_vectors = None
    exponent_vectors_actual = None
    exponent_as = None
    while True:
        a = x + root_n
        b = a**2 % n

        print(f"a: {a}, b: {b}, x: {x}")

        #do trial divison (replace with toneli shanks)
        factors = []
        b1 = b

        if b1 == 0:
            x = incriment_x(x)
            continue

        for p in factor_base:
            print(f"p: {p}, b1: {b1}")
            while b1 / p == b1 // p:
                b1 = b1 // p
                factors.append(p)

        print(f"factors: {factors}, b1: {b1}")
        
        if b1 == 1 and len(factors) > 0:
            exponent_vector = np.zeros(len(factor_base), dtype=int)
            exponent_vector_actual = np.zeros(len(factor_base), dtype=int)
            for i, p in enumerate(factor_base):
                exponent_vector[i] = factors.count(p) % 2
                exponent_vector_actual[i] = factors.count(p)
            # if exponent vector is not all 0, add it to the list of exponent vectors
            if not np.all(exponent_vector == 0):
                # if vector is already in the list, don't add it
                if exponent_vectors is None:
                    exponent_vectors = exponent_vector
                    exponent_vectors_actual = exponent_vector_actual
                    exponent_as = np.array([a])
                elif exponent_vectors.ndim == 1:
                    if not np.all(exponent_vectors == exponent_vector):
                        exponent_vectors = np.vstack((exponent_vectors, exponent_vector))
                        exponent_vectors_actual = np.vstack((exponent_vectors_actual, exponent_vector_actual))
                        exponent_as = np.append(exponent_as, a)
                elif not np.any(np.all(exponent_vectors == exponent_vector, axis=1)):
                    exponent_vectors = np.vstack((exponent_vectors, exponent_vector))
                    exponent_vectors_actual = np.vstack((exponent_vectors_actual, exponent_vector_actual))
                    exponent_as = np.append(exponent_as, a)
                print(f"exponent_vector: {exponent_vector}")
        
        print(f"exponent_vectors:\n {exponent_vectors}")

        if exponent_vectors is not None and exponent_vectors.ndim == 2 and exponent_vectors.shape[0] == len(factor_base) + 1:
            break

        x = incriment_x(x)
        
        # key = msvcrt.getch()
        # if key == b'q':
        #     return 0

    print(f"exponent_vectors:\n {exponent_vectors}")
    print(f"exponent_vectors_actual:\n {exponent_vectors_actual}")
    print(f"exponent_as:\n {exponent_as}")
    ns = null_space(exponent_vectors.T)
    ns = ns.T[0]
    threshold = 1e-10
    ns = np.where(np.abs(ns) < threshold, 0, 1).flatten()
    print(f"null_space:\n {ns}")

    as_product = 1
    primes_product = 1

    for i, a in enumerate(ns):
        if a == 1:
            as_product *= exponent_as[i]
            as_product = as_product % n

    prime_power_vector = np.zeros(len(factor_base), dtype=int)
    for i, a in enumerate(ns):
        if a == 1:
            prime_power_vector += exponent_vectors_actual[i]

    prime_power_vector = prime_power_vector // 2
    
    print(f"prime_power_vector: {prime_power_vector}")

    for i, p in enumerate(factor_base):
        primes_product *= p**prime_power_vector[i]
        primes_product = primes_product % n

    print(f"as_product: {as_product}")
    print(f"primes_product: {primes_product}")

    factor = math.gcd(as_product - primes_product, n)
    return factor

if __name__ == "__main__":
    f = sieve(589)
    print(f"factor: {f}")

    #TODO: if nontrivial factor repeat...