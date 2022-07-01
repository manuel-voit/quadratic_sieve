def factor(n):
    """Return the prime factors of n.

    Method: Trial divison

    Input:
        n: Number N that needs to be factorized

    Output:
        factors: Factorization of N
    """

    factors = []
    f = 2

    while n > 1:
        if n % f == 0:
            factors.append(f)
            n //= f
        else:
            f += 1

    return factors


if __name__ == "__main__":
        n = 330616
        factors = factor(n)
        print("\n")
        print("Factors of {}: {}".format(n, factors))
