import random


def is_prime(n, k=10):
    """Return if n is a probable prime number.

    Method: Miller-Rabin

    Input:
        n: Number n for that the primality test is performed
        k: Rounds
    """

    if n <= 1 or n == 4:
        return False

    if n <= 3:
        return True

    m = n - 1
    s = 0

    while m % 2 == 0:
        s += 1
        m //= 2

    for j in range(k):
        a = 2 + random.randint(1, n - 4)
        b = pow(a, m, n)

        if b != 1 and b != (n - 1):
            i = 1

            while i < s and b != (n - 1):
                b = (b * b) % n

                if b == 1:
                    return False

                i += 1

            if b != (n - 1):
                return False

        return True


if __name__ == "__main__":
    n = 53987
    b_prime = is_prime(n)
    print("\n")

    if b_prime:
        print("{} is very likely a prime number.".format(n))
    else:
        print("{} is not a prime number.".format(n))
