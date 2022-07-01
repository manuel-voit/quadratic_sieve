from math import sqrt, isqrt, ceil, prod
from itertools import chain

import miller_rabin
import trial_divison


def gen_f(b):
    """Generate the factor base.

    Method: Sieve of Eratosthenes

    Input:
        b: Smoothness bound B

    Output:
        f: Factor base F
    """
    if b < 2:
        return []

    is_prime = [True] * (b + 1)
    is_prime[0] = False
    is_prime[1] = False

    idx = 2

    while idx * idx <= b:
        if is_prime[idx]:
            for i in range(idx * idx, b + 1, idx):
                is_prime[i] = False
        idx += 1

    f = []
    for i in range(2, b + 1):
        if is_prime[i]:
            f.append(i)

    return f


def find_smooth(n, factor_base):
    """Find B-smooth numbers.

    Method: Trial devision on factor base

    Input:
        n: Number N that needs to be factorized
        factor_base: Factor base F

    Output:
        smooth_nums: B-smooth numbers
        fac_smooth_nums: Factorization of the B-smooth numbers
        xlist: x values for congruence of squares
    """
    smooth_nums = []
    fac_smooth_nums = []
    xlist = []

    num_relations = len(factor_base) + 1
    i = 0
    x = ceil(sqrt(n)) - 1

    while i < num_relations:
        x += 1
        psn = pow(x, 2, n)
        b_smooth, factors = trial_division_with_f(psn, factor_base)

        if b_smooth:
            smooth_nums.append(psn)
            fac_smooth_nums.append(factors)
            xlist.append(x)
            i += 1
    return smooth_nums, fac_smooth_nums, xlist


def trial_division_with_f(psn, factor_base):
    """Return wether psn is a B-smooth number or not. Return the factorization for B-smooth numbers. Return an empty list otherwise.

    Method: Trial division with a factor base F

    Input:
        psn: Potential B-smooth number; Needs to be checked and then factored
        factor_base: Factor base F

    Output:
        smooth: psn is smooth or not
        smooth == True:
            factors: Factorization of psn
        smooth == False:
            factors: Empty list
    """
    factors = []

    for prime in factor_base:
        while psn % prime == 0:
            factors.append(prime)
            psn //= prime
    if psn == 1:
        smooth = True
    else:
        smooth = False
        factors = []
    return smooth, factors


def gen_parity_matrix(fac_smooth_nums, factor_base):
    """Generate the parity matrix based on the parity of the exponents of the prime coefficients from the factor base.

    Input:
        fac_smooth_nums: Factorization of the B-smooth numbers
        factor_base: Factor base F

    Output:
        contains_square: Relations on the factor base contain square number
        contains_square == True:
            m: Empty list
            fac: Factorization of the square number
        contains_square == False:
            m: Parity matrix
            fac: Empty list
    """

    m = []
    len_f = len(factor_base)

    for i, fac in enumerate(fac_smooth_nums):
        exp_parity = [0] * len_f

        for j in range(len_f):
            if factor_base[j] in fac_smooth_nums[i]:
                exp_parity[j] = (fac_smooth_nums[i].count(factor_base[j])) % 2

        if 1 not in exp_parity:
            return True, [], fac
        else:
            pass

        m.append(exp_parity)

    return False, m, []


def gcd(a, b):
    """Return the greatest common divisor of a and b."""

    if b == 0:
        return a
    elif a >= b:
        return gcd(b, a % b)
    else:
        return gcd(b, a)


def transpose(m):
    """Transpase the matrix m."""

    m_t = []
    for i in range(len(m[0])):
        row_t = []
        for row in m:
            row_t.append(row[i])
        m_t.append(row_t)
    return m_t


def gauss_elim(m):
    """Solve the linear equation system given by m.

    Method: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    """

    m = transpose(m)
    marks = [False] * len(m[0])

    for i in range(len(m)):
        row = m[i]

        for num in row:
            if num == 1:
                j = row.index(num)
                marks[j] = True

                for k in chain(range(0, i), range(i + 1, len(m))):
                    if m[k][j] == 1:
                        for i in range(len(m[k])):
                            m[k][i] = (m[k][i] + row[i]) % 2
                break

    m = transpose(m)

    sol_rows = []
    for i in range(len(marks)):
        if not marks[i]:
            free_row = [m[i], i]
            sol_rows.append(free_row)

    if not sol_rows:
        return "No solution found. Try a different value for B."
    print("Found {} potential solutions.".format(len(sol_rows)))
    return sol_rows, marks, m


def solve_row(sol_rows, m, marks, k=0):
    solution_vec, indices = [], []
    free_row = sol_rows[k][0]
    for i in range(len(free_row)):
        if free_row[i] == 1:
            indices.append(i)
    for r in range(len(m)):
        for i in indices:
            if m[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break

    solution_vec.append(sol_rows[k][1])
    return solution_vec


def cong_sq(solution_vec, smooth_nums, xlist, n):
    """Solve the congruence of squares.

    Methods: Congruence of squares, gcd

    Input:
        solution_vec: Index numbers of the relations that need to be combined
        smooth_nums: B-smooth numbers
        xlist: x values for congruence of squares
        n: Number N that needs to be factorized

    Output:
        factor: One factor of N
    """

    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]

    y_2 = 1
    for num in solution_nums:
        y_2 *= num

    x = 1
    for num in x_nums:
        x *= num

    y = isqrt(y_2)
    x %= n
    y %= n

    factor = gcd(abs(x - y), n)
    return factor


def qs(n, b):
    """Return the prime factors of n.

    Method: Quadratic sieve

    Input:
        n: Number N that needs to be factorized
        b: Smoothness bound B
    Output:
        Factorization or error
    """

    if not isinstance(n, int) or n < 1:
        return "Input a positive integer for N."
    if not isinstance(b, int) or b < 2:
        return "Input a positive integer greater than or equal to 2 for B."

    print("\n")
    print("Checking if {} is a prime number...".format(n))
    if miller_rabin.is_prime(n):
        print("{} is a prime number!".format(n))
        return "Factor of {}: {}".format(n, n)
    print("{} is not a prime number.".format(n))
    print("\n")

    print("Checking if {} is a sqare number...".format(n))
    if sqrt(n).is_integer():
        print("{} is a sqare number!".format(n))
        return "Factors of {}: {}; {}".format(n, isqrt(n), isqrt(n))
    print("{} is not a sqaure number.".format(n))
    print("\n")

    print("Trying to find factors for {} with QS:".format(n))
    print("Generating {}-smooth factor base...".format(b))
    factor_base = gen_f(b)
    print(factor_base)
    print("\n")

    len_f = len(factor_base)
    num_relations = len_f + 1
    print("Looking for {} differen {}-smooth relations...".format(num_relations, b))
    smooth_nums, fac_smooth_nums, xlist = find_smooth(n, factor_base)
    print("Found {} B-smooth numbers:".format(len(smooth_nums)))
    print(smooth_nums)
    print("\n")

    print("Checking for square numbers among the relations...")
    contains_square, parity_matrix, fac = gen_parity_matrix(fac_smooth_nums, factor_base)

    if contains_square:
        print("Square number found!")
        print("Solving congruence of squares...")
        y_2 = prod(fac)
        idx = smooth_nums.index(y_2)
        x = xlist[idx] % n
        y = isqrt(prod(fac)) % n

        p = gcd(x - y, n)
        q = n // p

        if not miller_rabin.is_prime(p) or not miller_rabin.is_prime(q):
            print("{} is not a semi-prime!".format(n))
            print("\n")

            print("Trying to factor with trial division...")
            factors = trial_divison.factor(n)
            return "Factors of {}: {}".format(n, factors)

        return "Factors of {}: {}; {}".format(n, p, q)
    print("No square numbers were found among the relations.")
    print("\n")

    print("Generating exponent parity matrix...")
    for row in parity_matrix:
        print(row)
    print("\n")

    print("Performing Gaussian Elimination...")
    sol_rows, marks, m = gauss_elim(parity_matrix)
    solution_vec = solve_row(sol_rows, m, marks, 0)

    comb = [0] * len(smooth_nums)
    for i in solution_vec:
        comb[i] = 1
    print("Solution vector found: " + str(comb))
    print("\n")

    print("Solving congruence of squares...")
    p = cong_sq(solution_vec, smooth_nums, xlist, n)

    for k in range(1, len(sol_rows)):
        if p == 1 or p == n:
            print("Did not work. Trying different solution vector...")
            solution_vec = solve_row(sol_rows, m, marks, k)
            p = cong_sq(solution_vec, smooth_nums, xlist, n)
        else:
            print("Found non-trivial factors!")
            q = n // p
            return "Factors of {}: {}; {}".format(n, p, q)

    print("Could not find any non-trivial factors!")
    return "Try a different value for B."


if __name__ == "__main__":
    factorization = qs(539873, 17)
    print(factorization)
