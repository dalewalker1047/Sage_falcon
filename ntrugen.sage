"""
This file implements the section 3.8.2 of Falcon's documentation.
"""

from sage.all import ZZ, CDF, xgcd

load("fft.sage")
load("ntt.sage")
load("common.sage")
load("samplerz.sage")

#from fft import fft, ifft, add_fft, mul_fft, adj_fft, div_fft
#from fft import add, mul, div, adj
#from ntt import ntt
#from common import sqnorm, safe_int
#from samplerz import samplerz


q = ZZ(12 * 1024 + 1)


def karatsuba(a, b, n):
    """
    Karatsuba multiplication between polynomials.
    The coefficients may be either integer or real.
    """
    if n == 1:
        return [ZZ(a[0]) * ZZ(b[0]), ZZ(0)]
    else:
        n2 = n // 2
        a0 = a[:n2]
        a1 = a[n2:]
        b0 = b[:n2]
        b1 = b[n2:]
        ax = [ZZ(a0[i]) + ZZ(a1[i]) for i in range(n2)]
        bx = [ZZ(b0[i]) + ZZ(b1[i]) for i in range(n2)]
        a0b0 = karatsuba(a0, b0, n2)
        a1b1 = karatsuba(a1, b1, n2)
        axbx = karatsuba(ax, bx, n2)
        for i in range(n):
            axbx[i] -= (a0b0[i] + a1b1[i])
        ab = [ZZ(0)] * (2 * n)
        for i in range(n):
            ab[i] += a0b0[i]
            ab[i + n] += a1b1[i]
            ab[i + n2] += axbx[i]
        return ab


def karamul(a, b):
    """
    Karatsuba multiplication, followed by reduction mod (x ** n + 1).
    """
    n = len(a)
    ab = karatsuba(a, b, n)
    abr = [ab[i] - ab[i + n] for i in range(n)]
    return abr


def galois_conjugate(a):
    """
    Galois conjugate of an element a in Q[x] / (x ** n + 1).
    Here, the Galois conjugate of a(x) is simply a(-x).
    """
    return [((-1) ** i) * ZZ(a[i]) for i in range(len(a))]


def field_norm(a):
    """
    Project an element a of Q[x] / (x ** n + 1) onto Q[x] / (x ** (n // 2) + 1).
    Only works if n is a power-of-two.
    """
    n2 = len(a) // 2
    ae = [ZZ(a[2 * i]) for i in range(n2)]
    ao = [ZZ(a[2 * i + 1]) for i in range(n2)]
    ae_sq = karamul(ae, ae)
    ao_sq = karamul(ao, ao)

    res = ae_sq[:]
    for i in range(n2 - 1):
        res[i + 1] -= ao_sq[i]
    res[0] += ao_sq[n2 - 1]

    return res


def lift(a):
    """
    Lift an element a of Q[x] / (x ** (n // 2) + 1) up to Q[x] / (x ** n + 1).
    The lift of a(x) is simply a(x ** 2) seen as an element of Q[x] / (x ** n + 1).
    """
    n = len(a)
    res = [ZZ(0)] * (2 * n)
    for i in range(n):
        res[2 * i] = ZZ(a[i])
    return res


def bitsize(a):
    """
    Compute the bitsize of an element of Z (not counting the sign).
    The bitsize is rounded to the next multiple of 8.
    This makes the function slightly imprecise, but faster to compute.
    """
    val = abs(ZZ(a))
    res = 0
    while val:
        res += 8
        val >>= 8
    return res


def reduce(f, g, F, G):
    """
    Reduce (F, G) relatively to (f, g).

    This is done via Babai's reduction.
    (F, G) <-- (F, G) - k * (f, g), where k = round((F f* + G g*) / (f f* + g g*)).
    Corresponds to algorithm 7 (Reduce) of Falcon's documentation.
    """
    n = len(f)
    size = max(53, bitsize(min(f)), bitsize(max(f)), bitsize(min(g)), bitsize(max(g)))

    f_adj = [ZZ(elt) >> (size - 53) for elt in f]
    g_adj = [ZZ(elt) >> (size - 53) for elt in g]

    fa_fft = fft(f_adj)
    ga_fft = fft(g_adj)

    while(1):
        # Because we work in finite precision to reduce very large polynomials,
        # we may need to perform the reduction several times.
        Size = max(53, bitsize(min(F)), bitsize(max(F)), bitsize(min(G)), bitsize(max(G)))
        if Size < size:
            break

        F_adj = [ZZ(elt) >> (Size - 53) for elt in F]
        G_adj = [ZZ(elt) >> (Size - 53) for elt in G]

        Fa_fft = fft(F_adj)
        Ga_fft = fft(G_adj)

        den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)), mul_fft(ga_fft, adj_fft(ga_fft)))
        num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)), mul_fft(Ga_fft, adj_fft(ga_fft)))
        k_fft = div_fft(num_fft, den_fft)
        k = ifft(k_fft)
        k = [safe_int(elt) for elt in k]
        if all(elt == 0 for elt in k):
            break
        # The two next lines are the costliest operations in ntru_gen
        # (more than 75% of the total cost in dimension n = 1024).
        # There are at least two ways to make them faster:
        # - replace Karatsuba with Toom-Cook
        # - mutualized Karatsuba, see ia.cr/2020/268
        # For simplicity reasons, we didn't implement these optimisations here.
        fk = karamul(f, k)
        gk = karamul(g, k)
        for i in range(n):
            F[i] -= fk[i] << (Size - size)
            G[i] -= gk[i] << (Size - size)
    return F, G

"""
Replaced by built in xgcd
def xgcd(b, n):
    Compute the extended GCD of two integers b and n.
    Return d, u, v such that d = u * b + v * n, and d is the GCD of b, n.

    b = ZZ(b)
    n = ZZ(n)

    x0, x1 = ZZ(1), ZZ(0)
    y0, y1 = ZZ(0), ZZ(1)

    while n != 0:
        q_, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q_ * x1
        y0, y1 = y1, y0 - q_ * y1

    return b, x0, y0
"""


def ntru_solve(f, g):
    """
    Solve the NTRU equation for f and g.
    Corresponds to NTRUSolve in Falcon's documentation.
    """
    n = len(f)
    if n == 1:
        f0 = ZZ(f[0])
        g0 = ZZ(g[0])
        d, u, v = xgcd(f0, g0)
        if d != 1:
            raise ValueError
        else:
            return [- q * v], [q * u]
    else:
        fp = field_norm(f)
        gp = field_norm(g)
        Fp, Gp = ntru_solve(fp, gp)
        F = karamul(lift(Fp), galois_conjugate(g))
        G = karamul(lift(Gp), galois_conjugate(f))
        F, G = reduce(f, g, F, G)
        return F, G


def gs_norm(f, g, q):
    """
    Compute the squared Gram-Schmidt norm of the NTRU matrix generated by f, g.
    This matrix is [[g, - f], [G, - F]].
    This algorithm is equivalent to line 9 of algorithm 5 (NTRUGen).
    """
    sqnorm_fg = sqnorm([f, g])
    ffgg = add(mul(f, adj(f)), mul(g, adj(g)))
    Ft = div(adj(g), ffgg)
    Gt = div(adj(f), ffgg)
    Ft_int = [safe_int(x) for x in Ft]
    Gt_int = [safe_int(x) for x in Gt]
    sqnorm_FG = (q ** 2) * sqnorm([Ft_int, Gt_int])
    return max(sqnorm_fg, sqnorm_FG)


def gen_poly(n):
    """
    Generate a polynomial of degree at most (n - 1), with coefficients
    following a discrete Gaussian distribution D_{Z, 0, sigma_fg} with
    sigma_fg = 1.17 * sqrt(q / (2 * n)).
    """
    # 1.17 * sqrt(12289 / 8192)
    sigma = 1.43300980528773
    assert(n < 4096)
    f0 = [samplerz(0, sigma, sigma - 0.001) for _ in range(4096)]
    f = [ZZ(0)] * n
    k = 4096 // n
    for i in range(n):
        # We use the fact that adding k Gaussian samples of std. dev. sigma
        # gives a Gaussian sample of std. dev. sqrt(k) * sigma.
        f[i] = sum(ZZ(f0[i * k + j]) for j in range(k))
    return f


def ntru_gen(n):
    """
    Implement the algorithm 5 (NTRUGen) of Falcon's documentation.
    At the end of the function, polynomials f, g, F, G in Z[x]/(x ** n + 1)
    are output, which verify f * G - g * F = q mod (x ** n + 1).
    """
    for _ in range(10000):
    #while True:
        f = gen_poly(n)
        g = gen_poly(n)
        if gs_norm(f, g, q) > (1.17 ** 2) * q:
            continue
        f_ntt = ntt(f)
        if any((elem == 0) for elem in f_ntt):
            continue
        try:
            F, G = ntru_solve(f, g)
            F = [ZZ(coef) for coef in F]
            G = [ZZ(coef) for coef in G]
            return f, g, F, G
        # If the NTRU equation cannot be solved, a ValueError is raised
        # In this case, we start again
        except ValueError:
            continue
    raise RuntimeError("NTRU_Gen stuck")


def test_ntt(n, iterations=10):
    """Test the NTT."""
    for i in range(iterations):
        f = [randint(0, q - 1) for j in range(n)]
        g = [randint(0, q - 1) for j in range(n)]
        h = mul_zq(f, g)
        try:
            k = div_zq(h, f)
            if k != g:
                print("(f * g) / f =", k)
                print("g =", g)
                print("mismatch")
                return False
        except ZeroDivisionError:
            continue
    return True


def check_ntru(f, g, F, G):
    """Check that f * G - g * F = q mod (x ** n + 1)."""
    a = karamul(f, G)
    b = karamul(g, F)
    c = [a[i] - b[i] for i in range(len(f))]
    return ((c[0] == q) and all(coef == 0 for coef in c[1:]))


def test_ntrugen(n, iterations=10):
    """Test ntru_gen."""
    for i in range(iterations):
        f, g, F, G = ntru_gen(n)
        if check_ntru(f, g, F, G) is False:
            return False
    return True
"""
if __name__ == "__main__":
    load("ntt.sage")
    load("common.sage")
    #from ntt import mul_zq, sub_zq, ntt
    
    #from common import R

    def verify_ntru(n):
        print(f"\n=== Testing n = {n} ===")

        f, g, F, G = ntru_gen(n)

        print("Generated polynomials.")

        # -----------------------
        # 1. Check NTRU equation
        # -----------------------
        fG = mul_zq(f, G)
        gF = mul_zq(g, F)
        lhs = sub_zq(fG, gF)

        expected = [R(q)] + [R(0)] * (n - 1)

        if lhs == expected:
            print("✅ NTRU equation holds")
        else:
            print("❌ NTRU equation FAILED")
            print("lhs:", lhs[:10], "...")
            print("expected:", expected[:10], "...")
            return False

        # -----------------------
        # 2. Gram-Schmidt norm
        # -----------------------
        norm_val = gs_norm(f, g, q)
        bound = (1.17 ** 2) * q

        print(f"GS norm: {norm_val}")
        print(f"Bound:   {bound}")

        if norm_val <= bound:
            print("✅ Gram-Schmidt bound OK")
        else:
            print("❌ Gram-Schmidt bound FAILED")
            return False

        # -----------------------
        # 3. Check NTT(f)
        # -----------------------
        f_ntt = ntt(f)

        if any(elem == 0 for elem in f_ntt):
            print("❌ NTT(f) has zero")
            return False
        else:
            print("✅ NTT(f) invertible")

        print("🎉 ALL TESTS PASSED")
        return True

    
    # -----------------------
    # Run tests
    # -----------------------

    # Start small (fast debug)
    for n in [1024]:
        verify_ntru(n)

    # Uncomment when ready for full test
    # verify_ntru(256)
    # verify_ntru(512)
    # verify_ntru(1024)
    
    n = 8  # example polynomial length, adjust as needed
    iterations = 10  # number of test iterations

    print("\n=== Running NTT test ===")
    if test_ntt(n, iterations):
        print("✅ NTT test passed")
    else:
        print("❌ NTT test failed")

    print("\n=== Running NTRU generation test ===")
    if test_ntrugen(n, iterations):
        print("✅ NTRU gen test passed")
    else:
        print("❌ NTRU gen test failed")
"""