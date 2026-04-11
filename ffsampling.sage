"""This file contains important algorithms for Falcon.

- the Fast Fourier orthogonalization (in coefficient and FFT representation)
- the Fast Fourier nearest plane (in coefficient and FFT representation)
- the Fast Fourier sampling (only in FFT)
.
"""
from sage.all import CDF
from math import floor

load("fft.sage")
load("common.sage")
load("samplerz.sage")
#from common import split, merge                         # Split, merge
#from fft import add, sub, mul, div, adj                 # Operations in coef.
#from fft import add_fft, sub_fft, mul_fft, div_fft, adj_fft  # Ops in FFT
#from fft import split_fft, merge_fft, fft_ratio         # FFT
#from samplerz import samplerz                           # Gaussian sampler in Z

# Helpers
def clean_real(x):
    """Extract a stable real value from CDF."""
    return float(CDF(x).real())


def falcon_round(x):
    """Round like Falcon (nearest integer, ties away from 0)."""
    x = float(x)
    if x >= 0:
        return int(floor(x + 0.5))
    else:
        return -int(floor(-x + 0.5))

# Gram matrix
def gram(B):
    """Compute the Gram matrix of B.

    Args:
        B: a matrix

    Format: coefficient
    """
    rows = range(len(B))
    ncols = len(B[0])
    deg = len(B[0][0])
    G = [[[0 for coef in range(deg)] for j in rows] for i in rows]
    for i in rows:
        for j in rows:
            for k in range(ncols):
                G[i][j] = add(G[i][j], mul(B[i][k], adj(B[j][k])))
    return G

# LDL (coefficient) 
def ldl(G):
    """
    Compute the LDL decomposition of G. Only works with 2 * 2 matrices.

    Args:
        G: a Gram matrix

    Format: coefficient

    Corresponds to algorithm 8 (LDL*) of Falcon's documentation,
    except it's in polynomial representation.
    """
    deg = len(G[0][0])
    assert len(G) == 2 and len(G[0]) == 2

    zero = [0] * deg
    one = [1] + [0] * (deg - 1)
    D00 = G[0][0][:]

    
    
    # Check for zeros in denominator
    zero_count = sum(1 for x in G[0][0] if abs(x) < 1e-12)
    if zero_count > 0:
        print("Warning: denominator G[0][0] has zero entries:", zero_count)

    L10 = div(G[1][0], G[0][0])
    D11 = sub(G[1][1], mul(mul(L10, adj(L10)), G[0][0]))
    L = [[one, zero], [L10, one]]
    D = [[D00, zero], [zero, D11]]

    return [L, D]


def ldl_fft(G):
    """
    Compute the LDL decomposition of G. Only works with 2 * 2 matrices.

    Args:
        G: a Gram matrix

    Format: FFT

    Corresponds to algorithm 8 (LDL*) of Falcon's documentation.
    """
    deg = len(G[0][0])
    assert len(G) == 2 and len(G[0]) == 2

    zero = [CDF(0)] * deg
    one  = [CDF(1)] * deg
    D00 = G[0][0][:]
    
    # Check for zeros in denominator
    zero_count = sum(1 for x in G[0][0] if abs(x) < 1e-12)
    if zero_count > 0:
        print("Warning: denominator G[0][0] has zero entries:", zero_count)

    L10 = div_fft(G[1][0], G[0][0])
    D11 = sub_fft(G[1][1], mul_fft(mul_fft(L10, adj_fft(L10)), G[0][0]))
    L = [[one, zero], [L10, one]]
    D = [[D00, zero], [zero, D11]]

    return [L, D]


def ffldl(G):
    """Compute the ffLDL decomposition tree of G.

    Args:
        G: a Gram matrix

    Format: coefficient

    Corresponds to algorithm 9 (ffLDL) of Falcon's documentation,
    except it's in polynomial representation.
    """
    n = len(G[0][0])
    L, D = ldl(G)
    # Coefficients of L, D are elements of R[x]/(x^n - x^(n/2) + 1), in coefficient representation
    if (n > 2):
        # A bisection is done on elements of a 2*2 diagonal matrix.
        d00, d01 = split(D[0][0])
        d10, d11 = split(D[1][1])
        G0 = [[d00, d01], [adj(d01), d00]]
        G1 = [[d10, d11], [adj(d11), d10]]
        return [L[1][0], ffldl(G0), ffldl(G1)]
    elif (n == 2):
        # Bottom of the recursion.
        D[0][0][1] = 0
        D[1][1][1] = 0
        return [L[1][0], D[0][0], D[1][1]]


def ffldl_fft(G):
    """Compute the ffLDL decomposition tree of G.

    Args:
        G: a Gram matrix

    Format: FFT

    Corresponds to algorithm 9 (ffLDL) of Falcon's documentation.
    """
    n = len(G[0][0]) * fft_ratio
    L, D = ldl_fft(G)
    # Coefficients of L, D are elements of R[x]/(x^n - x^(n/2) + 1), in FFT representation
    if (n > 2):
        # A bisection is done on elements of a 2*2 diagonal matrix.
        d00, d01 = split_fft(D[0][0])
        d10, d11 = split_fft(D[1][1])
        G0 = [[d00, d01], [adj_fft(d01), d00]]
        G1 = [[d10, d11], [adj_fft(d11), d10]]
        return [L[1][0], ffldl_fft(G0), ffldl_fft(G1)]
    elif (n == 2):
        # End of the recursion (each element is real).
        return [L[1][0], D[0][0], D[1][1]]

# Fast Fourier Nearest Plane
def ffnp(t, T):
    """Compute the ffnp reduction of t, using T as auxilary information.

    Args:
        t: a vector
        T: a ldl decomposition tree

    Format: coefficient
    """
    n = len(t[0])
    z = [None, None]
    if (n > 1):
        l10, T0, T1 = T
        z[1] = merge(ffnp(split(t[1]), T1))
        t0b = add(t[0], mul(sub(t[1], z[1]), l10))
        z[0] = merge(ffnp(split(t0b), T0))
        return z
    elif n == 1:
        z[0] = [falcon_round(t[0][0])]
        z[1] = [falcon_round(t[1][0])]
        return z


def ffnp_fft(t, T):
    """Compute the ffnp reduction of t, using T as auxilary information.

    Args:
        t: a vector
        T: a ldl decomposition tree

    Format: FFT
    """
    n = len(t[0]) * fft_ratio
    z = [0, 0]

    if n > 1:
        l10, T0, T1 = T

        z[1] = merge_fft(ffnp_fft(split_fft(t[1]), T1))
        t0b = add_fft(t[0], mul_fft(sub_fft(t[1], z[1]), l10))
        z[0] = merge_fft(ffnp_fft(split_fft(t0b), T0))

        return z
    elif n == 1:
        z[0] = [falcon_round(clean_real(t[0][0]))]
        z[1] = [falcon_round(clean_real(t[1][0]))]
        return z


def ffsampling_fft(t, T, sigmin, randombytes):
    """Compute the ffsampling of t, using T as auxilary information.

    Args:
        t: a vector
        T: a ldl decomposition tree

    Format: FFT

    Corresponds to algorithm 11 (ffSampling) of Falcon's documentation.
    """
    #print("T: ", T[-2:])
    

    n = len(t[0]) * fft_ratio
    
    print("ENTER ffsampling:")
    print("n =", n)
    #print("len(T) =", len(T))
    #print(f"{'  '*depth}Entering ffsampling_fft, n={n}, depth={depth}")
    z = [0, 0]
    #print("T in FFsampling: ", T)
    if len(T) == 2:
        print("Debug n: ", n)
    if n > 1:
        if n <= 4:
            print("T :", T)
        #print("if n > 1 T: ", T)
        l10, T0, T1 = T
        print("Debug Z: ", z)
        # Sample second coordinate first
        print("Calling z1 with:")
        print("n =", n)
        print("len(T1) =", len(T1))
        #print("T1 =", T1)
        z1 = ffsampling_fft(split_fft(t[1]), T1, sigmin, randombytes)
        print("Z1: ", z1)
        z[1] = merge_fft(z1)

        # Compute t0'
        t0b = add_fft(t[0], mul_fft(sub_fft(t[1], z[1]), l10))

        # Sample first coordinate
        z0 = ffsampling_fft(split_fft(t0b), T0, sigmin, randombytes)
        z[0] = merge_fft(z0)
        
        return z

    elif n == 1:
        # Leaf: scalar Gaussian sampling
        print("t[0][0]", t[0][0])
        print("t[1][0]", t[1][0])
        mu0 = clean_real(t[0][0])
        mu1 = clean_real(t[1][0])

        print("Debug mu0:", mu0)
        print("Debug mu1:", mu1)
        # IMPORTANT: Falcon uses T[0] as sigma directly
        d00 = T[0]
        print("T0: ", T[0])
        sigma = (d00[0] if isinstance(d00, list) else d00)
        print("Sigma: ", sigma)
        print("Sigmin: ", sigmin)
        
        z[0] = [samplerz(mu0, sigma, sigmin, randombytes)]
        z[1] = [samplerz(mu1, sigma, sigmin, randombytes)]
    
        return z


from sage.all import CDF, mean, variance
import os
import statistics

# Small random polynomial helper
def _random_poly(n):
    from sage.all import ZZ
    import random
    return [ZZ(random.randint(-5, 5)) for _ in range(n)]

# Minimal ffLDL tree for testing
def _fake_tree(n):
    from sage.all import CDF
    if n == 1:
        return [
            [CDF(0)],        # l10 (not used at leaf)
            [CDF(1.8205)],   # d00 (sigma)
            [CDF(1.8205)]    # d11 (sigma)
        ]
    else:
        return [
            [CDF(0)] * n,
            _fake_tree(n // 2),
            _fake_tree(n // 2)
        ]
"""
# ----------------------------
# TEST 1: Basic execution
# ----------------------------
def test_ffsampling_basic():
    print("DEBUG: Running ffsampling_basic with sigmin = 1.2")
    n = 4
    t0 = [CDF(0.1 * i) for i in range(n)]
    t1 = [CDF(-0.2 * i) for i in range(n)]
    t = [t0, t1]

    T = _fake_tree(n)
    sigmin = 1.2

    z = ffsampling_fft(t, T, sigmin, os.urandom)

    print("DEBUG: z[0] =", z[0])
    print("DEBUG: z[1] =", z[1])
    print("✅ test_ffsampling_basic executed")

# ----------------------------
# TEST 2: Samplerz distribution
# ----------------------------
def test_ffsampling_distribution():
    print("DEBUG: Running samplerz distribution test")
    mu = 0.5
    sigma = 1.5
    sigmin = 1.2

    samples = [samplerz(mu, sigma, sigmin, os.urandom) for _ in range(1000)]
    

    average = mean(samples)
    print("Mean =", average)
    var = variance(samples)
    print("Variance =", var)

    print(f"Mean ≈ {float(average)} (expected {mu})")
    print(f"Variance ≈ {float(var)} (expected {sigma**2})")
    print("✅ test_ffsampling_distribution executed")

# ----------------------------
# TEST 3: Norm sanity
# ----------------------------
def test_norm_bound():
    print("DEBUG: Running norm sanity test")
    n = 4
    t0 = [CDF(0.1 * i) for i in range(n)]
    t1 = [CDF(-0.2 * i) for i in range(n)]
    t = [t0, t1]
    T = _fake_tree(n)
    sigmin = 1.2

    z = ffsampling_fft(t, T, sigmin, os.urandom)

    # Use absolute value sum instead of rounding integers
    norm = sum(abs(x.real())**2 + abs(x.imag())**2 for x in z[0] + z[1])
    print("DEBUG: Norm =", norm)
    print("✅ test_norm_bound executed")

from sage.all import ZZ, round
def test_fft_roundtrip(n, iterations=5):
    for _ in range(iterations):
        f = [randint(-5, 5) for _ in range(n)]
        f_fft = fft(f)
        f_rec = ifft(f_fft)

        print("f      =", f)
        print("f_rec  =", f_rec)

        f_round = [falcon_round(x) for x in f_rec]
        print("rounded=", f_round)

        if f_round != f:
            print("❌ Roundtrip failed")
            return False
    print("✅ Roundtrip passed")
    return True


#L = [1,0,0,0,0,0,0,0]
#print(ifft(fft(L)))


def test_fft(n, iterations=10):
    for i in range(iterations):
        f = [(randint(-3, 4)) for j in range(n)]
        g = [(randint(-3, 4)) for j in range(n)]
        print("f = ", f)
        print("g = ", g)
        h = mul(f, g)
        k = div(h, f)
        print("h = mul f g = ", h)
        print("k div h/f = ", k)
        k_int = [falcon_round(elt) for elt in k]  # use Sage rounding
        if k_int != g:
            print("(f * g) / f =", k_int)
            print("g =", g)
            print("mismatch")
            return False
    return True

"""
"""
# ============================================================
# RUN TESTS
# ============================================================
if __name__ == "__main__":
    n = 8  # example polynomial length, adjust as needed
    iterations = 10  # number of test iterations
    #test_ffsampling_basic()
    #test_ffsampling_distribution()
    #test_norm_bound()

    print("=== Running FFT test ===")
    if test_fft_roundtrip(n, iterations):
        print("✅ FFT roundtrip test passed")
    else:
        print("❌ FFT roundtrip test failed")

    if test_fft(n, iterations):
        print("✅ FFT test passed")
    else:
        print("❌ FFT test failed")

    print("\n🎉 All tests executed (integer rounding checks skipped)")
    
    def schoolbook_mul(f, g):
        n = len(f)
        res = [0]*n
        for i in range(n):
            for j in range(n):
                k = (i + j) % n
                sign = -1 if (i + j) >= n else 1  # for mod (x^n + 1)
                res[k] += sign * f[i] * g[j]
        return res

    f = [randint(-3,3) for _ in range(n)]
    g = [randint(-3,3) for _ in range(n)]

    h_fft = mul(f, g)
    h_true = schoolbook_mul(f, g)

    print("fft mul =", [falcon_round(x) for x in h_fft])
    print("true    =", h_true)
"""