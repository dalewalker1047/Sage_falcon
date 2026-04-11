"""This file contains an implementation of the FFT(SageMath version).

The FFT implemented here is for polynomials in R[x]/(phi), with:
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024

We use Sage's Complex Double Field (CDF) for numerical FFT values.
"""
#from sage.all import CDF

load("common.sage")
load("fft_constants.sage")
#from common import split, merge         # Import split and merge
#from fft_constants import roots_dict    # Import constants useful for the FFT


def split_fft(f_fft):
    """Split a polynomial f in two polynomials.

    Args:
        f: a polynomial

    Format: FFT

    Corresponds to algorithm 1 (splitfft_2) of Falcon's documentation.
    """
    n = len(f_fft)
    w = roots_dict[n]

    f0_fft = [CDF(0)] * (n // 2)
    f1_fft = [CDF(0)] * (n // 2)

    for i in range(n // 2):
        a = CDF(f_fft[2 * i])
        b = CDF(f_fft[2 * i + 1])

        f0_fft[i] = CDF(0.5) * (a + b)
        f1_fft[i] = CDF(0.5) * (a - b) * CDF(w[2 * i]).conjugate()

    return [f0_fft, f1_fft]


def merge_fft(f_list_fft):
    """Merge two or three polynomials into a single polynomial f.

    Args:
        f_list: a list of polynomials

    Format: FFT

    Corresponds to algorithm 2 (mergefft_2) of Falcon's documentation.
    """
    f0_fft, f1_fft = f_list_fft
    n = 2 * len(f0_fft)
    w = roots_dict[n]

    f_fft = [CDF(0)] * n

    for i in range(n // 2):
        w_i = CDF(w[2 * i])

        f_fft[2 * i]     = CDF(f0_fft[i]) + w_i * CDF(f1_fft[i])
        f_fft[2 * i + 1] = CDF(f0_fft[i]) - w_i * CDF(f1_fft[i])

    return f_fft


def fft(f):
    """Compute the FFT of a polynomial mod (x ** n + 1).

    Args:
        f: a polynomial (coefficient list)

    Format: input as coefficients, output as FFT (CDF elements)
    """
    n = len(f)
    if (n > 2):
        f0, f1 = split(f)
        f0_fft = fft(f0)
        f1_fft = fft(f1)
        return merge_fft([f0_fft, f1_fft])
    elif (n == 2):
        return [
            CDF(f[0]) + CDF(1j) * CDF(f[1]),
            CDF(f[0]) - CDF(1j) * CDF(f[1])
        ]
    else:
        raise ValueError("FFT input size must be >= 2 and a power of 2")

def ifft(f_fft):
    """Compute the inverse FFT of a polynomial mod (x ** n + 1).

    Args:
        f: a FFT of a polynomial (CDF)

    Format: input as FFT, output as coefficient list (real values)
    """
    n = len(f_fft)
    #print(f_fft[0], f_fft[1].conjugate())
    if (n > 2):
        f0_fft, f1_fft = split_fft(f_fft)
        f0 = ifft(f0_fft)
        f1 = ifft(f1_fft)
        return merge([f0, f1])
    elif (n == 2):
        a = CDF(f_fft[0])
        b = CDF(f_fft[1])
        #print("a: ", a)
        #print("b: ", b)
        return [
            (a + b) / 2,
            (a - b) / (2 * CDF(1j))
        ]
        #return [
        #    CDF(f_fft[0]).real(),
        #    CDF(f_fft[0]).imag()
        #]
    else:
        raise ValueError("FFT input size must be >= 2 and a power of 2")
    #return [x / 2 for x in f]
    
# Basic polynomial ops
def add(f, g):
    """Addition of two polynomials (coefficient representation)."""
    assert len(f) == len(g)
    return [f[i] + g[i] for i in range(len(f))]


def neg(f):
    """Negation of a polynomials (any representation)."""
    return [-f[i] for i in range(len(f))]


def sub(f, g):
    """Substraction of two polynomials (any representation)."""
    return add(f, neg(g))


def mul(f, g):
    """Multiplication of two polynomials (coefficient representation)."""
    return ifft(mul_fft(fft(f), fft(g)))


def div(f, g):
    """Division of two polynomials (coefficient representation)."""
    return ifft(div_fft(fft(f), fft(g)))


def adj(f):
    """Ajoint of a polynomial (coefficient representation)."""
    return ifft(adj_fft(fft(f)))

# FFT-domain ops

def add_fft(f_fft, g_fft):
    """Addition of two polynomials (FFT representation)."""
    return add(f_fft, g_fft)


def sub_fft(f_fft, g_fft):
    """Substraction of two polynomials (FFT representation)."""
    return sub(f_fft, g_fft)


def mul_fft(f_fft, g_fft):
    """Multiplication of two polynomials (coefficient representation). Pointwise multiplication"""
    return [CDF(f_fft[i]) * CDF(g_fft[i]) for i in range(len(f_fft))]

EPS = 1e-14
def div_fft(f_fft, g_fft):
    """Division of two polynomials (FFT representation). Pointwise division."""
    assert len(f_fft) == len(g_fft)
    deg = len(f_fft)
    out = []
    for i in range(len(f_fft)):
        denom = CDF(g_fft[i])
        if abs(denom) < EPS:  # check for zero-ish denominator
            print(f"⚠️ Zero-ish denominator at index {i}: {denom}")
            # optionally: print numerator too
            print(f"Numerator at index {i}: {f_fft[i]}")
        out.append(CDF(f_fft[i]) / denom)
    return out
    #return [CDF(f_fft[i]) / CDF(g_fft[i]) for i in range(len(f_fft))]


def adj_fft(f_fft):
    """Ajoint of a polynomial (FFT representation). Complex conjugation"""
    return [CDF(x).conjugate() for x in f_fft]


"""This value is the ratio between:
    - The degree n
    - The number of complex coefficients of the NTT
While here this ratio is 1, it is possible to develop a short NTT such that it is 2.
"""
fft_ratio = 1
