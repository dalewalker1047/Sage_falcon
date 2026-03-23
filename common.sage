"""This file contains methods and objects which are reused through multiple files."""
from sage.all import IntegerModRing, ZZ

"""q is the integer modulus which is used in Falcon."""
q = 12 * 1024 + 1 #12289

#define sage rings
R = IntegerModRing(q)
# Later define polynomial rings like:
# PR.<x> = PolynomialRing(R)

def split(f):
    """Split a polynomial f in two polynomials.

    Args:
        f: a polynomial

    Format: coefficient
    """
    n = len(f)
    f0 = [f[2 * i] for i in range(n // 2)]
    f1 = [f[2 * i + 1] for i in range(n // 2)]
    return [f0, f1]

def merge(f_list):
    """Merge two polynomials into a single polynomial f.

    Args:
        f_list: a list of polynomials

    Format: coefficient
    """
    f0, f1 = f_list
    n = 2 * len(f0)
    f = [R(0)] * n

    for i in range(n // 2):
        f[2 * i]     = f0[i]
        f[2 * i + 1] = f1[i]

    return f


def sqnorm(v):
    """Compute the square euclidean norm of the vector v."""
    res = ZZ(0)
    for elt in v:
        for coef in elt:
            c = ZZ(coef)  # lift out of mod ring if needed
            res += c * c
    return res

#helper function 
def safe_int(x):
        # if x is a Sage RealDoubleElement, convert via float
        if hasattr(x, 'real'):
            return int(round(float(x.real())))
        else:
            return int(x)