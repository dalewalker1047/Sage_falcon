# Importing dependencies
from math import floor

# Use high-quality randomness
# The "secrets" library could also work (Python >= 3.6)
from os import urandom


# Upper bound on all the values of sigma
MAX_SIGMA = 1.8205
INV_2SIGMA2 = 1 / (2 * (MAX_SIGMA ** 2))

# Precision of RCDT
RCDT_PREC = 72

# ln(2) and 1 / ln(2), with ln the natural logarithm
LN2 = 0.69314718056
ILN2 = 1.44269504089


# RCDT is the reverse cumulative distribution table of a distribution that
# is very close to a half-Gaussian of parameter MAX_SIGMA.
RCDT = [
    3024686241123004913666,
    1564742784480091954050,
    636254429462080897535,
    199560484645026482916,
    47667343854657281903,
    8595902006365044063,
    1163297957344668388,
    117656387352093658,
    8867391802663976,
    496969357462633,
    20680885154299,
    638331848991,
    14602316184,
    247426747,
    3104126,
    28824,
    198,
    1]


# C contains the coefficients of a polynomial that approximates exp(-x)
# More precisely, the value:
# (2 ** -63) * sum(C[12 - i] * (x ** i) for i in range(i))
# Should be very close to exp(-x).
# This polynomial is lifted from FACCT: https://doi.org/10.1109/TC.2019.2940949
C = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000]


def basesampler(randombytes=urandom):
    #might want randombytes = ChaCha20.randombytes is signatures don't match reference
    """
    Sample z0 in {0, 1, ..., 18} with a distribution
    very close to the half-Gaussian D_{Z+, 0, MAX_SIGMA}.
    Takes as (optional) input the randomness source (default: urandom).
    """
    u = int.from_bytes(randombytes(RCDT_PREC >> 3), "little")

    z0 = 0
    for elt in RCDT:
        z0 += int(u < elt)
    return z0


def approxexp(x, ccs):
    """
    Compute an approximation of 2^63 * ccs * exp(-x).

    Input:
    - a floating-point number x
    - a scaling factor ccs
    Both inputs x and ccs MUST be positive.

    Output:
    - an integral approximation of 2^63 * ccs * exp(-x).
    """
    # y, z are always positive
    y = C[0]
    # Since z is positive, int is equivalent to floor
    z = int(x * (1 << 63))
    for elt in C[1:]:
        y = elt - ((z * y) >> 63)
    z = int(ccs * (1 << 63)) << 1
    y = (z * y) >> 63
    return y


def berexp(x, ccs, randombytes=urandom):
    #might want randombytes = ChaCha20.randombytes is signatures don't match reference
    
    """
    Return a single bit, equal to 1 with probability ~ ccs * exp(-x).
    Both inputs x and ccs MUST be positive.
    Also takes as (optional) input the randomness source (default: urandom).
    """
    s = int(x * ILN2)
    r = x - s * LN2
    s = min(s, 63)
    z = (approxexp(r, ccs) - 1) >> s
    w = 0

    for i in range(56, -8, -8):
        p = int.from_bytes(randombytes(1), "little")
        w = p - ((z >> i) & 0xFF)
        if w:
            break
    return (w < 0)

    #rng = ChaCha20(seed)
    #samplerz(mu, sigma, sigmin, randombytes=rng.randombytes)
def samplerz(mu, sigma, sigmin, randombytes=urandom):
    #might want randombytes = ChaCha20.randombytes is signatures don't match reference
    """
    Given floating-point values mu, sigma (and sigmin),
    output an integer z according to the discrete
    Gaussian distribution D_{Z, mu, sigma}.

    Input:
    - the center mu
    - the standard deviation sigma
    - a scaling factor sigmin
    - optional: the randomness source randombytes (default: urandom)
      randombytes(k) should output k pseudorandom bytes
    The inputs MUST verify 1 < sigmin < sigma < MAX_SIGMA.

    Output:
    - a sample z from the distribution D_{Z, mu, sigma}.
    """
    

    #Force floats explicitly:
    sigma = float(sigma)
    mu = float(mu)

    s = int(floor(mu))
    r = mu - s
    dss = 1 / (2 * sigma * sigma)
    ccs = sigmin / sigma

    #go back to while(1) and no raise runtime error after debugging  
    for _ in range(100000):
        # Sampler z0 from a Half-Gaussian
        z0 = basesampler(randombytes=randombytes)
        # Convert z0 into a pseudo-Gaussian sample z
        b = int.from_bytes(randombytes(1), "little")
        b &= 1
        z = b + (2 * b - 1) * z0
        # Rejection sampling to obtain a true Gaussian sample
        x = ((z - r) ** 2) * dss
        x -= (z0 ** 2) * INV_2SIGMA2
        if berexp(x, ccs, randombytes=randombytes):
            return z + s
    raise RuntimeError("Sampler stuck")

"""
KATs TEST
"""

from samplerz_KAT512 import sampler_KAT512
from samplerz_KAT1024 import sampler_KAT1024

def KAT_randbytes(k):
    
    #Use a fixed bytestring 'octets' as a source of random bytes
    
    global octets
    oc = octets[: (2 * k)]
    if len(oc) != (2 * k):
        raise IndexError("Randomness string out of bounds")
    octets = octets[(2 * k):]
    return bytes.fromhex(oc)[::-1]


"""
debugging test
"""
"""
from sage.all import mean, variance
if __name__ == "__main__":
    import statistics
    import matplotlib.pyplot as plt

    # Parameters for the test
    mu = 0.0       # center
    sigma = 1.5    # standard deviation
    sigmin = 1.0   # minimal sigma
    nsamples = 10000

    # Use system randomness
    from os import urandom
    random_source = urandom

    # Generate samples
    samples = [samplerz(mu, sigma, sigmin, random_source) for _ in range(nsamples)]

    # Basic checks
    for s in samples:
        if not isinstance(s, int):
            print("Non-integer detected:", s, type(s))
    #assert all(isinstance(s, int) for s in samples), "Non-integer output detected!"

    # Statistics
    means = mean(samples)
    var = variance(samples)
    print(f"Sampler test with mu={mu}, sigma={sigma}")
    print(f"Mean = {means}, Expected ~ {mu}")
    print(f"Variance = {var}, Expected ~ {sigma**2}")
    print(f"Min = {min(samples)}, Max = {max(samples)}")

    # Histogram for visual inspection
    plt.hist(samples, bins=range(min(samples), max(samples)+2), density=True, alpha=0.7, color='skyblue')
    plt.title("Histogram of samplerz output")
    plt.xlabel("Sample value")
    plt.ylabel("Frequency")
    plt.show()

    # Check empirical 6-sigma range
    count_outside = sum(1 for s in samples if abs(s - mu) > 6*sigma)
    print(f"Samples outside ±6σ: {count_outside} / {nsamples}")    
"""