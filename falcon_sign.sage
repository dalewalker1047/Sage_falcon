from sage.all import *
from sage.all import vector as sage_vector
from hashlib import shake_256
import os
from os import urandom
from typing import List, Optional
from math import sqrt 

load('ffsampling.sage')
load('common.sage')
load('fft.sage')
load('samplerz.sage')
load('encoding.sage')
load('ntt.sage')
load('ntrugen.sage')

n = 512
q = 12289

Zx = PolynomialRing(ZZ, 'x')  # Define the integer polynomial ring
x = Zx.gen()
Fq = GF(q)  # Define the finite field
Fqx = PolynomialRing(Fq, 'x')  # Define the polynomial ring over Fq
Rq = Fqx.quotient(x**n + 1)  # Define the quotient ring

logn = {
    512: 9,
    1024: 10
}

sigma_min = {
    512: 1.277833697,
    1024: 1.298280334,
}

sigma = {
    512: 165.7366171829776,
    1024: 168.38857144654395,
}

HEAD_LEN = 1
SALT_LEN = 40
SEED_LEN = 56


def hash_to_point(message, nonce, n, q):

    hasher = shake_256()
    hasher.update(nonce)    
    hasher.update(message)
    stream = hasher.digest(int(2*n*10))
    i = 0

    coeffs = []
    limit = (1 << 16) - ((1 << 16) % q)

    while len(coeffs) < n:
        val = int.from_bytes(stream[i:i+2], 'little')
        i += 2

        if val < limit: 
            coeffs.append(val % q)

    return sage_vector(ZZ, coeffs)


def poly_norm_sq(f):
    lifted = f.lift()
    coeffs = lifted.list()
    return sum(int(c)**2 for c in coeffs)


def parse_secret_key(sk, n):
    #print("Debug parse_secret_key start")
    """
    Parse the secret key from a tuple of lists into sage_vectors, padding components to length n if necessary.

    Parameters:
        sk (tuple): The secret key as a tuple (f, g, F, G).
        n (int): Degree of the polynomial ring.

    Returns:
        tuple: A tuple (f, g, F, G) of sage_vectors.
    """
    if not isinstance(sk, tuple) or len(sk) != 4:
        raise ValueError("Invalid secret key format. Expected a tuple (f, g, F, G).")

    f, g, F, G = sk


    """
    print("Debug: Type of pad_to_n:", type(pad_to_n))
    print("Debug: Value of n:", n)
    print("Debug: Type and value of f:", type(f), f)
    print("Debug: Type and value of g:", type(g), g)
    print("Debug: Type and value of F:", type(F), F)
    print("Debug: Type and value of G:", type(G), G)
    """

    f = pad_to_n(f)
    g = pad_to_n(g)
    F = pad_to_n(F)
    G = pad_to_n(G)
    #print("Debug parse_secret_key end")
    #print("Type of sage_vector:", type(sage_vector))
    #print("Type of ZZ:", type(ZZ))
    return sage_vector(ZZ, f), sage_vector(ZZ, g), sage_vector(ZZ, F), sage_vector(ZZ, G)

def pad_to_n(poly):
    if len(poly) > n:
        raise ValueError("A component of the secret key is longer than n.")
    return poly + [0] * (n - len(poly))

def sign(sk, message, randombytes=urandom):

    f, g, F, G, B0_fft, T_fft = sk

    int_header = 0x30 + logn[n]
    header = int(int_header).to_bytes(1, 'little')

    salt = randombytes(SALT_LEN)

    hashed = hash_to_point(message, salt, n, q)
    
    while True:
        #print("T slice: ", T_fft[-2:])
        print(type(T_fft))
        print(len(T_fft))
        
        if randombytes == urandom:
            #print("Hashed: ", hashed)
            #print(T_fft)
            s = ffsampling_fft([hashed, hashed], T_fft, 1.2, randombytes)
        else:
            seed = randombytes(SEED_LEN)
            s = ffsampling_fft([hashed, hashed], T_fft, 1.2, seed)
        print("debug: got to norm sign")
        norm_sign = sum(Integer(coeff)**2 for coeff in s[0])
        norm_sign += sum(Integer(coeff)**2 for coeff in s[1])

        if norm_sign <= 1.36 * n:
            enc_s = compress(s[1], sig_bytelen - SALT_LEN - HEAD_LEN, n)

            if enc_s is not False:
                return header + salt + enc_s

def verify(vk, message, signature, n, q, sig_bound, sig_bytelen, HEAD_LEN, SALT_LEN):

    if (8*len(vk) % n != 0):
        raise ValueError("Invalid public key length")

    # Polynomial rings
    Zx.<x> = ZZ[]
    Rq.<x> = GF(q)['x'].quotient(x**n + 1)

    # Deserialize public key
    h_list = deserialize_to_poly(vk, n)
    h = Rq(h_list)

    # Get signature parts
    salt = signature[HEAD_LEN:HEAD_LEN+SALT_LEN]
    enc_s = signature[HEAD_LEN+SALT_LEN:]

    s1 = decompress(enc_s, sig_bytelen - SALT_LEN - HEAD_LEN, n)

    if s1 is False:
        print("Invalid encoding")
        return False

    s1 = sage_vector(ZZ, s1)

    # Hash to point
    hashed = hash_to_point(message, salt, n, q)
    c = Rq(hashed)

    s1_poly = Rq(s1.list())
    s0_poly = c - s1_poly * h

    # Lift
    s0_coeffs = []
    for coeff in s0_poly.lift().list():
        centered = (int(coeff) + (q>>1)) % q - (q>>1)
        s0+coeffs.append(centered)

    s0 = sage_vector(ZZ, s0_coeffs)

    norm_sq = sum(c*c for c in s0) + sum(c*c for c in s1)
    if norm_sq > sig_bound:
        print("Signature norm too large")
        return False

    return True

def deserialize_to_poly(bytestring, n):
    assert (8*len(bytestring) % n == 0)

    BITS_PER_COEFF = 14
    mask = (1 << BITS_PER_COEFF) - 1

    int_buffer = Integer(int.from_bytes(bytestring, 'little'))

    coeffs = []

    for _ in range(n):
        coeffs.append(int(int_buffer & mask))
        int_buffer >>= BITS_PER_COEFF

    return sage_vector(ZZ, coeffs)


def verify_signature(message, signature, public_key):
    """Wrapper for the verify function with predefined parameters."""
    n = 512
    q = 12289
    sig_bound = 1.36 * n 
    sig_bytelen = 666  
    HEAD_LEN = 40  
    SALT_LEN = 40  

    return verify(public_key, message, signature, n, q, sig_bound, sig_bytelen, HEAD_LEN, SALT_LEN)

def keygen(polys: Optional[List[List[int]]]=None):
    """
    Initialize a secret key.
    """
    # Public parameters

    # Compute NTRU polynomials f, g, F, G verifying fG - gF = q mod Phi
    if polys is None:
        f, g, F, G = ntru_gen(n)
    else:
        
        assert all((len(poly) == n) for poly in polys)
        [f, g, F, G] = [poly[:] for poly in polys]

    # From f, g, F, G, compute the basis B0 of a NTRU lattice
    # as well as its Gram matrix and their fft's.
    B0 = [[g, neg(f)], [G, neg(F)]]
    G0 = gram(B0)
    B0_fft = [[fft(elt) for elt in row] for row in B0]
    G0_fft = [[fft(elt) for elt in row] for row in G0]

    T_fft = ffldl_fft(G0_fft)

    #sigma=165.7366171829776
    # for 512 sigma=165.7366171829776
    # for 1024 sigma=168.38857144654395
    # Normalize Falcon tree
    normalize_tree(T_fft, sigma[n])

    # The public key is a polynomial such that h*f = g mod (Phi,q)
    h = div_zq(g, f)
    


    sk = (f, g, F, G, B0_fft, T_fft)
    #convert h into a python int before passing to serialize
    h = [safe_int(x) for x in h]
    
    vk = serialize_poly(h)
    return (sk, vk)

def serialize_poly(poly: List[int]) -> bytes:
    """
    Serialize a polynomial (a list of integers) into a bytestring.
    We assume that all entries are between 0 and q - 1.
    """
    n = len(poly)

    if (min(poly) < 0) or (max(poly) >= q):
        raise ValueError("The entries of poly are outside bounds")

    BITS_PER_COEF = int(14)
    int_buffer = int(0)
    for idx in range(n):
        shift = idx * BITS_PER_COEF
        int_buffer |= poly[idx] << shift
    # The "+ 7" allows to round to the nearest higher integer.
    bytelen = (n * BITS_PER_COEF + 7) >> 3

    return int_buffer.to_bytes(bytelen, 'little')

def normalize_tree(tree, sigma):
    """
    Normalize leaves of a LDL tree (from values ||b_i||**2 to sigma/||b_i||).

    Args:
        T: a LDL tree
        sigma: a standard deviation

    Format: coefficient or fft
    """
    if len(tree) == 3:
        normalize_tree(tree[1], sigma)
        normalize_tree(tree[2], sigma)
    else:
        val = clean_real(tree[0])
        #print(val)
        tree[0] = sigma / sqrt(val)
        tree[1] = 0

def pack_sk(self, sk) -> bytes:
        """
        Pack a Falcon secret key into a bytestring.
        To be used in conjunction with unpack_sk.
        
        Not used in this module but can be useful for exporting keys.
        """
        (f, g, F, G, B0_fft, T_fft) = sk
        sk_bytes = b""
        for poly in (f, g, F, G):
            # We reduce modulo q before serializing
            sk_bytes += serialize_poly([coef % q for coef in poly])
        return sk_bytes

def unpack_sk(self, sk_bytes: bytes):
    """
    Unpack a bytestring to a Falcon secret key.
    
    How to use:

    >>> falcon = Falcon(n)
    >>> sk, _ = falcon.keygen()
    >>> sk_bytes = falcon.pack(sk)
    >>> sk_2 = falcon.unpack_sk(sk_bytes)
    >>> assert (sk == sk2)
    
    """
    assert (len(sk_bytes) % 4 == 0)
    # There are four polys in sk
    len_poly = len(sk_bytes) // 4

    polys = [None, None, None, None]
    for i in range(4):
        polys[i] = deserialize_to_poly(sk_bytes[i * len_poly:(i + 1) * len_poly], self.param.n)
        polys[i] = [((coef + (q >> 1)) % q) - (q >> 1) for coef in polys[i]]
    sk, _ = self.keygen(polys)
    return sk