# Test harness for Falcon KAT tests
from hashlib import shake_256
from sage.all import ZZ, PolynomialRing, GF, vector
from ast import literal_eval
import json
import sys
import traceback  # Import traceback for detailed error reporting
from sign_KAT import sign_KAT2


load('falcon_sign.sage')
load('ffsampling.sage')
load('common.sage')
load('fft.sage')
load('samplerz.sage')
load('encoding.sage')
load('ntt.sage')


def load_kat_vectors_from_sign_KAT():
    """Load KAT vectors from the sign_KAT.py file."""
    kat_vectors = []
    for vector in sign_KAT2:
        kat_vectors.append({
            "count": vector.get("read_bytes"),
            "f": vector.get("f"),
            "g": vector.get("g"),
            "F": vector.get("F"),
            "G": vector.get("G"),
            "nonce": vector.get("nonce"),
            "prng_seed": vector.get("prng_seed"),
            "s2": vector.get("s2"),
            "sig": vector.get("sig"),
        })
    return kat_vectors

def verify_signature(message, signature, public_key):
    """Placeholder for signature verification logic."""
    # Implement the verification logic here
    pass

def test_kat_without_sk(count, seed, mlen, msg, pk, sk, smlen, sm):
    if not sk:  # If the secret key is empty
        print(f"Testing vector {count} with empty secret key.")
        # Perform checks using only the public key and other parameters
        # Example: Verify the signature using the public key
        try:
            is_valid = verify_signature(msg, sm, pk)  # Assuming `verify` is defined elsewhere
            if not is_valid:
                print(f"Test vector {count} failed verification with empty secret key.")
                return False
        except Exception as e:
            print(f"Error during verification for vector {count}: {e}")
            return False

    # Generate signature
    n = 512
    q = 12289
    beta_sq = 1.36  # Example value, adjust as needed
    slen = 256  # Example value, adjust as needed

    generated_sig = sign(sk, msg, seed)

    # Compare generated signature with expected signature
    if generated_sig == sm:
        print(f"KAT vector {count} passed.")
        return True
    else:
        print(f"KAT vector {count} failed.")
        print(f"Expected: {sm}")
        print(f"Generated: {generated_sig}")
        return False

def test_kat(count, seed, mlen, msg, pk, sk, smlen, sm):
    if not sk:
        return test_kat_without_sk(count, seed, mlen, msg, pk, sk, smlen, sm)

    # Generate signature
    n = 512
    q = 12289
    beta_sq = 1.36  # Example value, adjust as needed
    slen = 256  # Example value, adjust as needed

   
    #put list of sk into key gen. put that into sign job done 
    generated_sig = sign(sk, msg, seed)

    # Compare generated signature with expected signature
    if generated_sig == sm:
        print(f"KAT vector {count} passed.")
        return True
    else:
        print(f"KAT vector {count} failed.")
        print(f"Expected: {sm}")
        print(f"Generated: {generated_sig}")
        return False




if __name__ == "__main__":
    # Load KAT vectors from sign_KAT.py
    kat_vectors = load_kat_vectors_from_sign_KAT()
    
    # Run tests for each KAT vector
    for vector in kat_vectors:
        count = vector.get("count")
        f = vector.get("f")
        g = vector.get("g")
        F = vector.get("F")
        G = vector.get("G")
        nonce = vector.get("nonce")
        prng_seed = vector.get("prng_seed")
        s2 = vector.get("s2")
        sig = vector.get("sig")

        #raw_sk = (f, g, F, G)  # Combine f, g, F, G into the secret key tuple
        f = pad_to_n(f)
        g = pad_to_n(g)
        F = pad_to_n(F)
        G = pad_to_n(G)

        sk, vk = keygen([f, g, F, G])
        msg = b"data1"
        
        sm = sign(sk, msg)



        print(f"\nRunning test vector {count}...")
        print(f"f: {f}")
        print(f"g: {g}")
        print(f"F: {F}")
        print(f"G: {G}")
        print(f"Nonce: {nonce}")
        print(f"PRNG Seed: {prng_seed}")
        print(f"s2: {s2}")
        print(f"Signature: {sig}")

        # Run the test
        try:
            result = test_kat(count, prng_seed, len(msg), msg, None, sk, len(sm), sm)
            if result:
                print(f"Test vector {count} passed.")
            else:
                print(f"Test vector {count} failed.")
        except Exception as e:
            print(f"Error while testing vector {count}: {e}")
            traceback.print_exc()  # Print the full exception stack trace