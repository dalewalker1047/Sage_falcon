
def compress(v, slen):
    bits = []

    for coeff in v: 
        # Sign bit
        bits.append(1 if coeff < 0 else 0)
        abs_coeff = abs(coeff)

        # 7 lowest bits
        for i in range(7):
            bits.append((abs_coeff >> i) & 1)
        
        # high bits get unary encoding
        high = abs_coeff >> 7
        bits.extend([0]*high + [1])

    # Check length
    if len(bits) > slen*8:
        raise ValueError("Too many bits to fit in signature")
    
    # Pad with zeros
    bits.extend([0]*(slen*8 - len(bits)))

    # Pack bits into bytes
    byte_list = []
    for i in range (0, len(bits), 8):
        byte = 0
        for j in range(8): 
            byte |= (bits[i+j] << j)
        byte_list.append(byte)

    return bytes(byte_list)


def decompress(x, slen, n):
    if len(x) > slen:
        raise ValueError("Signature too long")
    
    # Unpack bytes into bits
    bits = []
    for b in x: 
        b_int = Integer(b)
        bits.extend([b_int >> i & 1 for i in range(8)])

    
    # Remove padding
    while bits and bits[-1] == 0: 
        bits.pop()

    v = []
    i = 0
    try: 
        while i < len(bits) and len(v) < n:
            sign = -1 if bits[i] == 1 else 1
            i += 1

            # 7 lowest bits
            low = sum(bits[i+j] << j for j in range(7))
            i += 7

            # High bits
            high = 0

            while bits[i] == 0: 
                high += 1
                i += 1
            i += 1 # skip the terminating 1

            coeff = sign * (low + (high << 7))
            if coeff == 0 and sign == -1: 
                return False
            v.append(coeff)

        if len(v) != n:
            return False
        
        return vector(ZZ, v)
    
    except IndexError: 
        return False

