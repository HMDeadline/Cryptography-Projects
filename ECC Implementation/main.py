class Elliptic_Curve:

    # None will be considered as the point in infinity, i.e. (0;1;0)
    # Also note that for simplicity of calculation, the curve will be considered in a finite field Fp
    # Curve is of the form Y^2 = X^3 + AX + B

    def __init__(self, A, B, prime):
        self.A = A
        self.B = B
        self.p = prime
        if ((4 * A ** 3 + 27 * B ** 2) % prime == 0):
            raise ValueError("Curve is not smooth! (Curve has points of singularity)")
        
    def satisfies_equation(self, point): # Check if a point lies on the curve or not
        if point is None: return True
        else:
            X , Y = point
            return (Y ** 2 - X ** 3 - self.A * X - self.B) % self.p == 0

    def inverse(self, point): # Calculates the additive inverse of a point, i.e. -P
        if point is None: return point
        else:
            X , Y = point
            return (X, (-Y) % self.p)

    def addition(self, p_1, p_2): # Function to calculate the addition of two points assumed to be lying on the curve
        if p_1 is None: return p_2
        elif p_2 is None: return p_1
        elif p_1 == self.inverse(p_2): return None
        elif p_1 == p_2: # In the case of calculating 2P we use the slope of the tangent line and find the x-coordinate of the point of addition
            X , Y = p_1
            partial_x, partial_y = (3 * X ** 2 + self.A) % self.p, (2 * Y) % self.p
            slope = (partial_x * pow(partial_y, -1, self.p)) % self.p
            add_x = (slope ** 2 - 2 * X) % self.p

        else: # Otherwise we use the usual slope of the line connecting the two points
            (X, Y), (X_, Y_) = (p_1, p_2)
            slope = ((Y-Y_)*pow(X-X_,-1,self.p)) % self.p
            add_x = (slope ** 2 - X - X_) % self.p
        
        return (add_x, (slope * (X - add_x) - Y) % self.p)

    def scalar_multiplication(self, k, point): # We use the efficient method described in Silverman's book section 6.3.1, i.e. using the base 2 expansion of k
        ans = None
        p = point
        while k:
            if k & 1:
                ans = self.addition(ans, p)
            p = self.addition(p, p)
            k >>= 1
        return ans
    

# As an example we use the curve in Example 6.15 of Silverman's book, i.e. Y^2 = X^3 + 8X + 7 over F_73
Curve = Elliptic_Curve(8, 7, 73)
P = (32, 53)
assert Curve.satisfies_equation(P)

# Private keys for Alice and Bob
nA = 37
nB = 28

# Public keys

QA = Curve.scalar_multiplication(nA, P) # Cross check with the book that the result is (35,47)
QB = Curve.scalar_multiplication(nB, P) # Cross check with the book that the result is (58,4)

print("Alice and Bob's private keys are the following, in order:", QA, QB, sep='\n')

# Finally calculate the shared secret they arrive at:
SA = Curve.scalar_multiplication(nA, QB)
SB = Curve.scalar_multiplication(nB, QA)

print("Alice and Bob both arrive at the following shared secrets, in order:", SA, SB, "which verifies the algorithm works!", sep='\n')