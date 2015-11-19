import math

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modInv(a, m):
	a %= m
	#print "%s,%s" % (a,m)
	if(a == 0):
		raise AssertionError
	g, x, y = egcd(a, m)
	if g != 1:
		raise Exception('modular inverse does not exist')
	else:
		return x % m

def qe_mod(a,b,c,p):
	roots = []
	num = modular_sqrt(b*b - 4*a*c,p)
	den = modInv(2*a,p)
	roots.append((-b+num)*den % p)
	roots.append((-b-num)*den % p)
	return roots
	
def modular_sqrt(a, p):
    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.

        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.

        0 is returned is no square root exists for
        these a and p.

        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases
    #
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) / 4, p)

    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s /= 2
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    x = pow(a, (s + 1) / 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in xrange(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m


def legendre_symbol(a, p):
    """ Compute the Legendre symbol a|p using
        Euler's criterion. p is a prime, a is
        relatively prime to p (if p divides
        a, then a|p = 0)

        Returns 1 if a has a square root modulo
        p, -1 otherwise.
    """
    ls = pow(a, (p - 1) / 2, p)
    return -1 if ls == p - 1 else ls		
	
#http://code.activestate.com/recipes/577821-integer-square-root-function/
def isqrt(x):
    if x < 0:
        raise ValueError('square root not defined for negative numbers')
    n = int(x)
    if n == 0:
        return 0
    a, b = divmod(n.bit_length(), 2)
    x = 2**(a+b)
    while True:
        y = (x + n//x)//2
        if y >= x:
            return x
        x = y	
		

		
	
def trialdivide(todivide, factorbase):
	if(todivide == 0):
		return False
	todivide = abs(todivide)
	
	for p in factorbase:
		if(len(p) == 2):
			p = p[1]
		while(todivide % p == 0):
			todivide = todivide/p
		if(todivide == 1):
			break
			
	if(todivide == 1):
		return True
	else: 
		return False


def jacobi(n, m):
    j = 1
 
    # rule 5
    n %= m
     
    while n:
        # rules 3 and 4
        t = 0
        while not n & 1:
            n /= 2
            t += 1
        if t&1 and m%8 in (3, 5):
            j = -j
 
        # rule 6
        if (n % 4 == m % 4 == 3):
            j = -j
 
        # rules 5 and 6
        n, m = m % n, n
 
    return j if m == 1 else 0
	

def calculateB(n):	
	top=math.pow(8.0/9.0,1.0/3.0)*math.pow(math.log(n),1.0/3.0)*math.pow(math.log(math.log(n)),2.0/3.0)
	return int(math.exp(top))
	
	
def calcRequiredPrimeLength(n, m, nfspoly, tomult):
	logest = math.log(nfspoly.degree(),2)*(nfspoly.degree()+5)*.5;
	logest += math.log(n,2)
	maxu = tomult[0][0]
	for smooth in tomult:
		if(maxu < abs(smooth[0])):
			maxu = abs(smooth[0])
		if(maxu < abs(smooth[1])):
			maxu = abs(smooth[1])
	b = 2*maxu*math.sqrt(nfspoly.degree())*m
	logest+= len(tomult)*.5*math.log(b,2)
	return logest
	
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
		
def getsr(p,e):
	q = p**e - 1
	counter = 0
	while (q%2 == 0):
		counter += 1
		q /= 2
	return (q,counter)
	