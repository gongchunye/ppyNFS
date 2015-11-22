import primemath
import poly
import math
import etcmath
	
def findDivisors(n):
	if(n < 0):
		raise AssertionError()
	if(n == 1):
		return [1]
	if(n == 2):
		return [1,2]
	if(primemath.isPrime(n)):
		return [1,n]
		
	probPrimes = primemath.generatePrimes(etcmath.isqrt(n))
	primes = []
	for p in probPrimes:
		if(n%p==0):
			primes.append(p)
	
	primeExp = [0]*len(primes)
	for p in primes:
		while(n%p==0):
			n /= p
			primeExp[primes.index(p)] += 1
	
	primeExpC = [0]*len(primes)
	
	divisors = [1]
	while(True):
		i = 0
		while(True):
			primeExpC[i] += 1
			if(primeExpC[i] <= primeExp[i]):
				break
			primeExpC[i] = 0
			i += 1
			if(i >= len(primes)):
				return divisors
		
		d = 1
		for p in primes:
			d *= p**primeExpC[primes.index(p)]
			
		divisors.append(d)
	
def shrinkNfsPoly(poly,m):
	# bound coeffs from -m/2 to m/2.
	for i in range(poly.degree()-1):
		coeff = poly.coeffs[i]
		if coeff > m/2:
			poly.coeffs[i] = -(m - coeff)
			poly.coeffs[i+1] = poly.coeffs[i+1]+1
			
def generateNFSPoly(n,d):
	m = int(math.pow(n,1.0/d))
	nfsPoly = poly.Poly(expansionBaseM(n,m))
	shrinkNfsPoly(nfsPoly,m)
	while(reducible(nfsPoly)):
		m -= 1
		nfsPoly = Poly(expansionBaseM(n,m))
		shrinkNfsPoly(nfsPoly,m)
		
	return (m,nfsPoly)
	
def reducible(poly):
	possibleRoots = findDivisors(abs(poly.coeffs[0]))
	for root in possibleRoots:
		if(poly.evaluate(root) == 0 or poly.evaluate(-root) == 0):
			return True
	return False

def expansionBaseM(n,m):
    q = n
    k = 0
    a = []
    i = len(str(n))

    while q != 0:
        a.append(q % m)
        q = q / m
        k += 1

    return a