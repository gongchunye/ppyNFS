import poly
import etcmath
import files
import primemath
import fractions
import math
import fractions
import sys
		
def getPrimeExponents(value,base):
	primeExponents = [0]*len(base)
	for prime in base:
		baseIndex = base.index(prime)
		while(value % prime == 0):
			primeExponents[baseIndex] += 1
			value /= prime
	
	if(abs(value) != 1):
		print value
		raise AssertionError
	return primeExponents
	
def sumExponents(e1,e2):
#	if(len(e1) != len(e2)):
#		raise AssertionError
	for i in range(len(e1)):
		e1[i] += e2[i]

def getRatPrimeExponents(smooths,base):
	primeExponents = [0]*len(base)
	for smooth in smooths:
		smoothRat = smooth[0]+smooth[1]*m
		smoothPrimeExponents = getPrimeExponents(smoothRat,base)
		sumExponents(primeExponents,smoothPrimeExponents)
		
	return primeExponents
	
def getAlgPrimeExponents(smooths,NF,base):
	primeExponents = [0]*len(base)
	for i in range(len(smooths)):
		smoothNorm = NF(poly.Poly(smooths[i])).norm()
		smoothPrimeExponents = getPrimeExponents(smoothNorm,base)
		sumExponents(primeExponents,smoothPrimeExponents)
		
	return primeExponents
	
def getSqrtModNFromPrimeExponents(primeExponents,base,n):
	prod = 1
	for i in range(len(primeExponents)):
#		if(primeExponents[i] % 2 != 0):
#			raise AssertionError
		#primePower = base[i]**(primeExponents[i]/2) % n
		primePower = pow(base[i],primeExponents[i]/2,n)
		prod = (prod * primePower) % n
		
	return prod
		
def generatePrimes(primes,nfsPoly,couvBound):		
	prodPrimes = 1
	sumBound = 0.0
	for prime in primes:
		sumBound += math.log(prime,2)
		prodPrimes *= prime
		
	while(sumBound <= couvBound):
		prime = primemath.generateLargePrime(64)
		if(not(poly.irreducibleModP(nfsPoly,prime))):
			continue
			
		sumBound += math.log(prime,2)
		primes.append(prime)
		prodPrimes *= prime
		print "%s/%s primes" % (int(sumBound/64)+1,int(couvBound/64)+1)
	return (primes,prodPrimes)
		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	NF = poly.NumberField(poly.Poly(nfsPoly))
	nfsPolyDerivative = poly.Poly(nfsPoly).derivative()
	
	rfBase = files.loadFileArray("rfbase-fil.txt")
	rfBase = zip(*rfBase)[1] #grab only the primes
	afBase = files.loadFileArray("afbase-fil.txt")
	afBase = zip(*afBase)[1]
	afBase = list(set(afBase)) # remove duplicates primes
	smooths = files.loadFileArray("smooths-fil.txt")
	deps = files.loadFileArray("deps.txt")
	
	
	tryAllDeps = False
	if(len(sys.argv) == 2 and sys.argv[1] == "-a"):
		tryAllDeps = True
		print "Trying all dependencies..."
	
	noncongruent = 0
	success = 0
	primes = []
	for dependency in deps:
		print "dependency %s/%s" % (deps.index(dependency)+1,len(deps))
		dependencySmooths = []
		ctr = 0
		for column in dependency:
			if(column == 1):
				smooth = smooths[ctr]
				dependencySmooths.append(smooth)
			ctr += 1
		
		if(len(dependencySmooths) == 0):
			continue
		
		print "dependency has %s entries." % len(dependencySmooths) 
		
		primeExponents = getRatPrimeExponents(dependencySmooths,rfBase)
		ratSide = getSqrtModNFromPrimeExponents(primeExponents,rfBase,n)
		ratSide = (ratSide * nfsPolyDerivative.evaluate(m)) % n
		
		couvBound = etcmath.calcRequiredPrimeLength(n,m,poly.Poly(nfsPoly),dependencySmooths)
		(primes,prodPrimes) = generatePrimes(primes,poly.Poly(nfsPoly),couvBound)
		primeExponents = getAlgPrimeExponents(dependencySmooths,NF,afBase)
		
		sum = 0
		for prime in primes:
			print "%s/%s sqrt" % (primes.index(prime)+1,len(primes))
			NFp = poly.NumberFieldModP(poly.Poly(nfsPoly),prime) 
			
			normModP = getSqrtModNFromPrimeExponents(primeExponents,afBase,prime)
			normModP = (normModP*NFp(nfsPolyDerivative).norm()) % prime
			
			prod = NFp(poly.Poly([1]))
			for smooth in dependencySmooths:
				prod = prod * NFp(poly.Poly(smooth))
				
			prod = prod * (NFp(nfsPolyDerivative)**2)
				
			sqrt = prod.sqrt()	
			if(sqrt.norm() != normModP):
				sqrt = -sqrt
			if(sqrt.norm() != normModP):				
				raise AssertionError
				
			q = prodPrimes/prime
			x = etcmath.modinv(q,prime)
			a = sqrt.getPoly().evaluate(m % prime) % prime
			sum = (sum + a*x*q) % prodPrimes
			
		if(math.log(sum) > (math.log(prodPrimes))/2):
			sum = -(-sum % prodPrimes)		
			
		algSide = sum % n
		
		if(algSide**2 % n != ratSide**2 % n):
			noncongruent += 1
			print "Warning: not congruent."
		
		possibleFactor = fractions.gcd(n,abs(algSide-ratSide))
		print "%s = %s*p" % (n,possibleFactor)
		if(possibleFactor != 1 and possibleFactor != n):
			success += 1
			if(not(tryAllDeps)):
				break
		else:
			print "Trivial factor found, trying next dependency..."
	
	print "%s/%s successful factorizations." % (success,len(deps))
	print "%s non-congruences." % noncongruent