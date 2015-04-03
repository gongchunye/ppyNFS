from poly import *
from etc_math import *
from files import *
from prime import *
import fractions
from fractions import *
		
def getPrimeExponents(value,base):
	primeExponents = [0]*len(base)
	for prime in base:
		baseIndex = base.index(prime)
		while(value % prime == 0):
			primeExponents[baseIndex] += 1
			value /= prime
			
	return primeExponents
	
def sumExponents(e1,e2):
#	if(len(e1) != len(e2)):
#		raise AssertionError
	for i in range(len(e1)):
		e1[i] += e2[i]

def getRatPrimeExponents(smooths,base):
	primeExponents = [0]*len(base)
	for i in range(len(smooths)):
		smoothRat = Poly(smooths[i]).evaluate(m)
		smoothPrimeExponents = getPrimeExponents(smoothRat,base)
		sumExponents(primeExponents,smoothPrimeExponents)
		
	return primeExponents
	
def getNormPrimeExponents(smooths,NF,base):
	primeExponents = [0]*len(base)
	for i in range(len(smooths)):
		smoothNorm = NF(Poly(smooths[i])).norm()
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
		prime = generateLargePrime(64)
		if(not(irreducibleModP(nfsPoly,prime))):
			continue
			
		sumBound += math.log(prime,2)
		primes.append(prime)
		prodPrimes *= prime
		print "%s/%s primes" % (int(sumBound/64)+1,int(couvBound/64)+1)
	return (primes,prodPrimes)
		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	nfsPolyDerivative = Poly(nfsPoly).derivative()
	
	rfBase = loadFileArray("rfbase.txt")
	smooths = loadFileArray("smooths.txt")
	deps = loadFileArray("deps.txt")
	
	primes = []
	for dependency in deps:
		dependencySmooths = []
		ctr = 0
		for column in dependency:
			if(column == 1):
				smooth = smooths[ctr]
				dependencySmooths.append(smooth)
			ctr += 1
			
		primeExponents = getRatPrimeExponents(dependencySmooths,rfBase)
		ratSide = getSqrtModNFromPrimeExponents(primeExponents,rfBase,n)
		ratSide = (ratSide * nfsPolyDerivative.evaluate(m)) % n
		
		couvBound = calcRequiredPrimeLength(n,m,Poly(nfsPoly),dependencySmooths)
		(primes,prodPrimes) = generatePrimes(primes,Poly(nfsPoly),couvBound)
		primeExponents = getNormPrimeExponents(dependencySmooths,NF,rfBase)
		
		sum = 0
		for prime in primes:
			print "%s/%s sqrt" % (primes.index(prime)+1,len(primes))
			NFp = NumberFieldModP(Poly(nfsPoly),prime) 
			
			normModP = getSqrtModNFromPrimeExponents(primeExponents,rfBase,prime)
			normModP = (normModP*NFp(nfsPolyDerivative).norm()) % prime
			
			prod = NFp(Poly([1]))
			for smooth in dependencySmooths:
				prod = prod * NFp(Poly(smooth))
				
			prod = prod * (NFp(nfsPolyDerivative)**2)
				
			sqrt = prod.sqrt()	
			if(sqrt.norm() != normModP):
				sqrt = -sqrt
			if(sqrt.norm() != normModP):				
				raise AssertionError
				
			q = prodPrimes/prime
			x = modinv(q,prime)
			a = sqrt.getPoly().evaluate(m % prime) % prime
			sum = (sum + a*x*q) % prodPrimes
			
		if(math.log(sum,10) > (math.log(prodPrimes,10))/2):
			sum = -(-sum % prodPrimes)		
			
		algSide = sum % n
		
		if(algSide**2 % n != ratSide**2 % n):
			print "Warning: not congruent."
		
		possibleFactor = fractions.gcd(n,abs(algSide-ratSide))
		print "%s = %s*p" % (n,possibleFactor)
		if(possibleFactor != 1 and possibleFactor != n):
			break
		else:
			print "Trivial factor found, trying next dependency..."
		