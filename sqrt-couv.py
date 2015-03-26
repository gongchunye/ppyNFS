from poly import *
from etc_math import *
from files import *
from sqrt_funcs import *
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
	
def getPrimeExponentsSqrtModN(primeExponents,base,n):
	prod = 1
	for i in range(len(primeExponents)):
		if(primeExponents[i] % 2 != 0):
			raise AssertionError
		primePower = base[i]**(primeExponents[i]/2) % n
		prod = (prod * primePower) % n
		
	return prod
		
def generatePrimes(primes,nfsPoly,sumBound,couvBound):		
	sumBound = 0
	prodPrimes = 1
	while(sumBound <= couvBound):
		prime = generateLargePrime(64)
		
		if(not(irreducibleModP(nfsPoly,prime))):
			continue
			
		(s,r) = getsr(prime,3)
		if(r > 5):
			continue
			
		sumBound += math.log(prime,2)
		primes.append(prime)
		prodPrimes *= prime
		print "%s/%s primes" % (int(sumBound/64)+1,int(couvBound/64)+2)
	return (primes,prodPrimes,sumBound)
		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	
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
		ratSide = getPrimeExponentsSqrtModN(primeExponents,rfBase,n)
		
		couvBound = calcRequiredPrimeLength(n,m,Poly(nfsPoly),dependencySmooths)
		sumBound = 0.0
		(primes,prodPrimes,sumBound) = generatePrimes(primes,Poly(nfsPoly),sumBound,couvBound)
		primeExponents = getNormPrimeExponents(dependencySmooths,NF,rfBase)
		
		r = Fraction()
		sum = 0
		for prime in primes:
			print "%s/%s sqrt" % (primes.index(prime)+1,len(primes))
			normModP = getPrimeExponentsSqrtModN(primeExponents,rfBase,prime)
			
			NFp = NumberFieldModP(Poly(nfsPoly),prime)
			prod = NFp(Poly([1]))
			for smooth in dependencySmooths:
				prod = prod * NFp(Poly(smooth))
				
			sqrt = prod.sqrt()	
			if(sqrt.norm() != normModP):
				sqrt = -sqrt
				
			x = modinv(prodPrimes/prime,prime)
			a = sqrt.getPoly().evaluate(m % prime) % prime
			sum += a*x*prodPrimes/prime
			r += Fraction((a*x),prime)
			
		algSide = (sum - int(r)*prodPrimes) % n
		
		possibleFactor = fractions.gcd(n,algSide-ratSide)
		if(possibleFactor > 1):
			print "%s = %s*p" % (n,possibleFactor)
			break
		else:
			print "Trivial factor found, trying next dependency..."
		