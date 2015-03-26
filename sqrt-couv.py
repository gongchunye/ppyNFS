from poly import *
from etc_math import *
from files import *
from sqrt_funcs import *
from poly_math import *
		
def getPrimeExponents(value,base):
	primeExponents = [0]*len(base)
	for prime in base:
		baseIndex = base.index(prime)
		while(value % prime == 0):
			primeExponents[baseIndex] += 1
			value /= prime
			
	return primeExponents
	
def sumExponents(e1,e2):
	if(len(e1) != len(e2)):
		raise AssertionError
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
		smoothRat = NF(Poly(smooths[i])).norm()
		smoothPrimeExponents = getPrimeExponents(smoothRat,base)
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
		
def generatePrimes(nfsPoly,couvBound):		
	sumBound = 0
	primes = []
	while(sumBound <= couvBound):
		prime = generateLargePrime(64)
		
		if(not(irreducibleModP(nfsPoly,prime))):
			continue
			
		(s,r) = getsr(prime,3)
		if(r > 5):
			continue
			
		sumBound += math.log(prime,2)
		primes.append(prime)
		
	return primes
		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	
	rfBase = loadFileArray("rfbase.txt")
	smooths = loadFileArray("smooths.txt")
	deps = loadFileArray("deps.txt")
	
	for dependency in deps:
	
		dependencySmooths = []
		for column in dependency:
			smooth = smooths[dependency.index(column)]
			dependencySmooths.append(smooth)
			
		primeExponents = getRatPrimeExponents(dependencySmooths,rfBase)
		ratSide = getPrimeExponentsSqrtModN(primeExponents,rfBase,n)
		
		couvBound = calcRequiredPrimeLength(n,m,Poly(nfsPoly),dependencySmooths)
		primes = generatePrimes(Poly(nfsPoly),couvBound)
		primeExponents = getNormPrimeExponents(dependencySmooths,NF,rfBase)
		
		for prime in primes:
			normModP = getPrimeExponentsSqrtModN(primeExponents,rfBase,prime)
			
			NFp = NumberFieldModP(Poly(nfsPoly),prime)
			prod = NFp(Poly([1]))
			for smooth in dependencySmooths:
				prod = prod * NFp(Poly(smooth))
				
			sqrt = prod.sqrt()	
			
			if(sqrt.norm() == normModP or (-sqrt).norm() == normModP):
				print "Norms match"
			else:
				print "Norms no match"