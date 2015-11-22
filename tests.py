import poly
import etcmath
import nfspolygen
import unittest
import math

def multAllElements(one,array):
	prod = one
	for entry in array:
		prod = prod * entry
		
	return prod
	

def randomPrimes(totalsize,size):
	primes = [] 
	sum = 0
	while(sum <= totalsize):
		p = generateLargePrime(size)
		primes.append(p)
		sum += math.log(p,2)
		
	return primes

class TestSequenceFunctions(unittest.TestCase):
	
	def test_degree(self):
		polynomial = poly.Poly([1,2,3])
		self.assertEqual(polynomial.degree(), 2)
		polynomial = poly.Poly([0,2,3])
		self.assertEqual(polynomial.degree(), 2)
		polynomial = poly.Poly([1,2,0])
		self.assertEqual(polynomial.degree(), 1)
		polynomial = poly.Poly([1,2,3,4,5])
		self.assertEqual(polynomial.degree(), 4)
		
	def test_equal(self):
		poly1 = poly.Poly([1,2,3,4])
		poly2 = poly.Poly([1,2,3,4])
		poly3 = poly.Poly([1,4,2,3])
		self.assertEqual(poly1==poly2,True)
		self.assertEqual(poly1,poly2)
		self.assertNotEqual(poly1,poly3)
		
	def test_numberfield_mult1(self):
		# Briggs' example
		nfspoly = poly.Poly([8,29,15,1])
		NF = poly.NumberField(nfspoly)
		
		nfspolyd = NF(nfspoly.derivative())
		correctProduct = NF(poly.Poly([22939402657683071224L, 54100105785512562427L, 22455983949710645412L]))
		
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		tomult = [NF(poly.Poly(x)) for x in tomult]
		
		prod = multAllElements(NF(poly.Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_mult2(self):
		# Spaan's example
		nfspoly = poly.Poly([161, 134, 2, 1])
		NF = poly.NumberField(nfspoly)
		nfspolyd = NF(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NF(poly.Poly(x)) for x in tomult]
		
		correctProduct = NF(poly.Poly([21124198049840950371210079793023892077432,18523314201045731615331644622444823801483,884477920457388669411401815623954662863]))
		
		prod = multAllElements(NF(poly.Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_square1(self):
		# Spaans's example
		nfspoly = poly.Poly([161, 134, 2, 1])
		NF = poly.NumberField(nfspoly)
		nfspolyd = NF(nfspoly.derivative())
		
		sqrt = NF(poly.Poly([-41757429265579073242,-34105727598423382475,1681812579256330563]))
		correctSquare = NF(poly.Poly([21124198049840950371210079793023892077432,18523314201045731615331644622444823801483,884477920457388669411401815623954662863]))
		
		square = sqrt * sqrt
		self.assertEqual(square,correctSquare)
		
	def test_numberfield_square2(self):
		# Brigg's example
		nfspoly = poly.Poly([8,29,15,1])
		NF = poly.NumberField(nfspoly)
		
		nfspolyd = NF(nfspoly.derivative())
		sqrt = NF(poly.Poly([3889976768, 3686043120, 599923511]))
		correctSquare = NF(poly.Poly([22939402657683071224L, 54100105785512562427L, 22455983949710645412L]))
		
		square = sqrt * sqrt
		self.assertEqual(square,correctSquare)
		
	def test_numberfield_modp1(self):
		# Briggs' example
		nfspoly = poly.Poly([8,29,15,1])
		NFp = poly.NumberFieldModP(nfspoly,9929)
		
		nfspolyd = NFp(nfspoly.derivative())
		correctProduct = NFp(poly.Poly([6659,3891,2027]))
		
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		tomult = [NFp(poly.Poly(x)) for x in tomult]
		prod = multAllElements(NFp(poly.Poly([1])),tomult)
		
		
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp2(self):
		# Spaans' examples
		nfspoly = poly.Poly([161, 134, 2, 1])
		NFp = poly.NumberFieldModP(nfspoly,2305843009213693951)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(poly.Poly(x)) for x in tomult]
		
		# Spaans provides the sqrt, so square that.
		sqrt = NFp(poly.Poly([2053587909481111827,481917539782026790,1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(poly.Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp3(self):
		# test again for a different prime
		nfspoly = poly.Poly([161, 134, 2, 1])
		NFp = poly.NumberFieldModP(nfspoly,2305843009213693967)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(poly.Poly(x)) for x in tomult]
		
		sqrt = NFp(poly.Poly([2053587909481112131, 481917539782027030, 1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(poly.Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp3(self):
		# again, for a different prime
		nfspoly = poly.Poly([161, 134, 2, 1])
		NFp = poly.NumberFieldModP(nfspoly,2305843009213693973)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(poly.Poly(x)) for x in tomult]
		
		sqrt = NFp(poly.Poly([2053587909481112245, 481917539782027120, 1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(poly.Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_power_modp1(self):
		# Briggs
		nfspoly = poly.Poly([8,29,15,1])
		NFp = poly.NumberFieldModP(nfspoly,9929)
		base = NFp(poly.Poly([6659,3891,2027]))
		s = 122356359011
		power = base ** s
		self.assertEqual(power,NFp(poly.Poly([9928])))
		
		power = base ** ((s+1)/2)
		prod = power * NFp(poly.Poly([7827]))
		self.assertEqual(prod,NFp(poly.Poly([3077,1160,3402])))

		
		base = NFp(poly.Poly([1,1]))
		power = base ** (2*s)
		self.assertEqual(power,NFp(poly.Poly([2102])))
		
	def test_numberfield_power_modp2(self):
		#Briggs
		polynomial = poly.Poly([8,29,15,1])
		p = 9923
		NFp = poly.NumberFieldModP(polynomial,p)
		g = NFp(poly.Poly([0,1])) ** p
		g = g - NFp(poly.Poly([0,1]))
		g = g.getPoly() 
		
		correctG = poly.Poly([7301,1477,7726])
		self.assertEqual(g,correctG)
		
		
		
	def test_numberfield_sqrt_modp1(self):
		# Briggs
		nfspoly = poly.Poly([8,29,15,1])
		primes = [9851,9907,9929]
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		
		for prime in primes:
			NFp = poly.NumberFieldModP(nfspoly,prime)
			tomultNFp = [NFp(poly.Poly(x)) for x in tomult]
			nfspolyd = NFp(nfspoly.derivative())
			
			prod = multAllElements(NFp(poly.Poly([1])),tomultNFp)
			prod = prod * nfspolyd * nfspolyd
			sqrt = prod.sqrt()
			
			self.assertEqual(sqrt*sqrt,prod)
			
	def test_numberfield_sqrt_modp2(self):
		# Briggs
		nfspoly = poly.Poly([161, 134, 2, 1])
		primes = [2305843009213693951,2305843009213693967,2305843009213693973,2305843009213694381]
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		
		for prime in primes:
			NFp = poly.NumberFieldModP(nfspoly,prime)
			tomultNFp = [NFp(poly.Poly(x)) for x in tomult]
			nfspolyd = NFp(nfspoly.derivative())
			
			prod = multAllElements(NFp(poly.Poly([1])),tomultNFp)
			prod = prod * nfspolyd * nfspolyd
			sqrt = prod.sqrt()
			self.assertEqual(sqrt*sqrt,prod)
		
	def test_numberfield_positivesquareroot_modp1(self):
		# Briggs' CRT example
		n = 45113
		m = 31
		nfspoly = poly.Poly([8,29,15,1])
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		primeSizeEst = etcmath.calcRequiredPrimeLength(n, m, nfspoly, tomult)
		#primes = randomPrimes(primeSizeEst,32)
		prime = 9929
		
		NFp = poly.NumberFieldModP(nfspoly,prime)
		nfspolyd = NFp(nfspoly.derivative())
		tomultNFp = [NFp(poly.Poly(x)) for x in tomult]
		prod = multAllElements(NFp(poly.Poly([1])),tomultNFp)
		prod = prod * nfspolyd * nfspolyd
			
		sqrt = prod.sqrt()
		posSqrt = NFp(poly.Poly([3077, 1160, 3402]))
		negSqrt = NFp(poly.Poly([6852, 8769, 6527]))
		
		self.assertTrue(sqrt == posSqrt or sqrt == negSqrt)

	def test_numberfield_sqrt_equality(self):
		nfspoly = poly.Poly([8,29,15,1])
		prime = 9929
		NFp = poly.NumberFieldModP(nfspoly,prime)
		posSqrt = NFp(poly.Poly([3077, 1160, 3402]))
		negSqrt = NFp(poly.Poly([6852, 8769, 6527]))
		
		self.assertEqual(posSqrt,-negSqrt)
		self.assertEqual(-posSqrt,negSqrt)
	
	def test_baseExpansion(self):
		#Briggs
		correctNfspoly = poly.Poly([8,29,15,1])
		testNfspoly = poly.Poly(nfspolygen.expansionBaseM(45113,31))
		self.assertEqual(correctNfspoly,testNfspoly)

	def test_reduciblePolynomial1(self):
		#Briggs
		polynomial = poly.Poly([8,29,15,1])
		self.assertFalse(nfspolygen.reducible(polynomial))
		
	def test_reduciblePolynomial(self):
		polynomial = poly.Poly([-6,11,-6,1])
		self.assertTrue(nfspolygen.reducible(polynomial))
		
	def test_nfsPolyGeneration(self):
		n = 45113
		d = 3
		(m,nfspoly) = nfspolygen.generateNFSPoly(n,d)
		self.assertFalse(nfspolygen.reducible(nfspoly))
		self.assertEqual(nfspoly.coeffs[-1],1)
		self.assertEqual(nfspoly.evaluate(m),n)
		
	def test_nfsPolyGeneration64bit(self):
		n = 8202545090182721807
		d = 3
		(m,nfspoly) = nfspolygen.generateNFSPoly(n,d)
		self.assertFalse(nfspolygen.reducible(nfspoly))
		self.assertEqual(nfspoly.coeffs[-1],1)
		self.assertEqual(nfspoly.evaluate(m),n)

	def test_polyRootModPSlow(self):
		nfspoly = poly.Poly([8,29,15,1])
		correctRoots = [2,44,6]
		testRoots = poly.getRootsModPSlow(nfspoly,67)
		
		self.assertEqual(len(correctRoots), len(testRoots))
		
		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polySubtract1(self):
		poly1 = poly.Poly([8,29,15,1])
		poly2 = poly.Poly([1,1])
		testDiff = poly1 - poly2
		correctDiff = poly.Poly([7,28,15,1])
		self.assertEqual(correctDiff,testDiff)
		
	def test_polySubtract2(self):
		poly1 = poly.Poly([])
		poly2 = poly.Poly([1,1])
		testDiff = poly1 - poly2
		
		correctDiff = poly.Poly([-1,-1])
		self.assertEqual(correctDiff,testDiff)
		
	
	def test_polyRootModPFast(self):
		# Briggs
		polynomial = poly.Poly([8,29,15,1])
		testRoots = poly.getRootsModPFast(polynomial,67)
		correctRoots = [2,44,6]
		correctRoots.sort()
		testRoots.sort()
		
		self.assertEqual(len(correctRoots),len(testRoots))
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyRootModPFastRandomDeg3(self):
		
		p = 503
		polynomial = poly.Poly([2034,234,24,123])
		correctRoots = poly.getRootsModPSlow(polynomial,p)
		
		testRoots = poly.getRootsModPFast(polynomial,p)

		correctRoots.sort()
		testRoots.sort()
		
		self.assertEqual(len(correctRoots),len(testRoots))
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyRootModPFastRandomZeroRoot(self):
		
		p = 503
		polynomial = poly.Poly([0,1])*poly.Poly([-11,1])*poly.Poly([-51,1])*poly.Poly([-231,1])
		correctRoots = poly.getRootsModPSlow(polynomial,p)
		testRoots = poly.getRootsModPFast(polynomial,p)

		correctRoots.sort()
		testRoots.sort()
		
		self.assertEqual(len(correctRoots),len(testRoots))
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
		
	def test_polyRootModPFastRandomDeg5(self):
		# no roots
		p = 503
		polynomial = poly.Poly([2034,234,24,123,101,1])
		correctRoots = poly.getRootsModPSlow(polynomial,p)
		testRoots = poly.getRootsModPFast(polynomial,p)
		
		correctRoots.sort()
		testRoots.sort()
		
		self.assertEqual(len(correctRoots),len(testRoots))
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
		
		
	def test_polyRootModPFastDeg5_2(self):
		p = 157
		polynomial = poly.Poly([11,1])*poly.Poly([-23,1])*poly.Poly([-1,1])*poly.Poly([-1,1])*poly.Poly([-1,1])
		correctRoots = poly.getRootsModPSlow(polynomial,p)
		testRoots = poly.getRootsModPFast(polynomial,p)

		correctRoots.sort()
		testRoots.sort()
		
		self.assertEqual(len(correctRoots),len(testRoots))
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyGCDModP(self):
		poly1 = poly.Poly([7301,1477,7726])
		poly2 = poly.Poly([8,29,15,1])
		testGCD = poly.polynomialGCDModP(poly1,poly2,9923)
		correctGCD = poly.Poly([9858, 7744])
		self.assertEqual(testGCD,correctGCD)
	
	
	def test_reduceToNFp(self):
		p = 9923
		poly1 = poly.Poly([7301,1477,7726])
		poly2 = poly.Poly([8,29,15,1])
		NFp = poly.NumberFieldModP(poly1,p)
		testReduction = NFp(poly2)
		correctReduction = NFp(poly.Poly([9858, 7744]))
		self.assertEqual(testReduction,correctReduction)
		
	def test_normPolyB1(self):
		NF = poly.NumberField(poly.Poly([8,29,15,1]))
		b = -5
		testPolyNormB = NF.getPolyNormB(b)
		correctPolyNormB = poly.Poly([1000, 725, 75, 1])
		self.assertEqual(testPolyNormB, correctPolyNormB)
	
	def test_normPolyB2(self):
		NF = poly.NumberField(poly.Poly([8,29,15,1]))
		a = -8
		b = 3
		NF = poly.NumberField(poly.Poly([8,29,15,1]))
		polyNormB = NF.getPolyNormB(b)
		self.assertEqual(polyNormB.evaluate(a),-5696)
		
	def test_norm(self):
		NF = poly.NumberField(poly.Poly([8,29,15,1]))
		smoothElement = NF(poly.Poly([-8,3]))
		self.assertEqual(smoothElement.norm(),-5696)
		
	def test_divisors(self):
		n = 120
		correctDivisors = [1,2,3,4,5,6,8,10,12,15,20,24,30,40,60,120]
		testDivisors = nfspolygen.findDivisors(n)
		
		correctDivisors.sort()
		testDivisors.sort()
		
		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])
			
	def test_divisors2(self):
		n = 72
		correctDivisors = [1, 2, 3, 4, 6, 8, 9, 12, 18, 24, 36, 72]
		testDivisors = nfspolygen.findDivisors(n)
		
		correctDivisors.sort()
		testDivisors.sort()
		
		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])
			
	def test_divisors_corner1(self):
		correctDivisors = [1]
		testDivisors = nfspolygen.findDivisors(1)
		
		correctDivisors.sort()
		testDivisors.sort()
		
		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])
		
	def test_divisors_corner2(self):
		correctDivisors = [1,2]
		testDivisors = nfspolygen.findDivisors(2)
		
		correctDivisors.sort()
		testDivisors.sort()

		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])

		
	def test_divisors_corner2(self):
		correctDivisors = [1,2]
		testDivisors = nfspolygen.findDivisors(2)
		
		correctDivisors.sort()
		testDivisors.sort()

		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])			
			
	def test_divisors_corner3(self):
		correctDivisors = [1,73]
		testDivisors = nfspolygen.findDivisors(73)
		
		correctDivisors.sort()
		testDivisors.sort()

		self.assertEqual(len(correctDivisors),len(testDivisors))
		for i in range(len(correctDivisors)):
			self.assertEqual(correctDivisors[i], testDivisors[i])			
			
if __name__ == '__main__':
	unittest.main()
	
