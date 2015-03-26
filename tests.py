from poly import *
from nfspolygen import *
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
		poly = Poly([1,2,3])
		self.assertEqual(poly.degree(), 2)
		poly = Poly([0,2,3])
		self.assertEqual(poly.degree(), 2)
		poly = Poly([1,2,0])
		self.assertEqual(poly.degree(), 1)
		poly = Poly([1,2,3,4,5])
		self.assertEqual(poly.degree(), 4)
		
	def test_equal(self):
		poly1 = Poly([1,2,3,4])
		poly2 = Poly([1,2,3,4])
		poly3 = Poly([1,4,2,3])
		self.assertEqual(poly1==poly2,True)
		self.assertEqual(poly1,poly2)
		self.assertNotEqual(poly1,poly3)
		
	def test_numberfield_mult1(self):
		# Briggs' example
		nfspoly = Poly([8,29,15,1])
		NF = NumberField(nfspoly)
		
		nfspolyd = NF(nfspoly.derivative())
		correctProduct = NF(Poly([22939402657683071224L, 54100105785512562427L, 22455983949710645412L]))
		
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		tomult = [NF(Poly(x)) for x in tomult]
		
		prod = multAllElements(NF(Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_mult2(self):
		# Spaan's example
		nfspoly = Poly([161, 134, 2, 1])
		NF = NumberField(nfspoly)
		nfspolyd = NF(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NF(Poly(x)) for x in tomult]
		
		correctProduct = NF(Poly([21124198049840950371210079793023892077432,18523314201045731615331644622444823801483,884477920457388669411401815623954662863]))
		
		prod = multAllElements(NF(Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_square1(self):
		# Spaans's example
		nfspoly = Poly([161, 134, 2, 1])
		NF = NumberField(nfspoly)
		nfspolyd = NF(nfspoly.derivative())
		
		sqrt = NF(Poly([-41757429265579073242,-34105727598423382475,1681812579256330563]))
		correctSquare = NF(Poly([21124198049840950371210079793023892077432,18523314201045731615331644622444823801483,884477920457388669411401815623954662863]))
		
		square = sqrt * sqrt
		self.assertEqual(square,correctSquare)
		
	def test_numberfield_square2(self):
		# Brigg's example
		nfspoly = Poly([8,29,15,1])
		NF = NumberField(nfspoly)
		
		nfspolyd = NF(nfspoly.derivative())
		sqrt = NF(Poly([3889976768, 3686043120, 599923511]))
		correctSquare = NF(Poly([22939402657683071224L, 54100105785512562427L, 22455983949710645412L]))
		
		square = sqrt * sqrt
		self.assertEqual(square,correctSquare)
		
	def test_numberfield_modp1(self):
		# Briggs' example
		nfspoly = Poly([8,29,15,1])
		NFp = NumberFieldModP(nfspoly,9929)
		
		nfspolyd = NFp(nfspoly.derivative())
		correctProduct = NFp(Poly([6659,3891,2027]))
		
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		tomult = [NFp(Poly(x)) for x in tomult]
		prod = multAllElements(NFp(Poly([1])),tomult)
		
		
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp2(self):
		# Spaans' examples
		nfspoly = Poly([161, 134, 2, 1])
		NFp = NumberFieldModP(nfspoly,2305843009213693951)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(Poly(x)) for x in tomult]
		
		# Spaans provides the sqrt, so square that.
		sqrt = NFp(Poly([2053587909481111827,481917539782026790,1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp3(self):
		# test again for a different prime
		nfspoly = Poly([161, 134, 2, 1])
		NFp = NumberFieldModP(nfspoly,2305843009213693967)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(Poly(x)) for x in tomult]
		
		sqrt = NFp(Poly([2053587909481112131, 481917539782027030, 1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_modp3(self):
		# again, for a different prime
		nfspoly = Poly([161, 134, 2, 1])
		NFp = NumberFieldModP(nfspoly,2305843009213693973)
		nfspolyd = NFp(nfspoly.derivative())
		
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		tomult = [NFp(Poly(x)) for x in tomult]
		
		sqrt = NFp(Poly([2053587909481112245, 481917539782027120, 1681812579256330563]))
		correctProduct = sqrt * sqrt
		
		prod = multAllElements(NFp(Poly([1])),tomult)
		prod = prod * nfspolyd * nfspolyd
		self.assertEqual(prod,correctProduct)
		
	def test_numberfield_power_modp1(self):
		# Briggs
		nfspoly = Poly([8,29,15,1])
		NFp = NumberFieldModP(nfspoly,9929)
		base = NFp(Poly([6659,3891,2027]))
		s = 122356359011
		power = base ** s
		self.assertEqual(power,NFp(Poly([9928])))
		
		power = base ** ((s+1)/2)
		prod = power * NFp(Poly([7827]))
		self.assertEqual(prod,NFp(Poly([3077,1160,3402])))

		
		base = NFp(Poly([1,1]))
		power = base ** (2*s)
		self.assertEqual(power,NFp(Poly([2102])))
		
	def test_numberfield_power_modp2(self):
		#Briggs
		poly = Poly([8,29,15,1])
		p = 9923
		NFp = NumberFieldModP(poly,p)
		g = NFp(Poly([0,1])) ** p
		g = g - NFp(Poly([0,1]))
		g = g.getPoly() 
		
		correctG = Poly([7301,1477,7726])
		self.assertEqual(g,correctG)
		
		
		
	def test_numberfield_sqrt_modp1(self):
		# Briggs
		nfspoly = Poly([8,29,15,1])
		primes = [9851,9907,9929]
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		
		for prime in primes:
			NFp = NumberFieldModP(nfspoly,prime)
			tomultNFp = [NFp(Poly(x)) for x in tomult]
			nfspolyd = NFp(nfspoly.derivative())
			
			prod = multAllElements(NFp(Poly([1])),tomultNFp)
			prod = prod * nfspolyd * nfspolyd
			sqrt = prod.sqrt()
			
			self.assertEqual(sqrt*sqrt,prod)
			
	def test_numberfield_sqrt_modp2(self):
		# Briggs
		nfspoly = Poly([161, 134, 2, 1])
		primes = [2305843009213693951,2305843009213693967,2305843009213693973,2305843009213694381]
		tomult = [[-92,-1],[-57,-1],[-23,-1],[-8,-1],[-7,-1],[2,-1],[10,-1],[17,-1],[29,-1],[35,-1],[84,-1],[115,-1],[139,-1],[-5,-2],[19,-2],[69,-2],[93,-2],[119,-2],[-542,-3],[-28,-3],[-23,-3],[-8,-3]]
		
		for prime in primes:
			NFp = NumberFieldModP(nfspoly,prime)
			tomultNFp = [NFp(Poly(x)) for x in tomult]
			nfspolyd = NFp(nfspoly.derivative())
			
			prod = multAllElements(NFp(Poly([1])),tomultNFp)
			prod = prod * nfspolyd * nfspolyd
			sqrt = prod.sqrt()
			self.assertEqual(sqrt*sqrt,prod)
		
	def test_numberfield_positivesquareroot_modp1(self):
		# Briggs' CRT example
		n = 45113
		m = 31
		nfspoly = Poly([8,29,15,1])
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		primeSizeEst = calcRequiredPrimeLength(n, m, nfspoly, tomult)
		#primes = randomPrimes(primeSizeEst,32)
		prime = 9929
		
		NFp = NumberFieldModP(nfspoly,prime)
		nfspolyd = NFp(nfspoly.derivative())
		tomultNFp = [NFp(Poly(x)) for x in tomult]
		prod = multAllElements(NFp(Poly([1])),tomultNFp)
		prod = prod * nfspolyd * nfspolyd
			
		sqrt = prod.sqrt()
		posSqrt = NFp(Poly([3077, 1160, 3402]))
		negSqrt = NFp(Poly([6852, 8769, 6527]))
		
		self.assertTrue(sqrt == posSqrt or sqrt == negSqrt)

	def test_numberfield_sqrt_equality(self):
		nfspoly = Poly([8,29,15,1])
		prime = 9929
		NFp = NumberFieldModP(nfspoly,prime)
		posSqrt = NFp(Poly([3077, 1160, 3402]))
		negSqrt = NFp(Poly([6852, 8769, 6527]))
		
		self.assertEqual(posSqrt,-negSqrt)
		self.assertEqual(-posSqrt,negSqrt)
	
	def test_baseExpansion(self):
		#Briggs
		correctNfspoly = Poly([8,29,15,1])
		testNfspoly = Poly(expansionBaseM(45113,31))
		self.assertEqual(correctNfspoly,testNfspoly)

	def test_reduciblePolynomial1(self):
		#Briggs
		poly = Poly([8,29,15,1])
		self.assertFalse(reducible(poly))
		
	def test_reduciblePolynomial(self):
		poly = Poly([-6,11,-6,1])
		self.assertTrue(reducible(poly))
		
	def test_nfsPolyGeneration(self):
		n = 45113
		d = 3
		(m,nfspoly) = generateNFSPoly(n,d)
		self.assertFalse(reducible(nfspoly))
		self.assertEqual(nfspoly.coeffs[-1],1)
		self.assertEqual(nfspoly.evaluate(m),n)
		
	def test_nfsPolyGeneration64bit(self):
		n = 8202545090182721807
		d = 3
		(m,nfspoly) = generateNFSPoly(n,d)
		self.assertFalse(reducible(nfspoly))
		self.assertEqual(nfspoly.coeffs[-1],1)
		self.assertEqual(nfspoly.evaluate(m),n)

	def test_polyRootModPSlow(self):
		nfspoly = Poly([8,29,15,1])
		correctRoots = [2,44,6]
		testRoots = getRootsModPSlow(nfspoly,67)
		
		self.assertEqual(len(correctRoots), len(testRoots))
		
		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polySubtract1(self):
		poly1 = Poly([8,29,15,1])
		poly2 = Poly([1,1])
		testDiff = poly1 - poly2
		correctDiff = Poly([7,28,15,1])
		self.assertEqual(correctDiff,testDiff)
		
	def test_polySubtract2(self):
		poly1 = Poly([])
		poly2 = Poly([1,1])
		testDiff = poly1 - poly2
		
		correctDiff = Poly([-1,-1])
		self.assertEqual(correctDiff,testDiff)
		
	
	def test_polyRootModPFast(self):
		# Briggs
		poly = Poly([8,29,15,1])
		testRoots = getRootsModPFast(poly,67)
		correctRoots = [2,44,6]
		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyRootModPFastRandomDeg3(self):
		
		p = 503
		poly = Poly([2034,234,24,123])
		correctRoots = getRootsModPSlow(poly,p)
		
		testRoots = getRootsModPFast(poly,p)

		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyRootModPFastRandomZeroRoot(self):
		
		p = 503
		poly = Poly([0,1])*Poly([-11,1])*Poly([-51,1])*Poly([-231,1])
		correctRoots = getRootsModPSlow(poly,p)
		testRoots = getRootsModPFast(poly,p)

		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
		
	def test_polyRootModPFastRandomDeg5(self):
		# no roots
		p = 503
		poly = Poly([2034,234,24,123,101,1])
		correctRoots = getRootsModPSlow(poly,p)
		testRoots = getRootsModPFast(poly,p)
		
		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
		
		
	def test_polyRootModPFastDeg5_2(self):
		p = 157
		poly = Poly([11,1])*Poly([-23,1])*Poly([-1,1])*Poly([-1,1])*Poly([-1,1])
		correctRoots = getRootsModPSlow(poly,p)
		testRoots = getRootsModPFast(poly,p)

		correctRoots.sort()
		testRoots.sort()
		
		for i in range(len(correctRoots)):
			self.assertEqual(correctRoots[i], testRoots[i])
			
	def test_polyGCDModP(self):
		poly1 = Poly([7301,1477,7726])
		poly2 = Poly([8,29,15,1])
		testGCD = polynomialGCDModP(poly1,poly2,9923)
		correctGCD = Poly([9858, 7744])
		self.assertEqual(testGCD,correctGCD)
	
	
	def test_reduceToNFp(self):
		p = 9923
		poly1 = Poly([7301,1477,7726])
		poly2 = Poly([8,29,15,1])
		NFp = NumberFieldModP(poly1,p)
		testReduction = NFp(poly2)
		correctReduction = NFp(Poly([9858, 7744]))
		self.assertEqual(testReduction,correctReduction)
		
	def test_normPolyB1(self):
		NF = NumberField(Poly([8,29,15,1]))
		b = -5
		testPolyNormB = NF.getPolyNormB(b)
		correctPolyNormB = Poly([1000, 725, 75, 1])
		self.assertEqual(testPolyNormB, correctPolyNormB)
	
	def test_normPolyB2(self):
		NF = NumberField(Poly([8,29,15,1]))
		a = -8
		b = 3
		NF = NumberField(Poly([8,29,15,1]))
		polyNormB = NF.getPolyNormB(b)
		self.assertEqual(polyNormB.evaluate(a),-5696)
		
	def test_norm(self):
		NF = NumberField(Poly([8,29,15,1]))
		smoothElement = NF(Poly([-8,3]))
		self.assertEqual(smoothElement.norm(),-5696)
	
if __name__ == '__main__':
	unittest.main()
	
