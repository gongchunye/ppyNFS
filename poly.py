import fractions
import etcmath
import random

def getRootsModPSlow(poly,p):
	roots = []
	for i in range(p):
		if(poly.evaluate(i)%p == 0):
			roots.append(i)
	return roots

def getRootsModPFast(poly,p):
	roots = []
	
	if(poly.degree() > p):
		return getRootsModPSlow(poly,p)
	
	g = getG(poly,p)
		
	roots.extend(rootsRecurse(g,p))
	return roots
	
def getG(poly,p):
	NFp = NumberFieldModP(poly,p)
		
	g = NFp(Poly([0,1])) ** p
	g = g - NFp(Poly([0,1]))
	g = g.getPoly()
	g = polynomialGCDModP(g,poly,p)
	
	return g
	
def irreducibleModP(poly,p):
	if(getG(poly,p).degree() == 0):
		return True
	else:
		return False

def rootsRecurse(g,p):
	if(g.degree() == 2):
		return etcmath.qe_mod(g.coeffs[2],g.coeffs[1],g.coeffs[0],p)
	elif(g.degree() == 1):
		return [(-g.coeffs[0]*etcmath.modInv(g.coeffs[1],p)) % p]
	elif(g.degree() == 0):
		return []
	
	roots = []
	h = Poly([1])
	NFp = NumberFieldModP(g,p)
	
	while(h.degree() == 0 or g == h):
		a = random.randint(0,p)
		base = NFp(Poly([a,1]))
		base = base ** ((p-1)/2)
		toGCD = base - NFp(Poly([1])) 
		h = polynomialGCDModP(toGCD.getPoly(),g,p)
		
	roots.extend(rootsRecurse(h,p))
	dividedPoly = NFp(g/h).getPoly()
	roots.extend(rootsRecurse(dividedPoly,p))
	return roots

def polynomialGCDModP(u,v,p):
	while(v.degree() >= 0):	
		NFp = NumberFieldModP(v,p)
		u = NFp(u) # reduce u(x) to u(x) mod v(x) mod p
		u = u.getPoly() # then extract poly
		(u,v) = (v,u)
	return u


def findnonresidue(NFp,e):
	poly = NFp(Poly([1,1]))
	while (True):
		test = poly**((NFp.prime**e-1)/2)
		if(test.degree() > 0):
			raise AssertionError
		if(test ==  NFp(Poly([NFp.prime-1]))):
			return poly
		
		nextPoly = poly.getPoly()
		nextPoly.coeffs[0] += 1
		poly = NFp(nextPoly)

class Poly():
	def __init__(self, coeffs):
		self.coeffs = coeffs
		
	def degree(self):
		while self.coeffs and self.coeffs[-1] == 0:
			self.coeffs.pop()   # normalize
		return len(self.coeffs)-1	
		
	def __str__(self):
		self.degree()
		return str(self.coeffs)
		
	def __eq__(self,other):
		if(self.degree() != other.degree()):
			return False
			
		for i in range(len(self.coeffs)):
			if(self.coeffs[i] != other.coeffs[i]):
				return False
			
		return True
		
	def __mul__ (self,other):
		res = [0]*(len(self.coeffs)+len(other.coeffs)-1)
		for o1,i1 in enumerate(self.coeffs):
			for o2,i2 in enumerate(other.coeffs):
				res[o1+o2] += i1*i2
				
		return Poly(res);
		
	def __str__(self):		
		return str(self.coeffs)
		
	def derivative(self):
		coeffs = [0]*len(self.coeffs)
		for i in range(len(coeffs)-1):
			coeffs[i] = (i+1)*self.coeffs[i+1]
		
		p = Poly(coeffs)
		p.degree()
		return Poly(coeffs)
		
	def modcoeffs(self,p):
		coeffs = [0]*len(self.coeffs)
		for i in range(len(self.coeffs)):
			if(type(self.coeffs[i]) == fractions.Fraction):
				coeffs[i] = etcmath.modInv(self.coeffs[i].denominator,p)*self.coeffs[i].numerator % p
			else:
				coeffs[i] = self.coeffs[i] % p
		return Poly(coeffs)
		
	def getCoeffs(self):
		return self.coeffs
		
		
	def evaluate(self,m):
		self.degree()
		
		x = 1
		sum = 0
		for coeff in self.coeffs:
			sum += coeff*x
			x *= m
			
		return sum
		
	def __sub__(self,other):
		selfCopy = self.coeffs[:]
		otherCopy = other.coeffs[:]
		if(len(selfCopy) > len(otherCopy)):
			otherCopy.extend([0]*(len(selfCopy)-len(otherCopy)))
		elif(other.degree() > self.degree()):
			selfCopy.extend([0]*(len(otherCopy)-len(selfCopy)))
	
		diffcoeffs = []
		
		for i in range(len(selfCopy)):
			diffcoeffs.append(selfCopy[i] - otherCopy[i])

		return Poly(diffcoeffs)
		
	def __mod__(self,other):
		(q,r) = self.divmod(other)
		return r
	def __div__(self,other):
		(q,r) = self.divmod(other)
		return q
	def __normalize(self,poly):
		while poly and poly[-1] == 0:
			poly.pop()
		if poly == []:
			poly.append(0)
			
	# this is from http://stackoverflow.com/questions/26173058/division-of-polynomials-in-python
	def divmod(self,other):
		num = self.coeffs[:]
		self.__normalize(num)
		den = other.coeffs[:]
		self.__normalize(den)
		quot = []
		
		if len(num) >= len(den):
			#Shift den towards right so it's the same degree as num
			shiftlen = len(num) - len(den)
			den = [0] * shiftlen + den
		else:
			return  (Poly(quot),Poly(num))

		
		divisor = den[-1]
		for i in xrange(shiftlen + 1):
			#Get the next coefficient of the quotient.
			
	
			if(num[-1] % divisor != 0):
				mult = fractions.Fraction(num[-1],divisor)
			else:
				mult = num[-1]/divisor
				
			quot = [mult] + quot

			#Subtract mult * den from num, but don't bother if mult == 0
			#Note that when i==0, mult!=0; so quot is automatically normalized.
			if mult != 0:
				d = [mult * u for u in den]
				num = [(u - v) for u, v in zip(num, d)]

			num.pop()
			den.pop(0)

		self.__normalize(num)
		self.__normalize(quot)
		
		return (Poly(quot),Poly(num))
		
'''
the following classes should only export the 
Poly class (i.e. non-field element) through getPoly().
'''
		
def NumberField(numberFieldPoly):
	class _NumberField():
	
		def __init__(self, __poly):
			self.__poly = __poly % numberFieldPoly #self.__modulo(__poly)
			
		def degree(self):
			return self.__poly.degree()
			
		def __mul__(self, other): 
			return _NumberField((self.__poly * other.__poly) % numberFieldPoly)
			
			
		def __str__(self):		
			return "%s (mod %s)" % (str(self.__poly.coeffs),str(self.numberFieldPoly.coeffs))
			
		def __eq__(self, other):		
			return self.__poly == other.__poly
			
		def getCoeffs(self):
			return self.__poly.getCoeffs()
			
		def __sub__(self,other):
			return _NumberField(self.__poly - other.__poly)
			
		def modcoeffs(self,p):
			return _NumberField(self.__poly.modcoeffs(p))	
			
		@staticmethod
		def getPolyNormB(b):
			d = numberFieldPoly.degree()
			coeffs = []
			for i in range(d+1):
				if(i%2 == 1):
					coeffs.append(numberFieldPoly.coeffs[i]*(b**(d-i)))
				else:
					coeffs.append(-numberFieldPoly.coeffs[i]*(b**(d-i)))
			poly = Poly(coeffs)
			return poly
			
		def norm(self):
			a=self.__poly.coeffs[0]
			b=self.__poly.coeffs[1]
			d = numberFieldPoly.degree()
			sum = 0
			for i in range(d+1):
				if(i%2 == 0):
					sum += numberFieldPoly.coeffs[d-i]*(a**(d-i))*(b**i)
				else:
					sum -= numberFieldPoly.coeffs[d-i]*(a**(d-i))*(b**i)
			return sum
			
		def getPoly(self):
			return Poly(self.__poly.coeffs[:])
			
	_NumberField.numberFieldPoly = numberFieldPoly
	return _NumberField
	
def NumberFieldModP(numberFieldPoly,prime):
	class _NumberFieldModP():
	
		def __init__(self, poly):
			poly = self.NF(poly)
			self.__polyInNF = poly.modcoeffs(prime)
			
		def degree(self):
			return self.__polyInNF.degree()
			
		def __mul__(self, other): 
			product = self.__polyInNF * other.__polyInNF
			return _NumberFieldModP(product.getPoly().modcoeffs(prime))
			
		def __str__(self):		
			return "%s (mod %s)" % (self.__polyInNF,self.prime)
			
		def __eq__(self, other):		
			return self.__polyInNF == other.__polyInNF
			
		def __neg__(self):
			negativeCoeffs = Poly([-x for x in self.__polyInNF.getPoly().coeffs])
			return _NumberFieldModP(negativeCoeffs)
			
		def __pow__(self,other):
			base = self
			power = other
			prod=_NumberFieldModP(Poly([1]))
			while(power > 0):
				if(power&1==1): prod = prod * base 
				base = base * base
				power >>= 1
			return prod
			
		def __sub__(self,other):
			difference = self.__polyInNF - other.__polyInNF
			return _NumberFieldModP(difference.getPoly()) 
			
		def sqrt(self):
			if(self.nonresidue == None):
				self.nonresidue = findnonresidue(_NumberFieldModP, 3)
				(self.s,self.r) = etcmath.getsr(prime, 3)
			
			'''
			sylowSquare = self ** self.s
			#print sylowSquare
			for i in range(2**self.r):
				sylowSquareTrial = self.nonresidue**(2*i*self.s)
				if(sylowSquareTrial == sylowSquare):
					sylowSqrt = self.nonresidue**(i*self.s)
					inverse = modinv(sylowSqrt.getCoeffs()[0],self.prime)
					inverse = _NumberFieldModP(Poly([inverse]))
					return inverse*(self ** ((self.s+1)/2))
			'''

			lamb = self ** self.s
			zeta = self.nonresidue**self.s
			w = self ** ((self.s+1)/2)
		
			while(lamb.getPoly().coeffs[0] != 1):
				m = 1
				while((lamb ** (2**m)).getPoly().coeffs[0] != 1):
					m += 1
				w = w*(zeta**(2**(self.r-m-1)))
				lamb = lamb*(zeta**(2**(self.r-m)))
				
			return w
		
		def getCoeffs(self):
			return self.__polyInNF.getCoeffs()
		
		def getPoly(self):
			return self.__polyInNF.getPoly()
		
		def norm(self):
			normPoly = self ** ((self.prime**3 - 1)/(self.prime - 1))
			return normPoly.getCoeffs()[0]
			
			
	_NumberFieldModP.numberFieldPoly = numberFieldPoly
	_NumberFieldModP.NF = NumberField(numberFieldPoly)
	_NumberFieldModP.prime = prime
	_NumberFieldModP.nonresidue = None
	(_NumberFieldModP.s,_NumberFieldModP.r) = (None,None)
	return _NumberFieldModP