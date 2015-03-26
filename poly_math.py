import math
import etc_math
from poly import *

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
	
	NFp = NumberFieldModP(poly,p)
	
	g = NFp(Poly([0,1])) ** p
	g = g - NFp(Poly([0,1]))
	g = g.getPoly()
	g = polynomialGCDModP(g,poly,p)
		
	roots.extend(rootsRecurse(g,p))
	return roots

def rootsRecurse(g,p):
	if(g.degree() == 2):
		return qe_mod(g.coeffs[2],g.coeffs[1],g.coeffs[0],p)
	elif(g.degree() == 1):
		return [(-g.coeffs[0]*modInv(g.coeffs[1],p)) % p]
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
			


	
