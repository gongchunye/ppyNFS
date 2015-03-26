from prime import *
from poly import *
import math


def findSmallFactors(n):
	n = abs(n)
	factors = [1,n]
	i = 2 
	while(i < n):
		if(n % i == 0):
			factors.append(i)
		i += 1	
		
	return factors
	
def shrinkNfsPoly(poly,m):
	# bound coeffs from -m/2 to m/2.
	for i in range(poly.degree()-1):
		coeff = poly.coeffs[i]
		if coeff > m/2:
			poly.coeffs[i] = -(m - coeff)
			poly.coeffs[i+1] = poly.coeffs[i+1]+1
			
def generateNFSPoly(n,d):
	m = int(math.pow(n,1.0/d))
	nfsPoly = Poly(expansionBaseM(n,m))
	shrinkNfsPoly(nfsPoly,m)
	while(reducible(nfsPoly)):
		m -= 1
		nfsPoly = Poly(expansionBaseM(n,m))
		shrinkNfsPoly(nfsPoly,m)
		
	return (m,nfsPoly)
	
def reducible(poly):
	possibleRoots = findSmallFactors(poly.coeffs[0])
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