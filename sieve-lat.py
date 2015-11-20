import poly
import primemath
import liblll
import etcmath
import files
import math
import fractions

def frankefy(basis):
	a0 = basis[0][0]
	b0 = basis[0][1]
	a1 = basis[1][0]
	b1 = basis[1][1]
	return (a0,b0,a1,b1)

def sieve_quadrant(q1,q2,I,J,basis,sieve,plog):
	(a0,b0,a1,b1) = frankefy(basis)
	l = 0
	while(True):
		l += q1
		k = 0
		i = a0*k+a1*l
		j = b0*k+b1*l
		
		if(i > I/2 or i < -I/2 or j > J/2 or j < -J/2):
			return
		while(True):
			k += q2
			i = a0*k+a1*l
			j = b0*k+b1*l
		
			if(i > I/2 or i < -I/2 or j > J/2 or j < -J/2):
				break
			
			sieve[i+I/2-1][j+J/2-1] += plog
	
def reduce_basis(p1,p2):
	mat_basis = liblll.create_matrix([[p1[0],p2[0]],[p1[1],p2[1]]]) 
	mat_reduced = liblll.lll_reduction(mat_basis)	
	a0 = mat_reduced[0][0].numerator # lll library creates fractions with denominator 1, so grab numerator
	a1 = mat_reduced[0][1].numerator
	b0 = mat_reduced[1][0].numerator
	b1 = mat_reduced[1][1].numerator
	return [[a0,b0],[a1,b1]]
	
	
def find_sublattice_basis(basis_q,p,r):
	(a0,b0,a1,b1) = frankefy(basis_q)
	j = 0
	while(True):
		j += 1
		#print "%s mod %s" % (a0+r*b0,p)
		if((a0+r*b0) % p == 0 or j > p):
			return False
		i = (-j*etcmath.modInv(a0+r*b0,p)*(a1+r*b1)) % p
		if(i == 0 or j == 0):
			continue			
		basis1 = [i,j]
		break
		
	while(True):
		j += 1
		i = (-j*etcmath.modInv(a0+r*b0,p)*(a1+r*b1)) % p

		if(i == 0 or j == 0):
			continue		
		if(float(basis1[0])/i - float(basis1[1])/j < 0.01):
			continue		
		basis2 = [i,j]
		break

	return [basis1,basis2]
	
def find_lattice_basis(q,s):
	b=0
	while(True):
		b += 1
		a = -b*s % q
		if(a == 0 or b == 0):
			continue			
		basis1 = [a,b]
		break
		
	while(True):
		b += 1
		a = -b*s % q
		if(a == 0 or b == 0):
			continue		
		if(float(basis1[0])/a - float(basis1[1])/b < 0.01):
			continue		
		basis2 = [a,b]
		break

	return [basis1,basis2]
	
def sieve_special_q(q,s,I,J,rfBase,afBase,(n,nfsPoly,m,B,M,K)):
	nfsPoly = poly.Poly(nfsPoly)
	NF = poly.NumberField(nfsPoly)
	basis_q = find_lattice_basis(q,s)
	basis_q = reduce_basis(basis_q[0],basis_q[1])
	(qa0,qb0,qa1,qb1) = frankefy(basis_q)
	logB = math.log(B)	
	lambd_a = 1.5
	lambd_r = 1.5
		
	sieve = [[0.0] * (J) for i in range(I)]
	print "sieving rational side..."
	for prime in rfBase:
		r = prime[0]
		p = prime[1]
		#print [r,p]
		basis_r = find_sublattice_basis(basis_q,p,r)
		
		if(basis_r == False):
			continue
			
		plog = math.log(p)
		basis_r = reduce_basis(basis_r[0],basis_r[1])
		sieve_quadrant(1,1,I,J,basis_r,sieve,plog)
		sieve_quadrant(1,-1,I,J,basis_r,sieve,plog)
		sieve_quadrant(-1,1,I,J,basis_r,sieve,plog)
		sieve_quadrant(-1,-1,I,J,basis_r,sieve,plog)

	
	for i in range(-I/2,I/2):
		for j in range(-J/2,J/2):
			ijgcd = abs(fractions.gcd(i,j))
			if(ijgcd > 1):
				continue
				
			a = qa0*i+qa1*j
			b = qb0*i+qb1*j
			if(a == 0 or b == 0):
				continue

			normRat = a+b*m
			elog_r = math.log(abs(normRat))
			sieveBound_r = elog_r-lambd_r*logB
			if ((sieve[i+I/2-1][j+J/2-1] < (sieveBound_r))): 
				sieve[i+I/2-1][j+J/2-1] = float('-inf')
					
	print "sieving algebraic side..."
	for prime in afBase:
		r = prime[0]
		p = prime[1]
		basis_r = find_sublattice_basis(basis_q,p,r)
		
		if(basis_r == False):
			continue
			
		plog = math.log(p)
		basis_r = reduce_basis(basis_r[0],basis_r[1])
		sieve_quadrant(1,1,I,J,basis_r,sieve,plog)
		sieve_quadrant(1,-1,I,J,basis_r,sieve,plog)
		sieve_quadrant(-1,1,I,J,basis_r,sieve,plog)
		sieve_quadrant(-1,-1,I,J,basis_r,sieve,plog)
	
		
	print "done sieving, doing trial division..."
	smooths = []
	count = 0
	candidates = 0
	candidatesprev = 0
	for i in range(-I/2,I/2):
		for j in range(-J/2,J/2):
			'''
			if((candidates - candidatesprev ) == 500):
				candidatesprev = candidates
				print "%s smooths found out of %s candidates so far. i=%s j=%s" % (count,candidates,i,j)	
			'''
			if (sieve[i+I/2-1][j+J/2-1]<0):
				continue
			ijgcd = abs(fractions.gcd(i,j))
			if(ijgcd != 1):
				continue
				
			a = qa0*i+qa1*j
			b = qb0*i+qb1*j
			if(a == 0 or b == 0):
				continue
				
			abgcd = fractions.gcd(a,b)
			if(abgcd != 1):
				continue
				
			normRat = a+b*m
			elog_r = math.log(abs(normRat))
			sieveBound_r = elog_r-lambd_r*logB
			
			if (sieve[i+I/2-1][j+J/2-1] > sieveBound_r): 		
			
				normAlg = NF(poly.Poly([a,b])).norm()/q
				elog_a = math.log(abs(normAlg))
				sieveBound_a = elog_a-lambd_a*logB
				
				if (sieve[i+I/2-1][j+J/2-1] > sieveBound_a+sieveBound_r): 		
					candidates += 1
					if(etcmath.trialdivide(normAlg,afBase) and etcmath.trialdivide(normRat,rfBase)):		
						count += 1
						smooths.append([a,b])
					
	print "%s smooths found out of %s candidates." % (count,candidates)
	return smooths
def main():
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	
	afBase = files.loadFileArray("afbase.txt")
	rfBase = files.loadFileArray("rfbase.txt")
	smoothsFile = open("smooths.txt", "w")
	smoothsCount = 0
	q = B
	K+=K/10
	
	specialqFile = open("specialq.txt", "w")
	while(smoothsCount < K):
		while(True):
			q = primemath.nextPrime(q)
			roots = poly.getRootsModPFast(poly.Poly(nfsPoly),q)
			if(len(roots) > 0):
				break
			
		roots.sort()
		for s in roots:
			specialqFile.write(str([s,q])+"\n")
			print "special_q: %s" % [q,s]
			
			I = 1024
			J = 512
			smooths = sieve_special_q(q,s,I,J,rfBase,afBase,(n,nfsPoly,m,B,M,K))	
			smoothsCount += len(smooths)
			print "total smooths: %s/%s" % (smoothsCount,K)
			smoothsFile = open("smooths.txt", "a")
			for smooth in smooths:
				smoothsFile.write(str(smooth)+"\n")
					
			smoothsFile.close()	
if __name__ == '__main__':
	main()