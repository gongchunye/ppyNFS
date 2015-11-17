import poly
import primemath
import liblll
import etcmath


def reduce_basis(p1,p2):
	mat_basis = liblll.create_matrix([[p1[0],p2[0]],[p1[1],p2[1]]]) 
	mat_reduced = liblll.lll_reduction(mat_basis)	
	a0 = mat_reduced[0][0].numerator # lll library creates fractions with denominator 1, so grab numerator
	a1 = mat_reduced[0][1].numerator
	b0 = mat_reduced[1][0].numerator
	b1 = mat_reduced[1][1].numerator
	return [[a0,b0],[a1,b1]]
	
def chinese_rem(a,b,p,q):
	return (a*q*etcmath.modInv(q,p) + b*p*etcmath.modInv(p,q)) % (p*q)
	
def find_sublattice_basis(q,r,subp,subr):
	basis1 = []
	b=0
	while(True):
		b += 1
		a = chinese_rem(-b*r,-b*subr,q,subp)
		if(a == 0 or b == 0):
			continue			
		basis1 = [a,b]
		break
		
	while(True):
		b += 1
		a = chinese_rem(-b*r,-b*subr,q,subp)
		if(a == 0 or b == 0):
			continue		
		if(float(basis1[0])/a - float(basis1[1])/b < 0.01):
			continue		
		basis2 = [a,b]
		break

	return [basis1,basis2]
	
def main():
	nfspoly = poly.Poly([8,29,15,1])
	NF = poly.NumberField(nfspoly)
	while(True):
		q = primemath.generateLargePrime(20)
		r = poly.getRootsModPFast(nfspoly,q)
		if(len(r) > 0):
			break
	r = r[0]

	print "special_q: f(%s) = 0 mod %s" %(r,q)
		
	while(True):
		subp = primemath.generateLargePrime(7)
		subr = poly.getRootsModPFast(nfspoly,subp)
		if(len(subr) > 0):
			break
	subr = subr[0]

	print "subp: f(%s) = 0 mod %s" %(subr,subp)
		
	basis = find_sublattice_basis(q,r,subp,subr)
	print "basis: %s" % basis
	basis_r = reduce_basis(basis[0],basis[1])
	a0 = basis_r[0][0]
	b0 = basis_r[0][1]
	a1 = basis_r[1][0]
	b1 = basis_r[1][1]
	print "reduced basis: %s" % [[a0,b0],[a1,b1]]

	print "generating lattice points..."
	for i in range(-2,2):
		for j in range(-2,2):
			a = a0*i+a1*j
			b = b0*i+b1*j
			if(a == 0 or b == 0):
				continue
			NFelement = NF(poly.Poly([a,b]))
			if(NFelement.norm()%q == 0 and NFelement.norm()%subp == 0):
				print "%s is q-divisible and subp-divisible" % [a,b]
			else:
				raise AssertionError("%s is not q-divisible or subp-divisible" % [a,b])

if __name__ == '__main__':
	main()