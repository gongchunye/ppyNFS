import poly
import primemath
import liblll
def main():
	nfspoly = poly.Poly([8,29,15,1])
	NF = poly.NumberField(nfspoly)
	q = 359
	r = 110
	if(nfspoly.evaluate(r)%q==0):
		print "r is a root of f() mod q"
	
	basis = []
	b = 1
	for i in range(2):
		a = -b*r + i*q
		basis.append([a,b])
	

	print "basis: %s" % basis
	mat_basis = liblll.create_matrix(basis) 
	mat_reduced = liblll.lll_reduction(mat_basis)
	#mat_reduced = mat_basis # uncomment to not reduce lattice basis
	a0 = mat_reduced[0][0].numerator # lll library creates fractions with denominator 1, so grab numerator
	b0 = mat_reduced[0][1].numerator
	a1 = mat_reduced[1][0].numerator
	b1 = mat_reduced[1][1].numerator
	print "reduced basis: %s" % [[a0,b0],[a1,b1]]
	
	for i in range(-2,2):
		for j in range(-2,2):
			a = a0*i+a1*j
			b = b0*i+b1*j
			if(a == 0 or b == 0):
				continue
			NFelement = NF(poly.Poly([a,b]))
			print [a,b]
			if(NFelement.norm()%q == 0):
				print "q-divisible"

if __name__ == '__main__':
	main()