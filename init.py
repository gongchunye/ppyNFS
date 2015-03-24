from poly import *
from nfspolygen import *
from etc_math import *
from poly_math import *
from prime import *

if __name__ == '__main__':
	semiprimeSize = 50
	n = generateLargePrime(semiprimeSize/2)*generateLargePrime(semiprimeSize/2)
	print "n = %s" % (n)
	
	d = 3
	(m,nfsPoly) = generateNFSPoly(n,d)	
	print "Using poly %s, m = %s" % (nfsPoly,m)
	B = M = calculateB(n)
	print "Using bound %s" % B
	prefactorBase = primes(B*2)
	
	qcBaseFile = open("qcbase.txt", "w")
	K = 0
	print "Generating qc base..."
	while(K < 3*math.log(n,10)):
		p = prefactorBase.pop()
		roots = getRootsModPFast(nfsPoly,p)
		for root in roots:
			if(nfsPoly.derivative().evaluate(root) % p != 0):
				qcBaseFile.write(str([root,p])+"\n")
				K += 1
	qcBaseFile.close()
	
	afBaseFile = open("afbase.txt", "w")
	rfBaseFile = open("rfbase.txt", "w")
	print "Generating af and rf bases..."	
	for p in prefactorBase:
		if(p > B): break
		rfBaseFile.write(str(p)+"\n")
		K += 1
		roots = getRootsModPFast(nfsPoly,p)
		for root in roots:
			afBaseFile.write(str([root,p])+"\n")
			K += 1
			
	afBaseFile.close()
	rfBaseFile.close()
	
	K += 10

	paramsFile = open("params.txt", "w")
	paramsFile.write("n = " +str(n)+"\n")
	paramsFile.write("B = " +str(B)+"\n")
	paramsFile.write("M = " +str(M)+"\n")
	paramsFile.write("nfsPoly = " +str(nfsPoly)+"\n")
	paramsFile.write("m = " +str(m)+"\n")
	paramsFile.write("K = " +str(K)+"\n")
	paramsFile.close()
	
	print "Run sieve.py."
	
				
