import poly
import nfspolygen
import etcmath
import primemath
import math

if __name__ == '__main__':
	semiprimeSize = 64
	n = primemath.generateLargePrime(semiprimeSize/2)*primemath.generateLargePrime(semiprimeSize/2)	
	print "n = %s" % (n)
	
	d = 3
	(m,nfsPoly) = nfspolygen.generateNFSPoly(n,d)	
	print "Using poly %s, m = %s" % (nfsPoly,m)
	B = M = etcmath.calculateB(n)
	print "Using bound %s" % B
	prefactorBase = primemath.generatePrimes(B)
	
	K = (int)(3*math.log(n,10)) # for quadratic characters
	
	afBaseFile = open("afbase.txt", "w")
	rfBaseFile = open("rfbase.txt", "w")
	print "Generating af and rf bases..."	
	for p in prefactorBase:
		if(p > B): break
		rfBaseFile.write(str([m%p,p])+"\n")
		K += 1
		roots = poly.getRootsModPFast(nfsPoly,p)
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
	
				
