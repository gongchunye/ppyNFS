from poly import *
from etc_math import *
from files import *
import fractions

if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	smooths = []
	
	rfBase = loadFileArray("rfbase.txt")
	
	a = -M
	b = -1
	while(len(smooths) < K):
	
		a += 1
		if(a > M):
			a = -M
			b -= 1
			
		if(fractions.gcd(abs(a),abs(b)) != 1):
			continue
		
		polyToTest = Poly([a,b])
		toTrialDivide = NF(polyToTest).norm()*(a+b*m)
		if(trialdivide(toTrialDivide,rfBase)):
			print "%s/%s: %s" % (len(smooths)+1,K,polyToTest)
			smooths.append(polyToTest)
			
	

			