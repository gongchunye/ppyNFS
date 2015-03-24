from poly import *
from etc_math import *
from files import *
import fractions

if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	
	rfBase = loadFileArray("rfbase.txt")
	
	smoothsFile = open("smooths.txt", "w")
	a = -M
	b = -1
	smoothsCount = 0
	while(smoothsCount < K):
	
		a += 1
		if(a > M):
			a = -M
			b -= 1
			
		if(fractions.gcd(abs(a),abs(b)) != 1):
			continue
		
		polyToTest = Poly([a,b])
		toTrialDivide = NF(polyToTest).norm()*(a+b*m)
		if(trialdivide(toTrialDivide,rfBase)):
			smoothsCount += 1
			print "%s/%s: %s" % (smoothsCount,K,polyToTest)
			smoothsFile.write(str([a,b])+"\n")
	
	
	smoothsFile.close()
			