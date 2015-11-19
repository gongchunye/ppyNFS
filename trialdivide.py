import poly
import etcmath
import files
import fractions

if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	NF = poly.NumberField(poly.Poly(nfsPoly))
	
	rfBase = files.loadFileArray("rfbase.txt")
	afBase = files.loadFileArray("afbase.txt")
	open("specialq.txt", "w")
	
	smoothsFile = open("smooths.txt", "w")
	a = -M
	b = 1
	smoothsCount = 0
	while(smoothsCount < K):
	
		a += 1
		if(a > M):
			a = -M
			b += 1
			
		if(fractions.gcd(a,b) != 1):
			continue
		
		polyToTest = poly.Poly([a,b])
		if(etcmath.trialdivide(NF(polyToTest).norm(),afBase) and etcmath.trialdivide((a+b*m),rfBase)):
			smoothsCount += 1
			print "%s/%s: %s" % (smoothsCount,K,polyToTest)
			smoothsFile.write(str([a,b])+"\n")
	
	
	smoothsFile.close()
			