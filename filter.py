import poly
import files
import matrixmath

def calcRatMatrixRowsSum(rfBase,smooths,m,NF):
	matrixRowSum = [0] * (len(rfBase))
	smoothsToRemove = list(smooths)
	for smoothPair in smooths:
		smoothPoly = poly.Poly(smoothPair)
		matrixRow = matrixmath.getMatrixRowRat(smoothPoly,m,rfBase)
		matrixRow.pop(0) # remove sign bit

		for i in range(len(matrixRow)):
			matrixRowSum[i] += matrixRow[i]
			if(matrixRowSum[i] > 1 and smoothPair in smoothsToRemove):
				smoothsToRemove.remove(smoothPair)
	
	return matrixRowSum
	
def calcAlgMatrixRowsSum(afBase,smooths,m,NF):
	matrixRowSum = [0] * (len(afBase))
	
	for smoothPair in smooths:
		smoothPoly = poly.Poly(smoothPair)
		matrixRow = matrixmath.getMatrixRowAlg(smoothPoly,NF,afBase)
		for i in range(len(matrixRow)):
			matrixRowSum[i] += matrixRow[i]
	
	return matrixRowSum

def findRowValueFB(fb,row,n):
	subfb = []
	for i in range(len(row)):
		if(row[i] == n):
			subfb.append(fb[i])
			
	return subfb
	
def hash(smooth):
	return 314159265358979323*smooth[0]+ 271828182845904523*smooth[1] % 2**64
	
def main():
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	NF = poly.NumberField(poly.Poly(nfsPoly))
	
	rfBase = files.loadFileArray("rfbase.txt")
	afBase = files.loadFileArray("afbase.txt")
	smooths = files.loadFileArray("smooths.txt")
	specialq = files.loadFileArray("specialq.txt")
	
	nOfSmooths = len(smooths)

	
	#smooths = sorted(smooths,key=hash)
	#smooths = sorted({hash(val): val for val in smooths}, key=hash)
	smooths = sorted({hash(val): val for val in smooths}.values(), key=hash)
	#for smooth in smooths:
		#print smooth
		
	print "%s duplicates removed."% (nOfSmooths-len(smooths))
	nOfSmooths = len(smooths)
	
	if(len(smooths) < K+len(specialq)):
		print "undersieved by %s." % ((K+len(specialq)) - len(smooths))
		raise AssertionError
	
	while(len(smooths) > K+len(specialq)):
		smooths.pop()
		
	print "%s extra relations removed."% (nOfSmooths-len(smooths))
	
	smoothsFile = open("smooths-fil.txt", "w")
	for smooth in smooths:
		smoothsFile.write(str(smooth)+"\n")
				
	rfBaseFile = open("rfbase-fil.txt", "w")
	for p in rfBase:
		rfBaseFile.write(str(p)+"\n")				
				
	afBaseFile = open("afbase-fil.txt", "w")
	for p in afBase:
		afBaseFile.write(str(p)+"\n")
	for p in specialq:
		afBaseFile.write(str(p)+"\n")
				

	'''
	matrixRowSum = calcRatMatrixRowsSum(rfBase,smooths,m,NF)
	print matrixRowSum
	print "\n"
	print findRowValueFB(rfBase,matrixRowSum,0)
	
	found = True
	while(found = True):
		found = False
		matrixRatRowSum = calcRatMatrixRowsSum(rfBase,smooths,m,NF)
		findRowValueIndices(matrixRatRowSum,0)
		findRowValueIndices(matrixRatRowSum,1)
	'''
	
	
if __name__ == '__main__':
	main()	

		