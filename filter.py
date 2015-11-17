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
	
def main():
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	NF = poly.NumberField(poly.Poly(nfsPoly))
	
	rfBase = files.loadFileArray("rfbase.txt")
	afBase = files.loadFileArray("afbase.txt")
	smooths = files.loadFileArray("smooths.txt")
	
	
	
	
	matrixRowSum = calcRatMatrixRowsSum(rfBase,smooths,m,NF)
	print matrixRowSum
	print "\n"
	print findRowValueFB(rfBase,matrixRowSum,0)
	'''
	found = True
	while(found = True):
		found = False
		matrixRatRowSum = calcRatMatrixRowsSum(rfBase,smooths,m,NF)
		findRowValueIndices(matrixRatRowSum,0)
		findRowValueIndices(matrixRatRowSum,1)
	'''
	
	
if __name__ == '__main__':
	main()	

		