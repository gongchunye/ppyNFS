import poly
import files
import matrixmath
import random
	
def gaussMatrix(matrix):
	# TODO: implement proper pivot algorithm.
	sideLen = len(matrix)
	found = True
	while(found):
		found = False
		for r in range(sideLen):
			for c in range(sideLen):
				if(matrix[r][c] == 1):
					for r2 in range(sideLen):
						if(r==r2): 
							continue
						if(matrix[r2][c] == 1):
							found = True
							matrixmath.addRow(matrix,r2,r)
					break
			
		
	found = True
	while(found):
		found = False
		for r in range(sideLen):
			for c in range(sideLen):
				if(matrix[r][c] == 1):
					if(c == r):
						break
					else:
						found = True
						matrixmath.swapRow(matrix,r,c)
						
def zeroRow(row):
	for i in range(len(row)):
		if(row[i] == 1):
			return False
	return True

def findDependency(matrix,n):
	sideLen = len(matrix)
	dep = [0]*sideLen
	for r in range(sideLen):
		if(zeroRow(matrix[r])):
			dep[r] = n & 1
			n >>= 1
			
	for r in range(sideLen):
		for c in range(r+1,sideLen):
			if(matrix[r][c] == 1):
				dep[r] = (dep[c] + dep[r]) % 2
	
	return dep

		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = files.loadParamsFile()
	NF = poly.NumberField(poly.Poly(nfsPoly))
	
	
	rfBase = files.loadFileArray("rfbase-fil.txt")
	rfBase = zip(*rfBase)[1] #grab only the primes
	afBase = files.loadFileArray("afbase-fil.txt")
	qcBase = matrixmath.generateQCBase(n,afBase[-1][1],poly.Poly(nfsPoly))
	
	smooths = files.loadFileArray("smooths-fil.txt")
	print "Building matrix..."
	matrix = []
	for smoothPair in smooths:
		smoothPoly = poly.Poly(smoothPair)
		matrixRow = matrixmath.getMatrixRowRat(smoothPoly,m,rfBase)
		matrixRow.extend(matrixmath.getMatrixRowAlg(smoothPoly,NF,afBase))
		matrixRow.extend(matrixmath.getMatrixRowQC(smoothPoly,qcBase))
		matrixRow.extend([0]*(K-len(rfBase)-len(afBase)-len(qcBase)-1)) # one extra for sign
		matrix.append(matrixRow)
	
	#matrixmath.printMatrix(matrix)
	print "Transposing matrix..."
	matrix = matrixmath.tranposeMatrix(matrix)
	print "Doing gaussian reduction..."
	gaussMatrix(matrix)
	#matrixmath.printMatrix(matrix)
	deps = open("deps.txt", "w")
	
	print "Writing dependencies..."
	for i in range(1,33):
		dep = findDependency(matrix,i)
		deps.write(str(dep)+"\n")
	deps.close()
		
	
	
