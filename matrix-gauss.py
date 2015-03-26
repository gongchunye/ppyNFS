from poly import *
from etc_math import *
from files import *

def getMatrixRowRat(smoothPoly,m,rfBase):
	value = smoothPoly.evaluate(m)
	matrixRow = []
	for prime in rfBase:
		ctr = 0
		while(value % prime == 0):
			value /= prime
			ctr = (ctr+1) % 2
		matrixRow.append(ctr)
		
	return matrixRow
	
def getMatrixRowAlg(smoothPoly,NF,afBase):
	smoothNorm = NF(smoothPoly).norm()
	matrixRow = []
	for prime in afBase:
		ctr = 0
		r = prime[0]
		p = prime[1]
		a = smoothPair[0]
		b = smoothPair[1]
		if(a%p == (-b*r) % p):
			while(smoothNorm % p == 0):
				ctr = (ctr+1) % 2
				smoothNorm /= p	
		matrixRow.append(ctr)
		
	return matrixRow
	
def getMatrixRowQC(smoothPoly,qcBase):
	matrixRow = []
	for qcPair in qcBase:
		r = qcPair[0]
		p = qcPair[1]
		c = legendre_symbol(smoothPoly.evaluate(r),p)
		if(c == 1):
			matrixRow.append(0)
		elif(c == -1):
			matrixRow.append(1)
			
	return matrixRow
			
def tranposeMatrix(matrix):
	return [list(i) for i in zip(*matrix)]
	
def swapRow(matrix,i,j):
	temp = matrix[i]
	matrix[i] = matrix[j]
	matrix[j] = temp
	
def addRow(matrix,i,j):
	for k in range(len(matrix[i])):
		matrix[i][k] = (matrix[i][k] + matrix[j][k]) % 2
	
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
							addRow(matrix,r2,r)
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
						swapRow(matrix,r,c)
						
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
	
def printMatrix(matrix):
	for row in matrix:
		print row
		
if __name__ == '__main__':
	(n,nfsPoly,m,B,M,K) = loadParamsFile()
	NF = NumberField(Poly(nfsPoly))
	
	
	rfBase = loadFileArray("rfbase.txt")
	afBase = loadFileArray("afbase.txt")
	qcBase = loadFileArray("qcbase.txt")
	
	smooths = loadFileArray("smooths.txt")
	print "Building matrix..."
	matrix = []
	for smoothPair in smooths:
		smoothPoly = Poly(smoothPair)
		matrixRow = getMatrixRowRat(smoothPoly,m,rfBase)
		matrixRow.extend(getMatrixRowAlg(smoothPoly,NF,afBase))
		matrixRow.extend(getMatrixRowQC(smoothPoly,qcBase))
		matrixRow.extend([0]*(K-len(rfBase)-len(afBase)-len(qcBase)))
		matrix.append(matrixRow)
	
	#print matrix
	print "Transposing matrix..."
	matrix = tranposeMatrix(matrix)
	print "Doing gaussian reduction..."
	gaussMatrix(matrix)
	
	deps = open("deps.txt", "w")
	
	print "Writing dependencies..."
	for i in range(16,32):
		dep = findDependency(matrix,i)
		deps.write(str(dep)+"\n")
	deps.close()
		
	
	
