import etcmath

def getMatrixRowRat(smoothPoly,m,rfBase):
	value = smoothPoly.evaluate(m)
	matrixRow = []
	if(value < 0):
		matrixRow.append(1)
	else:
		matrixRow.append(0)
	for prime in rfBase:
		ctr = 0
		while(value % prime == 0):
			value /= prime
			ctr = (ctr+1) % 2
		matrixRow.append(ctr)
	if(abs(value) != 1):
		print "unsmooth rational: %s" % value
		raise AssertionError
	return matrixRow
	
def getMatrixRowAlg(smoothPoly,NF,afBase):
	smoothNorm = NF(smoothPoly).norm()
	matrixRow = []
	for prime in afBase:
		ctr = 0
		r = prime[0]
		p = prime[1]
		a = smoothPoly.coeffs[0]
		b = smoothPoly.coeffs[1]
		if(a%p == (-b*r) % p):
			while(smoothNorm % p == 0):
				ctr = (ctr+1) % 2
				smoothNorm /= p	
		matrixRow.append(ctr)
		
	if(abs(smoothNorm) != 1):
		print "unsmooth norm: %s" % smoothNorm
		raise AssertionError
	return matrixRow

def getMatrixRowQC(smoothPoly,qcBase):
	matrixRow = []
	for qcPair in qcBase:
		r = qcPair[0]
		p = qcPair[1]
		c = etcmath.legendre_symbol(smoothPoly.evaluate(r),p)
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
		
def printMatrix(matrix):
	for row in matrix:
		print row		