from prime import *
import math

def multAllElements(one,array):
	prod = one
	for entry in array:
		prod = prod * entry
		
	return prod
	

def randomPrimes(totalsize,size):
	primes = [] 
	sum = 0
	while(sum <= totalsize):
		p = generateLargePrime(size)
		primes.append(p)
		sum += math.log(p,2)
		
	return primes
	
def calcRequiredPrimeLength(n, m, nfspoly, tomult):
	logest = math.log(nfspoly.degree(),2)*(nfspoly.degree()+5)*.5;
	logest += math.log(n,2)
	maxu = tomult[0][0]
	for smooth in tomult:
		if(maxu < abs(smooth[0])):
			maxu = abs(smooth[0])
		if(maxu < abs(smooth[1])):
			maxu = abs(smooth[1])
	b = 2*maxu*math.sqrt(nfspoly.degree())*m
	logest+= len(tomult)*.5*math.log(b,2)
	return logest
	
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
		
def getsr(p,e):
	q = p**e - 1
	counter = 0
	while (q%2 == 0):
		counter += 1
		q /= 2
	return (q,counter)
	
def expansionBaseM(n,m):
    q = n
    k = 0
    a = []
    i = len(str(n))

    while q != 0:
        a.append(q % m)
        q = q / m
        k += 1

    return a

