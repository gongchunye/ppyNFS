# ppyNFS
Number Field Sieve on python

The Number Field Sieve (NFS) is the asymptotically fastest integer factorization algorithm known. It is magnitudes more complicated than its younger brother, the Quadratic Sieve. The fact that a very complicated algorithm solves a (conceptually) very simple problem is what give its appeal (similar to FLT). This algorithm requires LOTS of sub-algorithms, and some code is copied from the implementations scattered throughout the internet. These will eventually be rewritten.

The goal of this project is to write an easy-to-understand implementation to ease learning of NFS-related concepts. Python allows overriding of operators, supports bignums natively, and is moderately popular for beginners. Here is a snippet of code that shows readability (see tests.py):
		
		# Brigg's example
		nfspoly = Poly([8,29,15,1])
		NF = NumberField(nfspoly)
		
		tomult = [[-1,1],[3,1],[13,1],[104,1],[3,2],[25,2],[-8,3],[48,5],[54,5],[-43,6],[-8,7],[11,7],[856,11]]
		tomult = [NF(Poly(x)) for x in tomult] # Transform each (a,b) into an element in the number field.
		nfspolyd = NF(nfspoly.derivative())
		
		prod = NF(Poly([1]))
		for x in tomult: # Multiply each (a,b) element.
			prod = prod * x
			
		prod = prod * nfspolyd * nfspolyd
		
		correctProduct = NF(Poly([22939402657683071224L, 54100105785512562427L, 22455983949710645412L])) # correct product, as per Briggs'
		self.assertEqual(prod,correctProduct)
		
A complete (dirty) implementation of NFS in python is being rewritten to ppyNFS, with the aims of making the code cleaner. Though, I don't claim that this code is absolutely tidy (cleanup is still needed in lots of places). The polynomial and number field logic is mostly rewritten. The main NFS stages that has been rewritten are:

[X] init

[X] trialdivide

[ ] sieve-line

[ ] filter

[X] matrix-gauss

[X] sqrt-couv

Also, stages that were not in the old code but would be nice to have in ppyNFS are:

[ ] matrix-lanc

[X] sieve-lat

[ ] sqrt-nguy

The completed parts is enough to run a basic NFS factorization. To do this, run:
		
		init.py
		sieve-lat.py
		filter.py
		matrix-gauss.py
		sqrt-couv.py
		
LLL code is from https://github.com/kutio/liblll. Unless otherwise stated, the code is in public domain.
