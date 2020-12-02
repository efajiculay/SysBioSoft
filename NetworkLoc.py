from scrollable_text import *
from propensity import *
from sympy import *

def subs2(Z,cval):
	for x in cval:
		Z = Z.subs(x,cval[x])
	return Z
	
def NetLoc_symbolic(Sp,Ks,conc,Rr,Rp,V,items=None,molar=False,mode=None,numer=False):

	if items:
		text = prepare_scroll_text(items)
		ffprint = lambda x: text.insert(INSERT," ".join([str(y) for y in x]))
	else:
		ffprint = lambda x: print(" ".join([str(y) for y in x]),end="")	
	
	Cs = {}
	Cso = {}
	equivals = []
	equiCo = []
	equiKs = []
	
	for x in Sp:
		Cs[x] = Symbol(x, real = True, negative=False)
		Cso[x] = Symbol(x+'o', real = True, negative=False)*(0 if conc[x] == 0 else 1)		
		equivals.append((Cso[x],conc[x]))
		equiCo.append((Cso[x],conc[x]))
	
	KCs = []
	for i in range(len(Ks)):
		row = []
		if len(Ks[i])==1:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, real = True, negative=False))
			equivals.append((row[0],Ks[i][0]))
			equiKs.append((row[0],Ks[i][0]))
		else:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, real = True, negative=False))
			equivals.append((row[0],Ks[i][0]))
			equiKs.append((row[0],Ks[i][0]))
			key = 'kb'+str(i+1)
			row.append(Symbol(key, real = True, negative=False))
			equivals.append((row[1],Ks[i][1]))	
			equiKs.append((row[1],Ks[i][1]))			
		KCs.append(row)  
		
	if not molar:
		f = Matrix(propensity_vec(KCs,Cs,Rr,Rp))
	else:
		f = Matrix(propensity_vec_molar(KCs,Cs,Rr,Rp))
	
	if mode=="Numeric":
		f = f.subs(equivals)
	elif mode == "fofks":
		f = f.subs(equiCo)
	elif mode == "fofCo":
		f = f.subs(equiKs)
		
	S = Matrix(V)
	#Cs might have change after call
	for x in Sp:
		Cs[x] = Symbol(x, real = True, negative=False)
	Si = [x for x in Cs]

	Ss = []
	nz = []
	for row in range(S.shape[0]):
		if sum(abs(S[row,:])) != 0 and Si[row][0]!="-":
			Ss.append(list(S[row,:]))
			nz.append(row)	   
	Ss = Matrix(Ss)
	  
	ccs = [Cs[x] for x in Cs]
	js = [ ccs[x] for x in nz ]		
	A = f.jacobian(js)
	kerV = Ss.nullspace()
	
	As = A.col_insert(A.shape[0],-kerV[0])
	for ih in range(1, len(kerV)):
		As = As.col_insert(As.shape[0],-kerV[ih])


	KNAMES = []
	for i in range(len(KCs)):
		if len(KCs[i])==1:
			KNAMES.append(KCs[i][0])
		else:
			KNAMES.append(KCs[i][0])
			KNAMES.append(KCs[i][1])
			
	ffprint(["When k decrease, the sign of the sensitivity tells what happen to species activity\n"])
	ffprint(["When k increase, the reverse of the sign of the sensitivity tells what happen to species activity\n\n"])
	
	if numer:
		As = As.subs(equivals)

	try:
		S = simplify(As.inv()).T
		if numer:
			S = S.subs(equivals)
		for ih in range(S.shape[0]):
			for ij in range(S.shape[1]):
				try:
					ffprint([KNAMES[ih],"\t",js[ij],"\t",S[ih,ij],"\n"])		
				except:
					pass	
			ffprint(["\n"])
	except NonSquareMatrixError:
		ffprint(["Non square matrix error"])
	finally:
		pass

	return [0,0]