from scrollable_text import *
from propensity import *
from sympy import *
import mglobals as globals2
from recalculate_globals import *
from func_timeout import func_timeout, FunctionTimedOut 
import time 
  
def solve_with_timeout(Fe,js):
	sol = None
	sol = solve(Fe,js)
	return sol

def get_sets(Rr,Rp):
	tosum = set()
	sets = [tosum]
	used = set()
	for ih in range(len(Rr)):
		used2 = set()
		if len(Rr[ih]) == 1:
			for s in sets:
				key = list(Rr[ih].keys())[0]
				if key not in used and Rr[ih][key]!=0:
					s.add(key)
					used2.add(key)
		else:
			for j in range(len(sets)):
				s = sets[j]
				sets.append(s.copy())
				key = list(Rr[ih].keys())[0]
				if key not in used and Rr[ih][key]!=0:
					s.add(key)
					used2.add(key)
				key = list(Rr[ih].keys())[1]
				if key not in used and Rr[ih][key]!=0:
					sets[-1].add(key)
					used2.add(key)   
		
		if len(Rp[ih]) == 1:
			for s in sets:
				key = list(Rp[ih].keys())[0]
				if key not in used and Rp[ih][key]!=0:
					s.add(key)
					used2.add(key)
		else:
			for j in range(len(sets)):
				s = sets[j]
				sets.append(s.copy())
				key = list(Rp[ih].keys())[0]
				if key not in used and Rp[ih][key]!=0:
					s.add(key)
					used2.add(key)
				key = list(Rp[ih].keys())[1]
				if key not in used and Rp[ih][key]!=0:
					sets[-1].add(key)
					used2.add(key)	
		for zs in used2:
			used.add(zs)	
	return sets

def grab_steady_state(Rr,Rp,Cs,Cso,dA_dt,js):
	sets = get_sets(Rr,Rp)
	Fez = [dA_dt]
	for s in sets:
		Fez.append(sum([Cs[x]-Cso[x] for x in s]))
	val2 = solve(Fez,{x for x in js})
	xs = Matrix([val2[js[ih]] for ih in range(len(js))])
	return xs
	
def CstoCsR(fx,Cs,CsR,not_semi):
	if not_semi:
		rr = []
		for x in Cs:
			rr.append((Cs[x],CsR[x]))
		return fx.subs(rr)
	else:
		return fx
	
def Analyt_soln(Sp,Ks,conc,Rr,Rp,V,items=None,rfile="",not_semi=True,mode=None): 
	get_globals(rfile)
	MyEquations = {}

	if items:
		text = prepare_scroll_text(items)
		ffprint = lambda x: text.insert(INSERT," ".join([str(y) for y in x]))
	elif items == 0:
		ffprint = lambda x: print("",end="")
	else:
		ffprint = lambda x: print(" ".join([str(y) for y in x]),end="")	
		
	if not_semi:
		ffprint(["\n\nComplex Analytical expressions\n"])
		ffprint(["\n\nThe complex expression is because sympy do not know how you want to simplify the expression\n\n"])
	else:
		ffprint(["Simple semi-analytical expression\n\n"])
		
	Cs = {}
	Cso = {}
	t    = Symbol('t', negative=False,real=True)
	
	equivals = [(t,100)]
	CsR = {}
	for x in Sp:
		Cs[x] = Function(x,negative=False,real=True)(t)
		CsR[x] = Symbol(x, negative=False,real=True)
		if mode == "ftxo":
			Cso[x] = Symbol(x+"o", negative=False,real=True)
		elif mode == "ftks":
			Cso[x] = conc[x] 
		else:
			Cso[x] = Symbol(x+"o", negative=False,real=True) if not_semi else conc[x] 
		equivals.append((Cso[x],conc[x]))
	
	KCs = []
	for i in range(len(Ks)):
		row = []
		if len(Ks[i])==1:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, negative=False,real=True))
			equivals.append((row[0],Ks[i][0]))
		else:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, negative=False,real=True))
			equivals.append((row[0],Ks[i][0]))
			key = 'kb'+str(i+1)
			row.append(Symbol(key, negative=False,real=True))  
			equivals.append((row[1],Ks[i][1]))
		KCs.append(row)  
	
	if mode == "ftks":
		f = Matrix(propensity_vec_molar(KCs,Cs,Rr,Rp,True))
	elif mode == "ftxo":
		f = Matrix(propensity_vec_molar(Ks,Cs,Rr,Rp,True))
	else:
		f = Matrix(propensity_vec_molar(KCs,Cs,Rr,Rp,True)) if not_semi else Matrix(propensity_vec_molar(Ks,Cs,Rr,Rp,True))
		
		
	S = Matrix(V)
	#Cs might have change after call
	for x in Sp:
		Cs[x] = Function(x,negative=False,real=True)(t)
	Si = [x for x in Cs]	
	stoich = lambda x : max(1,abs(V[Si.index(x)][0]))

	Ss = []
	nz = []
	for row in range(S.shape[0]):
		if sum(abs(S[row,:])) != 0 and Si[row][0]!="-":
			Ss.append(list(S[row,:]))
			nz.append(row)	   
	Ss = Matrix(Ss)

	dA_dt = Ss*f
	ccs = [Cs[x] for x in Cs]
	ccso = [Cso[x] for x in Cso]
	js = [ ccs[x] for x in nz ]		
	jso = [ ccso[x] for x in nz ]		

	dA_dt2 = []
	for x in js:
		dA_dt2.append(x.diff(t))
	
	Fe = [Eq(dA_dt2[i],dA_dt[i]) for i in range(len(dA_dt2)) ]
	
	LHS = dA_dt.jacobian(js)
	#print(LHS)
	RHS = dA_dt-LHS*Matrix(js)	
	#print(RHS)
	#print(Fe)
	#print(dA_dt)
	#print(js)
	
	check = True
	work = False
	for i in range(LHS.shape[0]):
		if sum(abs(LHS.col(i))) == 0:
			check = False

	if check:
		print("Inside Check")
		LHSj = np.array(LHS)
		#print(LHSj)
		LHSj = Matrix(LHSj.flatten()).jacobian(js)	
		#print(LHSj)
		if LHSj != zeros(LHSj.shape[0],LHSj.shape[1]):	
			sets = get_sets(Rr,Rp)
			Fe2 = []
			for s in sets:
				Fe2.append(Eq(sum([Cs[x]/stoich(x)-Cso[x]/stoich(x) for x in s]),0))	
			xret = solve(Fe2,js)
			#print(xret)
			#print(Fe)
			
			for s in xret:
				#print(Fe[-1].atoms(Function))
				Fe[-1] = Eq(Fe[-1].lhs,Fe[-1].rhs.subs(s,xret[s]))	
			denom = CstoCsR(Fe[-1].rhs,Cs,CsR,not_semi)
			itvar = CstoCsR(js[-1],Cs,CsR,not_semi)
			#potAn = integrate(1/denom,(itvar,jso[-1],js[-1]))-t
			potAn = integrate(1/Fe[-1].rhs,(js[-1],jso[-1],js[-1]))-t
			Fe2.append(potAn)
			try: 
				sol = func_timeout( 20,solve, args=(Fe2,js) ) 
				for cc in sol:
					rval = Matrix(cc).subs(equivals)
					if np.sum(np.array(rval)<0)>0:
						pass
					else:
						x_sol = [cc[x] for x in range(len(js))]
						ind = 0
						for x in x_sol:
							varSp = str(js[ind])
							answer = str(simplify(x))
							ffprint([varSp," = ",answer,"\n\n"])	
							ind = ind + 1		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 
						break
				work = True
			except FunctionTimedOut: 
				work = False
		else:
			print("Attempting Matrix method")
			s    = Symbol('s',negative=False,real=True)
			xo = Matrix(jso)
			if sum(abs(RHS)) == 0:
				#print("llll",1)
				try:					
					xs = grab_steady_state(Rr,Rp,Cs,Cso,dA_dt,js)
					try:
						expD = func_timeout( 200,exp, args=(LHS*t,) )
						x_sol = xs + expD*(xo - xs)	
						ind = 0
						for x in x_sol:
							varSp = str(js[ind])
							answer = str(simplify(x))
							ffprint([varSp," = ",answer,"\n\n"])	
							ind = ind + 1		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 						
						work = True						
					except FunctionTimedOut: 
						work = False
				except:
					pass
			else:
				#print("llll",2)
				try:
					try:
						xs = simplify(-LHS**-1*RHS)
						x_sol = xs + exp(LHS*t)*(xo - xs)
					except:
						x_sol = simplify(exp(LHS*t)*xo+integrate(exp(LHS*t)*RHS,(s,0,t)))
					ind = 0
					for x in x_sol:
						varSp = str(js[ind])
						answer = str(simplify(x))
						ffprint([varSp," = ",answer,"\n\n"])	
						ind = ind + 1		
						if varSp not in MyEquations:
							MyEquations[varSp] = answer 											
					work = True
				except:
					pass

	if not work:
		used3 = {}
		try:
			print(1)				
			try:
				#sol = dsolve(Fe,js,ics={js[j].subs(t,0):jso[j] for j in range(len(js))})
				sol = func_timeout( 200, lambda : dsolve(Fe,js,ics={js[j].subs(t,0):jso[j] for j in range(len(js))}) )
				for x in range(len(sol)):				
					varSp = str(sol[x].lhs)
					answer = str(simplify(sol[x].rhs))
					ffprint(["\n",varSp," = ",answer,"\n\n"])		
					if varSp not in MyEquations:
						MyEquations[varSp] = answer 					
					used3[varSp] = True		
					
			except FunctionTimedOut: 
				float("jjj")			
		except:
			try:
				print(2)				
				ind = 0
				notF = []
				good = []
				for Fex in Fe:
					ics={js[ind].subs(t,0):jso[ind]}
					spc =  list(Fex.lhs.atoms(Function))
					spr =  list(Fex.rhs.atoms(Function))
					if len(spr) == 1 and spc[0] == spr[0]:
						#print(3)
						sol = [dsolve(Fex,ics=ics)]
						good.append([sol[0].lhs,sol[0].rhs])
						if str(sol[0].lhs) not in used3:													
							varSp = str(sol[0].lhs)
							answer = str(simplify(sol[0].rhs))
							ffprint(["\n",varSp," = ",answer,"\n\n"])		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 
							used3[varSp] = True						
					else:
						#print(4)
						#print([Fex,ics])
						notF.append([Fex,ics])
					ind = ind + 1
				
				#print(good)				
				#print(notF)
				ffprint(["\n\n","Attempting to evaluate","\n\n"])
				flen = 2*len(notF)
				cc = 0
				
				sets = get_sets(Rr,Rp)
				Fez = []
				#print(Si,V)
				for s in sets:
					Fez.append(Eq(sum([Cs[x]/stoich(x)-Cso[x]/stoich(x) for x in s]),0))	
				xret = solve(Fez,js)	
				#print(xret)
				
				while len(notF)>0 and cc<flen:
					x = notF.pop(0)
					#print(x,100000000000000000000000000000000000)
					spc =  list(x[0].lhs.atoms(Function))
					spr =  list(x[0].rhs.atoms(Function))
					for s in xret:
						#print(s not in spc)
						if s not in spc:
							subx = Eq(Symbol("Something_not_needed"),xret[s]).rhs.atoms(Function)
							if spc[0] in subx:
								#print(x[0],x[1],111)
								x[0] = x[0].subs(s,xret[s])	
								#print(x[0],x[1],222)
					#print(x,100000000000000000000000000000000000)			
					spr =  list(x[0].rhs.atoms(Function))
					terms = len(Add.make_args(x[0].rhs))
					if len(spr) == 1 and spc[0] == spr[0]:
						try:
							#print(x[0],x[1],333)
							sol = [dsolve(x[0],ics=x[1])]
							Fez.append(Eq(sol[0].lhs[0]-sol[0].rhs,0))
							good.append([sol[0].lhs,sol[0].rhs])
							if str(sol[0].lhs) not in used3:						
								varSp = str(sol[0].lhs)
								answer = str(simplify(sol[0].rhs))
								ffprint(["\n",varSp," = ",answer,"\n\n"])		
								if varSp not in MyEquations:
									MyEquations[varSp] = answer 
								used3[varSp] = True								
						except:
							denom = CstoCsR(x[0].rhs,Cs,CsR,not_semi)
							itvar = CstoCsR(spc[0],Cs,CsR,not_semi)
							val = integrate(1/denom,(itvar,x[1][spc[0].subs(t,0)],spc[0]))-t
							Fez.append(val)
							sol = solve(val,spc[0])
							if type(sol) == dict:
								sol = simplify(sol[spc[0]])
							else:
								sol = simplify(sol[0])
							good.append([spc[0],sol])						
							varSp = str(spc[0])
							answer = str(sol)
							ffprint(["\n",varSp," = ",answer,"\n\n"])		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 
							used3[varSp] = True																				
					elif len(spr) == 2 and terms == 1 and (len(Rr) == 1 or len(Rp) == 1) and not_semi:
						subt = Eq(spr[0]-spr[1],jso[js.index(spr[0])]-jso[js.index(spr[1])])
						if spr[0]!=spc[0]:
							vari = spr[0]
						else:
							vari = spr[1]
						sol = solve(subt,vari)
						if type(sol) == dict:
							sol = sol[vari]
						else:
							sol = sol[0]
						subx = Eq(Symbol("Something_not_needed"),sol).rhs.atoms(Function)
						if spc[0] in subx:			
							#print(x[0],x[1],444)
							x[0] = x[0].subs(vari,sol)	
							#print(x[0],x[1],555)
						notF.append([x[0],x[1]])			
					else:
						notF.append([x[0],x[1]])
					cc = cc + 1
				
				if len(notF)>0 and len(used3)<len(js):
					ffprint(["\n","The following are difficult to solve using dsolve","\n\n"])
					for x in notF:
						ffprint(["\n",x[0],"\n\n"])
					
					ffprint(["\n","The following expression may be wrong , please compare with ODE_int result :","\n\n"])			
					cc = 0		
					while len(notF)>0 and cc<flen:
						x = notF.pop(0)
						#print(x[0])
						for s in xret:
							if s not in x[0].lhs.atoms(Function):
								x[0] = x[0].subs(s,xret[s])	
						#print(x[0])
						
						#for y in good:
							#if y[0] not in x[0].lhs.atoms(Function):
								#x[0] = x[0].subs(y[0],y[1])	
						try:				
							spc =  list(x[0].lhs.atoms(Function))
							spr =  list(x[0].rhs.atoms(Function))
							if len(spr) == 1 and spc[0] == spr[0]:
								spc = spc[0]							
								try:
									denom = CstoCsR(x[0].rhs,Cs,CsR,not_semi)
									itvar = CstoCsR(spc,Cs,CsR,not_semi)
									val = integrate(1/denom,(itvar,x[1][spc.subs(t,0)],spc))-t
								except:
									denom = CstoCsR(x[0].rhs,Cs,CsR,False)
									itvar = CstoCsR(spc,Cs,CsR,False)
									val = integrate(1/denom,(itvar,x[1][spc.subs(t,0)],spc))-t							
															
								ffprint(["\n","Not simplified solution","\n\n"])
								ffprint(["need to solve for ",spc," to simplify :\n",val,"\n\n"])
								Fez.append(val)
								sol = solve(val,spc)
								if type(sol) == dict:
									sol = sol[spc]
								else:
									sol = simplify(sol[0])
								good.append([spc,sol])
								#varSp = str(spc)
								#answer = str(sol)
								#ffprint(["\n",varSp," = ",answer,"\n\n"])		
								#if varSp not in MyEquations:
									#MyEquations[varSp] = answer 	
							else:
								notF.append([x[0],x[1]])	
						except:
							notF.append([x[0],x[1]])
						cc = cc + 1
					
					ffprint(["\n","Attempting to simplify","\n\n"])		
					#print(Fez)
					sol = solve(Fez,js)
					for cc in sol:
						rval = Matrix(cc).subs(equivals)
						if np.sum(np.array(rval)<0)>0:
							pass
						else:
							x_sol = [cc[x] for x in range(len(js))]
							ind = 0
							for x in x_sol:							
								varSp = str(js[ind])
								answer = str(simplify(x))
								ffprint([varSp," = ",answer,"\n\n"])		
								if varSp not in MyEquations:
									MyEquations[varSp] = answer 									
								ind = ind + 1					
							break
			except:
				print(3)
				ind = 0
				notF = []
				good = []
				for Fex in Fe:
					ics={js[ind].subs(t,0):jso[ind]}
					spc =  list(Fex.lhs.atoms(Function))
					spr =  list(Fex.rhs.atoms(Function))
					if len(spr) == 1 and spc[0] == spr[0]:
						#print(3)
						sol = [dsolve(Fex,ics=ics)]
						good.append([sol[0].lhs,sol[0].rhs])
						if str(sol[0].lhs) not in used3:							
							varSp = str(sol[0].lhs)
							answer = str(sol[0].rhs)
							ffprint(["\n",varSp," = ",answer,"\n\n"])		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 
							used3[varSp] = True						
					else:
						#print(4)
						#print([Fex,ics])
						notF.append([Fex,ics])
					ind = ind + 1
					
				#print(notF,1111111111111111111)
				ffprint(["\n\n","Attempting to evaluate","\n\n"])
				flen = 2*len(notF)
				cc = 0
				
				sets = get_sets(Rr,Rp)
				Fez = []
				#print(Si,V)
				
				while len(notF)>0 and cc<flen:
					x = notF.pop(0)
					#print(x,100000000000000000000000000000000000)
					spc =  list(x[0].lhs.atoms(Function))							
					spr =  list(x[0].rhs.atoms(Function))
					terms = len(Add.make_args(x[0].rhs))
					rightOK = True
					for ss in spr:
						if ss in js:
							rightOK = False
							break
					#print(x,100000000000000000000000000000000000)		
					if (len(spr) == 1 and spc[0] == spr[0]) or rightOK:
						try:
							#print(x[0],x[1],333)
							sol = [dsolve(x[0],ics=x[1])]
							Fez.append(Eq(sol[0].lhs[0]-sol[0].rhs,0))
							good.append([sol[0].lhs,sol[0].rhs])
							if str(sol[0].lhs) not in used3:						
								varSp = str(sol[0].lhs)
								answer = str(sol[0].rhs)
								ffprint(["\n",varSp," = ",answer,"\n\n"])		
								if varSp not in MyEquations:
									MyEquations[varSp] = answer 
								used3[varSp] = True								
						except:
							denom = CstoCsR(x[0].rhs,Cs,CsR,not_semi)
							itvar = CstoCsR(spc[0],Cs,CsR,not_semi)
							val = integrate(1/denom,(itvar,x[1][spc[0].subs(t,0)],spc[0]))-t
							Fez.append(val)
							sol = solve(val,spc[0])
							if type(sol) == dict:
								sol = simplify(sol[spc[0]])
							else:
								sol = simplify(sol[0])
							good.append([spc[0],sol])
							varSp = str(spc[0])
							answer = str(sol)
							ffprint(["\n",varSp," = ",answer,"\n\n"])		
							if varSp not in MyEquations:
								MyEquations[varSp] = answer  
							used3[varSp] = True																			
					else:
						notF.append([x[0],x[1]])
					cc = cc + 1
				
				#print(notF)
				if len(notF)>0 and len(used3)<len(js):
					ffprint(["\n","The following are difficult to solve using dsolve","\n\n"])
					for x in notF:
						ffprint(["\n",x[0],"\n\n"])
					
					ffprint(["\n","The following expression may be wrong , please compare with ODE_int result :","\n\n"])			
					cc = 0		
					while len(notF)>0 and cc<flen:
						x = notF.pop(0)
						#print(x[0],1111)
						for y in good:
							if y[0] not in x[0].lhs.atoms(Function):
								x[0] = x[0].subs(y[0],y[1])	
						#print(x[0],2222)
						try:				
							spc =  list(x[0].lhs.atoms(Function))
							spr =  list(x[0].rhs.atoms(Function))
							rightOK = True
							for ss in spr:
								if ss in js:
									rightOK = False
									break
							
							#print((len(spr) == 1 and spc[0] == spr[0]) or rightOK, (len(spr) == 1 and spc[0] == spr[0]),rightOK)						
							if (len(spr) == 1 and spc[0] == spr[0]) or rightOK:
								#print(111111111)
								spc = spc[0]
								try:
									denom = CstoCsR(x[0].rhs,Cs,CsR,not_semi)
									itvar = CstoCsR(spc,Cs,CsR,not_semi)
									#print(1/denom,(itvar,x[1][spc.subs(t,0)],spc),1111)
									val = integrate(1/denom,(itvar,x[1][spc.subs(t,0)],spc))-t
								except:
									denom = CstoCsR(x[0].rhs,Cs,CsR,False)
									itvar = CstoCsR(spc,Cs,CsR,False)
									#print(1/denom,(itvar,x[1][spc.subs(t,0)],spc),2222)
									val = integrate(1/denom,(itvar,x[1][spc.subs(t,0)],spc))-t								
									
								ffprint(["\n","Not simplified solution","\n\n"])
								ffprint(["need to solve for ",spc," to simplify :\n",val,"\n\n"])
								Fez.append(val)
								sol = solve(val,spc)
								if type(sol) == dict:
									sol = sol[spc]
								else:
									sol = simplify(sol[0])
								good.append([spc,sol])
								varSp = str(spc)
								answer = str(sol)
								ffprint(["\n",varSp," = ",answer,"\n\n"])		
								if varSp not in MyEquations:
									MyEquations[varSp] = answer 								
							else:
								notF.append([x[0],x[1]])	
						except:
							notF.append([x[0],x[1]])
						cc = cc + 1

					ffprint(["\n","Attempting to simplify","\n\n"])		
					#print(Fez)
					sol = solve(Fez,js)
					if type(sol) ==dict:
						for cc in sol:
							varSp = str(cc)
							answer = str(sol[cc])
							ffprint(["\n",varSp," = ",answer,"\n\n"])	
							if varSp not in MyEquations:
								MyEquations[varSp] = answer 
					else:
						for cc in sol:
							rval = Matrix(cc).subs(equivals)
							if np.sum(np.array(rval)<0)>0:
								pass
							else:
								x_sol = [cc[x] for x in range(len(js))]
								ind = 0
								for x in x_sol:							
									varSp = str(js[ind])
									answer = str(simplify(x))
									ffprint(["\n",varSp," = ",answer,"\n\n"])		
									if varSp not in MyEquations:
										MyEquations[varSp] = answer 									
									ind = ind + 1					
								break	
	#print(Si)
	return MyEquations
