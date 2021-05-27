from sympy import *
import numpy as np

done = set()
def process(x):
	global done
	xx =  x.strip("*").strip().replace("*/","/")
	if xx[0] == "/":
		xx = "1"+xx
	if xx[-1] == "/":
		xx = xx.strip("/")
	val = sympify(xx)
	if val in done:
		return None
	else:
		done.add(val)
	return val

def propExt(expr,prop):
	ex = str(expr)
	open = 0
	close = 0
	diff = 0
	
	collect = ""
	for v in ex:
		if v in ["+","-"] and diff == 0:
			d = process(collect)
			if d != None:
				prop.append(d)
			collect = ""
			open = 0
			close = 0
		elif v == "(":
			open = open + 1
			collect = collect + v
		elif v == ")":
			close = close + 1
			collect = collect + v
		elif (v.isnumeric() or v == ".") and (diff == 0):
			if len(collect)>=4:
				if collect[-2:] == "**" or collect[-1] not in ["*","/"]:
					collect = collect + v
			elif len(collect)>=1:
				if collect[-1] not in ["*","/"]:
					collect = collect + v
			else:
				pass
		else:
			collect = collect + v	
		diff = open - close

	d = process(collect)
	if d != None:
		prop.append(d)
		
def termExt(expr):
	term = []
	ex = str(expr)
	open = 0
	close = 0
	diff = 0
	
	collect = ""
	for v in ex:
		if v in ["+","-"] and diff == 0:
			term.append(collect)
			collect = ""
			open = 0
			close = 0
		elif v == "(":
			open = open + 1
			collect = collect + v
		elif v == ")":
			close = close + 1
			collect = collect + v
		else:
			collect = collect + v	
		diff = open - close

	term.append(collect)
	return term
		
def get_prop_stoich(dxdt):
	prop = []
	dAdt = []
	for expr in dxdt:
		propExt(expr,prop)
		dAdt.append(expr)
		
	w = prop
	V = [[0 for x in range(len(w))] for y in range(len(dAdt))]
	
	for i in range(len(dAdt)):
		for j in range(len(w)):
			s = 0
			for xx in termExt(dAdt[i]):
				x = sympify(xx)/w[j]
				try:
					s = s + float(x)
				except:
					pass
			V[i][j] = s
	
	return Matrix(V), Matrix(w)

"""
dxdt = [
	"2*A*B+k1*A/(2+B)+C*(D+A)/(G+W)",
	"24*A*B+k1*A*3/(2+B)+C*(D+A)/(G+W)",
	"24*A*B+k1*A*3/(2+B)+0.5*C*(D+A)/(G+W)",
	"24*A*B+k1*A**2/(2+B**2)+0.5*C*(D+A)/(G+W)",
	"24*A*B+k1*24/(2+B**2)+0.5*C*(D+A)/(G+W)",
	"24*A*B+24/(2+B**2)/2+0.5*C*(D+A)/(G+W)",
	"24*A*B+h123*245/(2+B**2)/2+0.5*C*(D+A)/(G+W)"
]

V, w = get_prop_stoich(dxdt)
for t in V*w:
	print(t)
print()

for c in np.array(V):
	print(c)
"""

