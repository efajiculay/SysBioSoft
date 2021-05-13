import sympy

#def abs(x):
	#return sympy.fabs(x)
	
def acos(x):
	return sympy.acos(x).evalf()
	
def arccos(x):
	return sympy.acos(x).evalf()
	
def acosh(x):
	return sympy.acosh(x).evalf()
	
def arccosh(x):
	return sympy.acosh(x).evalf()
	
def acot(x):
	return sympy.acot(x).evalf()
	
def arccot(x):
	return sympy.acot(x).evalf()
	
def acoth(x):
	return sympy.acoth(x).evalf()
	
def arccoth(x):
	return sympy.acoth(x).evalf()
	
def acsc(x):
	return sympy.acsc(x).evalf()
	
def arccsc(x):
	return sympy.acsc(x).evalf()
	
def acsch(x):
	return sympy.acsch(x).evalf()
	
def arccsch(x):
	return sympy.acsch(x).evalf()
	
def arcsec(x):
	return sympy.asec(x).evalf()
	
def asech(x):
	return sympy.asech(x).evalf()
	
def arcsech(x):
	return sympy.asech(x).evalf()
	
def asin(x):
	return sympy.asin(x).evalf()
	
def asinh(x):
	return sympy.asinh(x).evalf()
	
def arcsinh(x):
	return sympy.asinh(x).evalf()
	
def arcsin(x):
	return sympy.asin(x).evalf()	
	
def atan(x):
	return sympy.atan(x).evalf()
	
def arctan(x):
	return sympy.atan(x).evalf()
	
def atanh(x):
	return sympy.atanh(x).evalf()
	
def arctanh(x):
	return sympy.atanh(x).evalf()
	
def ceil(x):
	return sympy.ceiling(x).evalf()
	
def ceiling(x):
	return sympy.ceiling(x).evalf()
	
def cos(x):
	return sympy.cos(x).evalf()
	
def cosh(x):
	return sympy.cosh(x).evalf()
	
def cot(x):
	return sympy.cot(x).evalf()	
	
def coth(x):
	return sympy.coth(x).evalf()	
	
def csc(x):
	return sympy.csc(x).evalf()
	
def csch(x):
	return sympy.csch(x).evalf()	
	
def factorial(x):
	return sympy.factorial(x).evalf()
	
def exp(x):
	return sympy.exp(x).evalf()
	
def floor(x):
	return sympy.floor(x).evalf()
	
def ln(x):
	return sympy.ln(x).evalf()
	
def log(x):
	return sympy.log(x).evalf()
	
def log10(x):
	return sympy.log(x,10).evalf()
	
def piecewise(*x):
	if len(x) == 3:	
		return x[0] if x[1] else x[2]
	else:
		ans = []
		xlen = len(x) - 1
		for i in range(0,xlen,2):
			if x[i+1]:
				return x[i]
		return x[-1]
					
def pow(x,y):
	return sympy.Pow(x,y).evalf()
	
def power(x,y):
	return sympy.Pow(x,y).evalf()
	
def root(n,x):
	return sympy.root(x,n).evalf()
	
def sec(x):
	return sympy.sec(x).evalf()
	
def sech(x):
	return sympy.sech(x).evalf()
	
def sqr(x):
	return sympy.sqrt(x).evalf()
	
def sqrt(x):
	return sympy.sqrt(x).evalf()
	
def sin(x):
	return sympy.sin(x).evalf()
	
def sinh(x):
	return sympy.sinh(x).evalf()
	
def tan(x):
	return sympy.tan(x).evalf()
	
def tanh(x):
	return sympy.tanh(x).evalf()

def And(*x):
	for y in x:
		if y == False:
			return False
	return True

def Not(x):
	return not x
	
def Or(*x):
	for y in x:
		if y == True:
			return True
	return False
	
def xor(*x):
	odd = 0
	for y in x:
		if y == True:
			odd = odd + 1
	if odd % 2 == 0:
		return False
	else:
		return True
		
def eq(*x):
	last = None
	for y in x:
		if last == None:
			last = y
		else:
			if y != last:
				return False
	return True
	
def geq(*x):
	last = None
	for y in x:
		if last == None:
			last = y
		else:
			if last < y:
				return False
	return True
	
def gt(*x):
	last = None
	for y in x:
		if last == None:
			last = y
		else:
			if last <= y :
				return False
	return True
	
def leq(*x):
	last = None
	for y in x:
		if last == None:
			last = y
		else:
			if last > y:
				return False
	return True
	
def lt(*x):
	last = None
	for y in x:
		if last == None:
			last = y
		else:
			if last >= y :
				return False
	return True
	
def neq(x,y):
	if x == y:
		return False
	return True
	
def plus(*x):
	return sympy.fsum(x).evalf()
	
def times(*x):
	p = 1
	for y in x:
		p = p*y
	return p
	
def minus(x,y):
	return x - y
	
def divide(x,y):
	return x/y
	
def multiply(*x):
	return times(*x)
	
exponentiale = exp(1)