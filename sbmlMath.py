import sympy

#def abs(x):
	#return sympy.fabs(x)

def acos(x):
	return sympy.acos(x)
	
def arccos(x):
	return sympy.acos(x)
	
def acosh(x):
	return sympy.acosh(x)
	
def arccosh(x):
	return sympy.acosh(x)
	
def acot(x):
	return sympy.acot(x)
	
def arccot(x):
	return sympy.acot(x)
	
def acoth(x):
	return sympy.acoth(x)
	
def arccoth(x):
	return sympy.acoth(x)
	
def acsc(x):
	return sympy.acsc(x)
	
def acsch(x):
	return sympy.acsch(x)
	
def arccsch(x):
	return sympy.acsch(x)
	
def arcsec(x):
	return sympy.acsch(x)
	
def asech(x):
	return sympy.asech(x)
	
def arcsech(x):
	return sympy.asech(x)
	
def asin(x):
	return sympy.asin(x)
	
def arcsin(x):
	return sympy.asin(x)	
	
def atan(x):
	return sympy.atan(x)
	
def arctan(x):
	return sympy.atan(x)
	
def atanh(x):
	return sympy.atanh(x)
	
def arctanh(x):
	return sympy.atanh(x)
	
def ceil(x):
	return sympy.ceiling(x)
	
def ceiling(x):
	return sympy.ceiling(x)
	
def cos(x):
	return sympy.cos(x)
	
def cosh(x):
	return sympy.cosh(x)
	
def cot(x):
	return sympy.cot(x)	
	
def coth(x):
	return sympy.coth(x)	
	
def csc(x):
	return sympy.csc(x)
	
def csch(x):
	return sympy.csch(x)	
	
#def delay(x,y):
	
def factorial(x):
	return sympy.factorial(x)
	
def exp(x):
	return sympy.exp(x)
	
def floor(x):
	return sympy.floor(x)
	
def ln(x):
	return sympy.ln(x)
	
def log(x):
	return sympy.log(x,2)
	
def log10(x):
	return sympy.log(x,10)
	
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
	return sympy.Pow(x,y)
	
def power(x,y):
	return sympy.Pow(x,y)
	
def root(n,x):
	return sympy.root(x,n)
	
def sec(x):
	return sympy.sec(x)
	
def sech(x):
	return sympy.sech(x)
	
def sqr(x):
	return sympy.sqrt(x)
	
def sqrt(x):
	return sympy.sqrt(x)
	
def sin(x):
	return sympy.sin(x)
	
def sinh(x):
	return sympy.sinh(x)
	
def tan(x):
	return sympy.tan(x)
	
def tanh(x):
	return sympy.tanh(x)

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
	return sympy.fsum(x)
	
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
	
