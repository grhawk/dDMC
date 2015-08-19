import pylab as pl
from sympy import *
init_printing(use_unicode = True)
x,bij,s,a,R0 = symbols('x bij s a R0')

s = 23
a = 1
bij = 2
r0 = 5

# By definition fd > 0.99 at 1.1*a*Rvdw so the damping will be turned off when 
# the tt > 1-tol

fd = 0.5*(1+tanh(s*(x/(a*R0)-1)))
tt = 1 - (exp(-bij*x) * (1 + bij*x + (bij*x)**2/2 + (bij*x)**3/6 + (bij*x)**4/24 + (bij*x)**5/120 + (bij*x)**6/720))

solve(tt - 0.99 < 0, x)
