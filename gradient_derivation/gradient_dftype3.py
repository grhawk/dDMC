import pylab as pl
from sympy import *
init_printing(use_unicode = True)
x,bij,s,a,R0 = symbols('x bij s a R0')

# s = 23
# a = 1
# b = 2
# r0 = 5


fd = 0.5*(1+tanh(s*(x/(a*R0)-1)))
tt = 1 - (exp(-bij*x) * (1 + bij*x + (bij*x)**2/2 + (bij*x)**3/6 + (bij*x)**4/24 + (bij*x)**5/120 + (bij*x)**6/720))

df = fd*tt*x**(-6)*-1
dfp = diff(df,x)
dfp2 = diff(df,x,2)

print(latex(factor(dfp)))
print factor(dfp)

# fdp = diff(fd,t)
# fdp2 = diff(fd,t,2)

# ttp = diff(tt,x)
# ttp2 = diff(tt,x,2)

# f = lambdify(t,fd,"numpy")
# fp = lambdify(t,fdp,"numpy")
# fp2 = lambdify(t,fdp2,"numpy")

# t = lambdify(x,tt,"numpy")
# tp = lambdify(x,ttp,"numpy")
# tp2 = lambdify(x,ttp2,"numpy")

r = pl.linspace(0,8,10000)

# pl.plot(r,f(r))
# pl.plot(r,fp(r))
# pl.plot(r,fp2(r))

# pl.plot(r,t(r))
# pl.plot(r,tp(r))
# pl.plot(r,tp2(r))


d = lambdify(x,df,"numpy")
dp = lambdify(x,dfp,"numpy")
dp2 = lambdify(x,dfp2,"numpy")

#pl.plot(r,d(r))
#pl.plot(r,dp(r))
#pl.plot(r,dp2(r))

#pl.show()
