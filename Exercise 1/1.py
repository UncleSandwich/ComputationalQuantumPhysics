import numpy as np
from pylab import *
import cmath as cm

# Define constants
dx = 0.0001
a = 1
q = 1
m = 1
h = 1

# Define potential
b = int(a/dx)  # region that V_x is nonzero
n = 20*b
V_x = [0]*(2*n)

# Construct coordinates
x = []
for i in xrange(n):
    x.append(i*dx)
xp = []  # minus part of x axis
for i in xrange(n-1):
    xp.append(x[i]*-1)
xp.append(0)
xp.sort()
x = xp + x

# Construct potential
for i in xrange(n+1,n+1+b):
    # V_x[i] = 4*(((x[i]-a/2)/(a/2))**2-((x[i]-a/2)/(a/2))**4)
    # V_x[i] = 4*(x[i]/a-(x[i]*x[i])/(a*a))
    V_x[i]=10

# Define energy and the parameter of potential
E = 1
k = []
for i in xrange(len(V_x)):
    k.append(2*m*(E-V_x[i])/h**2)
# print k

# Initialize the wave function
fi = [0]*len(V_x)
C = 1
fi_0 = C * cm.exp(1j*q*a)

fi[n-1+b] = fi_0 # at a
fi[n-1+b+1] = fi_0*cm.exp(1j*q*dx)

# Iterate backward using Numerov method
for i in reversed(xrange(n+b)):
    fi_n_0=2*(1-5*dx*dx*k[i]/12)*fi[i]
    fi_n_1=-(1+dx*dx*k[i+1]/12)*fi[i+1]
    fi[i-1]=(fi_n_0+fi_n_1)/(1+dx*dx*k[i-1]/12)
 
# Calculate forward free wave
for i in xrange(n+b, len(V_x)-1):
    fi_n_0=2*(1-5*dx*dx*k[i]/12)*fi[i]
    fi_n__1=-(1+dx*dx*k[i-1]/12)*fi[i-1]
    fi[i+1]=(fi_n_0+fi_n__1)/(1+dx*dx*k[i+1]/12)

# Match the coefficients of free wave x<0
A = (fi[n-b]*cm.exp(1j*q*x[n-b])-fi[n-int(1.01*b)]*cm.exp(1j*q*x[n-int(1.01*b)]))/(cm.exp(2j*q*x[n-b])-cm.exp(2j*q*x[n-int(1.01*b)]))
B = (fi[n-b]*cm.exp(-1j*q*x[n-b])-fi[n-int(1.01*b)]*cm.exp(-1j*q*x[n-int(1.01*b)]))/(cm.exp(-2j*q*x[n-b])-cm.exp(-2j*q*x[n-int(1.01*b)]))

# Calculate the Transmit probability
T = 1/(abs(A))**2
# print T

# Compare T with exact solution
T_exact = 1/(1+((V_x[n+int(b/2)])**2*np.sinh(np.sqrt(2*m*(V_x[n+int(b/2)]-E)/h*a)))/(4*E*(V_x[n+int(b/2)]-E)))
print A, B, T, T_exact

# Visualize the wave and potential
fi = np.array(fi)

figure(figsize=(8,6), dpi=80)

# real part of wave
subplot(2,1,1)
plot(x,fi.real, color="blue", linewidth=1.0, linestyle="-")
plot(x,V_x, color="green", linewidth=1.0, linestyle="-")
xlim(-50, 50)
# ylim(-2,2)

# imaginary part of wave
subplot(2,1,2)
plot(x,fi.imag, color="blue", linewidth=1.0, linestyle="-")
plot(x,V_x, color="green", linewidth=1.0, linestyle="-")
xlim(-50, 50)
# ylim(-2,2)

show()

# Note for debugging:
# 1. increasing V_x and a will cause the free wave (x<0) to exponentially increase its amplitude. Don't know why.
# 2. T and T_exact do not match even in magnitude, but hasn't found any wrong during the calculation.
# 3. V_x or a, with T or T_exact, does not have well negative correlation.
