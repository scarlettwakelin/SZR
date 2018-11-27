import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 1000
step = 10

# The SZR model differential equations.
def deriv(S, Z, R, alpha, beta, gamma):
    S, Z, R = y
    dSdt = -beta * S * Z
    dZdt = beta * S * Z + gamma *R - alpha * S * Z
    dRdt = alpha * S * Z - gamma * R
    return dSdt, dZdt, dRdt
    
def SZR(t, s, z, r, alpha, beta, gamma, delta):
    t[0] = 0
    s[0] = N
    z[0] = 0
    r[0] = 0
    dt = 1
    for i in range(1,step):
       print("i=%u" % i)
       t[i] = i
       s[i] = s[i-1] + dt* (-beta * s[i-1] * z[i-1])
       z[i] = z[i-1] + dt* (beta * s[i-1] * z[i-1] + delta *r[i-1] - alpha * s[i-1] * z[i-1])
       r[i] = r[i-1] + dt*(alpha * s[i-1] * z[i-1] - delta * r[i-1])

t = np.linspace(0, 160, 160)

t = [0]*step
s = [0]*step
z = [0]*step
r = [0]*step

#print(t)

SZR(t, s, z, r, 0.005, 0.0095, 0, 0.0001)

print(s) 

plt.clf() #clears plot
plt.xlabel('Time/ days')
plt.ylabel('Number')
title = 'SZR_model'
plt.title(title)
plt.plot(t, s, color='b', linestyle='-', label='Susceptible')
plt.plot(t, z, color='r', linestyle='-', label='Zombie')
plt.plot(t, r, color='g', linestyle='-', label='Removed')
plt.savefig(title+'.png')
plt.legend()
plt.show()  