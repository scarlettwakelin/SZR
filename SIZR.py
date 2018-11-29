import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 500
step = 10

def SIZR(t, s, i, z, r, alpha, beta, gamma, zeta, rho):
    t[0] = 0
    s[0] = N
    i[0] = 0
    z[0] = 0
    r[0] = 0
    dt = 1
    for x in range(1,step):
       print("x=%u" % x)
       t[x] = x
       s[x] = s[x-1] + dt* (-beta * s[x-1] * z[x-1]-gamma * s[x-1])
       i[x] = i[x-1] + dt* (beta * s[x-1]* z[x-1]- rho* i[x-1]- gamma * i[x-1])
       z[x] = z[x-1] + dt* (rho* i[x-1] + zeta * r[x-1]  - alpha * s[x-1] * z[x-1])
       r[x] = r[x-1] + dt*(gamma * s[x-1] + gamma * i[x-1] + alpha * s[x-1] * z[x-1] - zeta * r[x-1])

#t = np.linspace(0, 100, 100)

t = [0]*step
s = [0]*step
i = [0]*step
z = [0]*step
r = [0]*step

#print(t)

SIZR(t, s, i, z, r, 0.005, 0.0095, 0.0001, 0.0001, 0.005)

print(s) 

plt.clf() #clears plot
plt.xlabel('Time/ days')
plt.ylabel('Number')
title = 'SIZR_model'
plt.title(title)
plt.plot(t, s, color='b', linestyle='-', label='Susceptible')
plt.plot(t, i, color='c', linestyle='-', label='Infected')
plt.plot(t, z, color='r', linestyle='-', label='Zombie')
plt.plot(t, r, color='g', linestyle='-', label='Removed')
plt.savefig(title+'.png')
plt.legend()
plt.show()  