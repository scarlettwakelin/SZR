import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 500
step = 200

# make sure x is between minimum and maximum
def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

def SIZRQ(t, s, i, z, r, q, alpha, beta, gamma, zeta, rho, kappa, sigma, delta):
    t[0] = 0
    s[0] = N
    i[0] = 0
    z[0] = 0
    r[0] = 0
    q[0] = 0
    dt = 1
    for x in range(1,step):
       t[x] = x

#       s[x] = s[x-1] + dt* (-beta * s[x-1] * z[x-1]-gamma * s[x-1])
#       i[x] = i[x-1] + dt* (beta * s[x-1]* z[x-1]- rho* i[x-1]- gamma * i[x-1])
#       z[x] = z[x-1] + dt* (rho* i[x-1] + zeta * r[x-1]  - alpha * s[x-1] * z[x-1])
#       r[x] = r[x-1] + dt*(gamma * s[x-1] + gamma * i[x-1] + alpha * s[x-1] * z[x-1] - zeta * r[x-1])

       StoI = clamp(0, beta * s[x-1] * z[x-1], s[x-1])
       StoR = max(0, delta * s[x-1])
       ItoR = max(0, delta * i[x-1])
       ItoZ = max(0, rho * i[x-1])
       RtoZ = max(0, zeta * r[x-1])
       ZtoR = clamp(0, alpha * s[x-1] * z[x-1], z[x-1])
       ItoQ = max(0, kappa * i[x-1])
       ZtoQ = max(0, sigma * z[x-1])
       QtoR = max(0, gamma * q[x-1] )

       s[x] = s[x-1] + dt* (- StoI - StoR)
       i[x] = i[x-1] + dt* (StoI - ItoZ - ItoR - ItoQ)
       z[x] = z[x-1] + dt* (ItoZ + RtoZ - ZtoR - ZtoQ)
       r[x] = r[x-1] + dt* (StoR + ItoR + ZtoR - RtoZ + QtoR)
       q[x] = q[x-1] + dt* (ItoQ + ZtoQ - QtoR)

#       print("x=%u, total=%2.2f s=%f=>%f, i=%f=>%f, z=%f=>%f, r=%f=>%f, alpha * s[x-1] * z[x-1]=%2.2f" % (x, s[x]+i[x]+z[x]+r[x], s[x-1], s[x], i[x-1], i[x], z[x-1], z[x], r[x-1], r[x], alpha * s[x-1] * z[x-1]))
       if s[x] < 0 or i[x] < 0 or z[x] < 0 or r[x] < 0 or q[x] < 0:
           print("BUG s=%f, i=%f, z=%f, r=%f, q=%f" % (s[x], i[x], z[x], r[x], q[x]))

#t = np.linspace(0, 100, 100)

t = [0]*step
s = [0]*step
i = [0]*step
z = [0]*step
r = [0]*step
q = [0]*step

#SIZR(t, s, i, z, r, 0.005, 0.0095, 0.0001, 0.0001, 0.005)
SIZRQ(t, s, i, z, r, q, 0.05, 0.0095, 0.0001, 0.0001, 0.05, 0.001, 0.001, 0.001)

#print(t)

plt.clf() #clears plot
plt.xlabel('Time/ days')
plt.ylabel('Number')
title = 'SIZRQ_model'
plt.title(title)
plt.plot(t, s, color='b', linestyle='-', label='Susceptible')
plt.plot(t, i, color='c', linestyle='-', label='Infected')
plt.plot(t, z, color='r', linestyle='-', label='Zombie')
plt.plot(t, r, color='g', linestyle='-', label='Removed')
plt.plot(t, q, color='m', linestyle='-', label='Quarantined')
plt.savefig(title+'.png')
plt.legend()
plt.show()  