import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 500
step = 200

# make sure x is between minimum and maximum
def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

def SIZR(t, s, i, z, r, alpha, beta, delta, zeta, rho, c):
    t[0] = 0
    s[0] = N
    i[0] = 0
    z[0] = 0
    r[0] = 0
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
       ZtoS = max(0, c * z[x-1])

       s[x] = s[x-1] + dt* (- StoI - StoR + ZtoS)
       i[x] = i[x-1] + dt* (StoI - ItoZ - ItoR)
       z[x] = z[x-1] + dt* (ItoZ + RtoZ - ZtoR - ZtoS)
       r[x] = r[x-1] + dt* (StoR + ItoR + ZtoR - RtoZ)

#       print("x=%u, total=%2.2f s=%f=>%f, i=%f=>%f, z=%f=>%f, r=%f=>%f, alpha * s[x-1] * z[x-1]=%2.2f" % (x, s[x]+i[x]+z[x]+r[x], s[x-1], s[x], i[x-1], i[x], z[x-1], z[x], r[x-1], r[x], alpha * s[x-1] * z[x-1]))
       if s[x] < 0 or i[x] < 0 or z[x] < 0 or r[x] < 0:
           print("BUG s=%f, i=%f, z=%f, r=%f" % (s[x], i[x], z[x], r[x]))

#t = np.linspace(0, 100, 100)

t = [0]*step
s = [0]*step
i = [0]*step
z = [0]*step
r = [0]*step

#SIZR(t, s, i, z, r, 0.005, 0.0095, 0.0001, 0.0001, 0.005)
SIZR(t, s, i, z, r, 0.00075, 0.0055, 0.0001, 0.09, 0.05, 0.001)

#print(t)

plt.clf() #clears plot
plt.xlabel('Time/ days')
plt.ylabel('Number')
title = 'SIZR_treatment_model'
plt.title(title)
plt.plot(t, s, color='b', linestyle='-', label='Susceptible')
plt.plot(t, i, color='c', linestyle='-', label='Infected')
plt.plot(t, z, color='r', linestyle='-', label='Zombie')
plt.plot(t, r, color='g', linestyle='-', label='Removed')
plt.savefig(title+'.png')
plt.legend()
plt.show()  