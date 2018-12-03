import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 500
step = 20

# make sure x is between minimum and maximum
def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

def SZR(t, s, z, r, alpha, beta, delta, zeta):
    t[0] = 0
    s[0] = N
    z[0] = 0
    r[0] = 0
    dt = 1
    for x in range(1,step):
       t[x] = x

#       s[x] = s[x-1] + dt* (-beta * s[x-1] * z[x-1]-gamma * s[x-1])
#       i[x] = i[x-1] + dt* (beta * s[x-1]* z[x-1]- rho* i[x-1]- gamma * i[x-1])
#       z[x] = z[x-1] + dt* (rho* i[x-1] + zeta * r[x-1]  - alpha * s[x-1] * z[x-1])
#       r[x] = r[x-1] + dt*(gamma * s[x-1] + gamma * i[x-1] + alpha * s[x-1] * z[x-1] - zeta * r[x-1])

       StoZ = clamp(0, beta * s[x-1] * z[x-1], s[x-1])
       StoR = max(0, delta * s[x-1])
       RtoZ = max(0, zeta * r[x-1])
       ZtoR = clamp(0, alpha * s[x-1] * z[x-1], z[x-1])

       s[x] = s[x-1] + dt* (- StoZ - StoR)
       z[x] = z[x-1] + dt* (StoZ + RtoZ - ZtoR)
       r[x] = r[x-1] + dt* (StoR + ZtoR - RtoZ)

#       print("x=%u, total=%2.2f s=%f=>%f, i=%f=>%f, z=%f=>%f, r=%f=>%f, alpha * s[x-1] * z[x-1]=%2.2f" % (x, s[x]+i[x]+z[x]+r[x], s[x-1], s[x], i[x-1], i[x], z[x-1], z[x], r[x-1], r[x], alpha * s[x-1] * z[x-1]))
       if s[x] < 0 or z[x] < 0 or r[x] < 0:
           print("BUG s=%f, z=%f, r=%f" % (s[x], z[x], r[x]))

#t = np.linspace(0, 100, 100)

t = [0]*step
s = [0]*step
z = [0]*step
r = [0]*step

#SIZR(t, s, i, z, r, 0.005, 0.0095, 0.0001, 0.0001, 0.005)
SZR(t, s, z, r, 0.005, 0.0095, 0.0001, 0.0001)

#print(t)

plt.clf() #clears plot
plt.xlabel('Time/ days')
plt.ylabel('Number')
title = 'SZR_model'
plt.title(title)
plt.plot(t, s, color='b', linestyle='-', label='Susceptible')
plt.plot(t, z, color='r', linestyle='-', label='Zombie')
#plt.plot(t, r, color='g', linestyle='-', label='Removed')
plt.savefig(title+'.png')
plt.legend()
plt.show()  