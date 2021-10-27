# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t, s):
    # constants
    k_0 = 0
    k_1 = 0.05
    k_2 = 0.1
    k_2prime = 0.5
    k_3 = 1
    k_4 = 0.2
    J_3 = 0.05
    J_4 = 0.05
    S = s
    # S = S + 2*np.sin(t/10)

    # assign each ODE to a vector element
    R = x[0]

    def E(R):
        return goldbeter_koshland(k_3, k_4 * R, J_3, J_4)

    # define each ODE: mutual inhibition model (f)
    dRdt = k_0 + k_1 * S - k_2 * R - k_2prime * E(R) * R

    return [dRdt]


def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G


# declare time vector
t = np.linspace(0, 50, 100)

# solve system of diff. eq.
S_values = np.linspace(0, 10, 11)
# initial condition
R_0_values = [0]
for S in S_values:
    for R_0 in R_0_values:
        x_0 = [R_0]
        x = odeint(odes, x_0, t, args=(S,))
        R = x[:, 0]
        plt.plot(t, R, label=f'$S={S}$')

plt.title('$R(t)$ for different $S$ and $R_0=0$')
plt.xlabel('$t$')
plt.ylabel('$R(t)$')
plt.legend(loc='lower right')
plt.show()
