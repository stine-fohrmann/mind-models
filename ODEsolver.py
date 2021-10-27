# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t, s):
    # constants
    k_0 = 0.01
    k_1 = 1
    k_2 = 1
    K_m1 = 0.05
    K_m2 = 0.05
    S = s
    S = S + 0.5*np.sin(t/10)

    # assign each ODE to a vector element
    R_P = x[0]
    R_T = 1

    # define each ODE: sigmoidal model (c)
    dR_Pdt = (k_1 * S * (R_T - R_P)) / (K_m1 + R_T - R_P) - (k_2 * R_P) / (K_m2 + R_P)

    return [dR_Pdt]


# declare time vector
t = np.linspace(0, 100, 100)

# solve system of diff. eq. for 0 < S < 4, R_0=1,2
# S_values = np.linspace(0, 2, 11)
S_values = [1]
# initial condition
R_P_values = [0.1]
for S in S_values:
    for R_P0 in R_P_values:
        x_0 = [R_P0]
        x = odeint(odes, x_0, t, args=(S,))
        R_P = x[:, 0]
        plt.plot(t, R_P, label=f'$S={S}$, ${R_P0=}$')

plt.title('$R_P(t)$ for oscillating $S$ and $R_{P0}=0.1$')
plt.xlabel('$t$')
plt.ylabel('$R_P(t)$')
plt.legend()
plt.show()
