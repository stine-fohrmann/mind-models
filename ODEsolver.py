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
    # S = 2 + 0.1*np.sin(t)

    # assign each ODE to a vector element
    R = x[0]
    R_P = x[1]
    R_T = 1

    # define each ODE: sigmoidal model (c)
    dRdt = k_0 + k_1 * S - k_2 * R
    dR_Pdt = (k_1 * S * (R_T - R_P)) / (K_m1 + R_T - R_P) - (k_2 * R_P) / (K_m2 + R_P)

    return [dRdt, dR_Pdt]


# initial condition
x_0 = [1, 1]

# declare time vector
t = np.linspace(0, 10, 100)

# solve system of diff. eq. for 0 < S < 5, R_0=1
S_values = range(5)
for S in S_values:
    x = odeint(odes, x_0, t, args=(S,))
    R = x[:, 0]
    R_P = x[:, 1]
    plt.plot(t, R)

plt.title('$R(t)$ for different values of $S$')
plt.xlabel('$t$')
plt.ylabel('$R(t)$')
plt.show()
