# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t):
    # constants
    k_0 = 0.01
    k_1 = 1
    k_2 = 5
    K_m1 = 0.05
    K_m2 = 0.05
    a = 1
    S = 2 + 0.1*np.sin(t)

    # assign each ODE to a vector element
    R = x[0]
    #R_P = x[2]
    #R_T = 1

    # define each ODE
    # sigmoidal:
    # dRdt = k_0 + k_1*S - k_2*R
    # dR_Pdt = (k_1*S * (R_T - R_P))/(K_m1 + R_T - R_P) - (k_2*R_P)/(K_m2+R_P)

    # linear:
    dRdt = k_0 + k_1 * S - k_2 * R
    #R_T = R + R_P
    #dR_Pdt = k_1 * S * (R_T - R_P) - k_2 * R_P

    return [dRdt]


# initial condition
x_0 = [1]

# print(odes(x=x_0, t=0))

# declare a time vector
t = np.linspace(0, 10, 100)

# solve system of diff. eq.
x = odeint(odes, x_0, t)

R = x[:, 0]
#R_P = x[:, 1]
#S = x[:, 2]

# plot results
plt.plot(t, R)
plt.xlabel('t')
plt.ylabel('$R$')
plt.show()
