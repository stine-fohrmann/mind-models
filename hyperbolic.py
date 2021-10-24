# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t):
    # constants
    k_0 = 0.01
    k_1 = 1
    k_2 = 1
    R_T = 1
    a = 1
    S = 2 + 0.1*np.sin(t)
    #S = 10
    #S = 2*t
    # assign each ODE to a vector element
    R = x[0]
    R_P = x[1]

    # define each ODE
    # sigmoidal:

    dRdt = k_0 + k_1 * S - k_2 * R
    dR_Pdt = k_1 * S * (R_T - R_P) - k_2 * R_P

    return [dRdt, dR_Pdt]


# initial condition
x_0 = [1,1]

# print(odes(x=x_0, t=0))

# declare a time vector
t = np.linspace(0, 20, 500)

# solve system of diff. eq.
x = odeint(odes, x_0, t)

R = x[:, 0]
R_P = x[:, 1]
#S = x[:, 2]

# plot results
plt.plot(t, R_P)
plt.xlabel('t')
plt.ylabel('$R_P$')
plt.show()
