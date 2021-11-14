# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes2(x, t, s):
    # constants
    k_0 = 0
    k_1 = 1
    k_2 = 0.01
    k_2prime = 10
    k_3 = 0.2
    k_4 = 0.2
    k_5 = 0.1
    k_6 = 0.04
    Y_T = 1
    R_T = 1
    K_m3 = 0.01
    K_m4 = 0.01
    K_m5 = 0.01
    K_m6 = 0.01
    S = s
    # S = S + 2*np.sin(t/10)

    # assign each ODE to a vector element
    X = x[0]
    Y_P = x[1]
    R_P = x[2]
    # define each ODE: mutual inhibition model (f)
    dXdt = k_0 + k_1 * S - k_2 * X - k_2prime * R_P * X
    dY_pdt = (k_3 * X * (Y_T-Y_P))/(K_m3 + Y_T - Y_P) - (k_4 * Y_P)/(K_m4 + R_P)
    dR_Pdt = (k_5 * Y_P * (R_T - R_P))/(K_m5 + R_T - R_P) - (k_6 * R_P)/(K_m6 + R_P)

    return [dXdt, dY_pdt, dR_Pdt]

def odes(x, t, s):


    # constants
    k_0 = 0
    k_1 = 1
    k_2 = 0.01
    k_2prime = 10
    k_3 = 0.2
    k_4 = 0.2
    k_5 = 0.1
    k_6 = 0.04
    Y_T = 1
    R_T = 1
    K_m3 = 0.01
    K_m4 = 0.01
    K_m5 = 0.01
    K_m6 = 0.01
    S = s

    # S = S + 2*np.sin(t/10)
    K_0 = 0
    K_1 = 0.05
    K_2 = 0.1
    K_2prime = 0.5
    K_3 = 0.2
    K_4 = 1
    J_3 = 0.05
    J_4 = 0.05
    S = s
    a = 0.6   # amplitude
    b = 0.1   # frequency
    #S = S + a*np.sin(b*t)


    # assign each ODE to a vector element

    def E(R):
        return goldbeter_koshland(K_3, K_4 * R, J_3, J_4)

    # define each ODE: mutual inhibition model (f)

    # assign each ODE to a vector element
    X = x[0]
    Y_P = x[1]
    R_P = x[2]
    R = x[3]


    # define each ODE: mutual inhibition model (f)
    dXdt = k_0 + k_1 * S - k_2 * X - k_2prime * R_P * X
    dY_pdt = (k_3 * X * (Y_T-Y_P))/(K_m3 + Y_T - Y_P) - (k_4 * Y_P)/(K_m4 + R_P)
    dR_Pdt = (k_5 * Y_P * (R_T - R_P))/(K_m5 + R_T - R_P) - (k_6 * R_P)/(K_m6 + R_P)
    dRdt = K_0 + K_1 * R_P - K_2 * R - K_2prime * E(R) * R

    return [dXdt, dY_pdt, dR_Pdt, dRdt]



# function that returns dy/dt
def odesT(x, t, s = 1):
    # constants
    #r = odeint(odes2, x, t, args=(S,))
    #S = r[:, 0]
    #S = odes2(x,t,1)

    # ------- NEGATIVE FEEDBACK ---------
    k_0 = 0
    k_1 = 1
    k_2 = 0.01
    k_2prime = 10
    k_3 = 0.2
    k_4 = 0.2
    k_5 = 0.1
    k_6 = 0.04
    Y_T = 1
    R_T = 1
    K_m3 = 0.01
    K_m4 = 0.01
    K_m5 = 0.01
    K_m6 = 0.01
    S = 1
    # S = S + 2*np.sin(t/10)

    # assign each ODE to a vector element
    # define each ODE: mutual inhibition model (f)

    # ---- NEGATIVE FEEDBACK END ------------
    K_0 = 0
    K_1 = 0.05
    K_2 = 0.1
    K_2prime = 0.5
    K_3 = 0.2
    K_4 = 1
    J_3 = 0.05
    J_4 = 0.05
    #S = 1
    # S = S + 2*np.sin(t/10)
    # assign each ODE to a vector element
    X = x[0]
    Y_P = x[1]
    R_P = x[2]
    R = x[3]

    def E(R):
        return goldbeter_koshland(K_3, K_4 * R, J_3, J_4)

    # define each ODE: mutual inhibition model (f)
    #dRdt = K_0 + K_1 * S - K_2 * R - K_2prime * E(R) * R
    dXdt   = k_0 + k_1 * S - k_2 * X - k_2prime * R_P * X
    dY_pdt = (k_3 * X * (Y_T-Y_P))/(K_m3 + Y_T - Y_P) - (k_4 * Y_P)/(K_m4 + R_P)
    dR_Pdt = (k_5 * Y_P * (R_T - R_P))/(K_m5 + R_T - R_P) - (k_6 * R_P)/(K_m6 + R_P)
    dRdt   = K_0 + K_1 * R_P - K_2 * R - K_2prime * E(R) * R

    return [dXdt, dY_pdt, dR_Pdt,dRdt]


def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G

t = np.linspace(0, 100, 100)
S_values = [2]

# initial condition
X_0 = 5
Y_P0 = 0.9
R_P0 = 0.1
init_cond = [X_0, Y_P0, R_P0,0]

# R_asymptote = []

for S in S_values:
    x = odeint(odes, init_cond, t, args=(S,))
    X   = x[:, 0]
    Y_P = x[:, 1]
    R_P = x[:, 2]
    R   = x[:, 3]
# declare time vector

"""
#Rt_plot = plt.figure(1)

# solve system of diff. eq.
#S_values = np.linspace(0, 2, 10)
#S_values = np.linspace(0, 2, 10)
S = 1
# initial condition
R_0_values = [1]
x_0 = [1,1,1,1]
#signal = odeint(odes2, x_0, t, args=(S,))
#x = odeint(odes, x_0, t, args=(signal,))
x = odeint(odes, x_0, t)
#x = odeint(odes2, x_0, t,args=(1,)
# )"""
#R = x[:, 3]
#R_P = x[:,2]

#R_asymptote = []
"""
for S in S_values:
    for R_0 in R_0_values:
        x_0 = [R_0]
        x = odeint(odes, x_0, t, args=(S,))
        R = x[:, 0]
        R_asymptote.append(R[-1])
        plt.plot(t, R, label=f'$S={round(S,2)}$')
"""

RS_plot = plt.figure(1)
plt.plot(t,R,   color = 'm', label = "R")
plt.plot(t,R_P, color = 'r', label = "$R_P$")
plt.plot(t,X,   color = 'b', label = "X")
plt.plot(t,Y_P, color = 'g', label = "$Y_P$")
plt.title('Stable values that $R$ approaches over time, depending on $S$')
plt.xlabel('$S$')
plt.ylabel(r'$R$ as $t \rightarrow \infty$')
plt.legend(loc='lower right')
plt.show()