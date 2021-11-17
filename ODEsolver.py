# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t, s, a = 0.6, b = 0.1):


 

    # constants for mutual inhibition
    K_0 = 0
    K_1 = 0.05
    K_2 = 0.1
    K_2prime = 0.5
    K_3 = 0.2
    K_4 = 1
    J_3 = 0.05
    J_4 = 0.05
    #a = 0.6   # amplitude
    #b = 0.1   # frequency
    S = s + a*np.sin(b*t)
    S = s

    k_0 = 4
    k_1 = 1
    k_2 = 1
    k_2prime = 1
    k_3 = 1
    k_4 = 1
    k_5 = 0.1
    k_6 = 0.075
    J_3 = 0.3
    J_4 = 0.3

    # assign each ODE to a vector element


 
    # assign each ODE to a vector element
    R1 =  x[0]
    X   = x[1]
    R   = x[2]

    def E(R):
        return goldbeter_koshland(k_3*R1, k_4, J_3, J_4)

    # negative feedback oscillator

    # define each ODE: mutual inhibition model (f)
    dR1dt = k_0*E(R1) + k_1*S - k_2 * R1 -k_2prime * X * R1
    dXdt = k_5*R1 - k_6*X

    # mutual inhibition with negative feedback (R_P) as signal
    dRdt   = K_0 + K_1 * R - K_2 * R - K_2prime * E(R) * R

    return [dR1dt, dXdt, dRdt]

def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G


a_values = [1] # first value for amplibute to plot for
b_values = [0.07]
a_step = 0        # stepwise change of a
b_step = 0.005    # stepwise change of b (ish)

S = 3

plots = 1
for mul in range(plots):
    mul += 1           #accumulator
    # appending the different frequences and amplitudes
    a_values.append(a_values[0] + a_step * mul)
    b_values.append(b_values[0] + b_step * mul**1.2)


t = np.linspace(0, 500, 300)

# initial condition
X_0  = 5
Y_P0 = 0.9
R_P0 = 0.1
R_0  = 0.001
init_cond = [X_0, Y_P0, R_P0, R_0]

fig,ax = plt.subplots(plots+1,1,figsize = (11,5))
fig.suptitle("Response from mutual inhibition when \"negative feedback oscillator\" is the signal: for different freq and ampl")
fig.legend(loc='lower right')
colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(a_values))))

for i in range(len(a_values)):
    # iterate through the different a and b values to plot them
    a = a_values[i]
    b = b_values[i]
    x = odeint(odes, init_cond, t, args=(S,a,b))
    X   = x[:, 0]
    Y_P = x[:, 1]
    R_P = x[:, 2]
    R   = x[:, 3]
    #ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = next(colors))
    ax[i].plot(t,R_P, label = f"a = {round(a,3)}, b = {round(b,3)}", color = next(colors))
    ax[i].set(ylabel="R")
    ax[i].set(xlabel="T")



plt.show()