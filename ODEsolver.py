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
    #S = s + a*np.cos(b*t)
    #S = s

    # activator inhibitor
    k_0 = 4
    k_1 = 1
    k_2 = 1
    #k_2 = 3
    k_2prime = 1
    k_3 = 1
    k_4 = 1
    k_5 = 0.1
    k_6 = 0.075
    j_3 = 0.3
    j_4 = 0.3

    # assign each ODE to a vector element


 
    # assign each ODE to a vector element
    R1 =  x[0]
    X   = x[1]
    R   = x[2]

    def E1(R1):
        return goldbeter_koshland(k_3*R1, k_4, j_3, j_4)
    
    def E(R):
        return goldbeter_koshland(K_3, K_4 * R, J_3, J_4)

    # negative feedback oscillator

    # activator inhibitor
    dR1dt = k_0*E1(R1) + k_1*S - k_2 * R1 -k_2prime * X * R1
    dXdt = k_5*R1 - k_6*X

    # mutual inhibition with negative feedback (R_P) as signal
    dRdt   = K_0 + K_1 * S - K_2 * R - K_2prime * E(R) * R
    #dRdt   = K_0 + K_1 * R1 - K_2 * R - K_2prime * E(R) * R

    return [dR1dt, dXdt, dRdt]

def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G


a_values = [2] # first value for amplibute to plot for
b_values = [1]
s_values = [1.25]
r_values  = [0.20748]
a_step = 0        # stepwise change of a
b_step = 0    # stepwise change of b (ish)
s_step = 0    # stepwise change of b (ish)
r_step = 0.00001    # stepwise change of b (ish)


plots = 2
for mul in range(plots):
    mul += 1           #accumulator
    # appending the different frequences and amplitudes
    a_values.append(a_values[0] + a_step * mul)
    b_values.append(b_values[0] + b_step * mul**1.2)
    s_values.append(s_values[0] + s_step * mul)
    r_values.append(r_values[0] + r_step * mul)



t = np.linspace(0, 1000,1300)

S = 1.3
# initial condition
X_0  = 5
R1_0 = 0
R_0  = 0.25
R_0_values  = 0.25
init_cond = [R1_0,X_0, R_0]
#signal = S + a_values[0]*np.sin(b_values[0]*t)
#print(f'max: [{max(signal)}')
#print(f'min: [{min(signal)}')

fig,ax = plt.subplots(plots+1,1,figsize = (11,5))
fig.suptitle("Response from mutual inhibition when \"negative feedback oscillator\" is the signal: for different freq and ampl")
colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(a_values))))

for i in range(len(a_values)):
    # iterate through the different a and b values to plot them
    init_cond[-1]=r_values[i]
    a = a_values[i]
    b = b_values[i]
    S = s_values[i]
    r = r_values[i]
    signal = S + a*np.sin(b*t)
    #ax[i].plot(t,signal, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "g")
    x = odeint(odes, init_cond, t, args=(S,0,b))
    R2   = x[:, 2]
    asymp = R2[-1]
    ax[i].plot(t,R2, color = "k")
    ax[i].plot(t[-1],R2[-1], "*", color = "b")
    x = odeint(odes, init_cond, t, args=(S,a,b))
    R1  = x[:, 0]
    X   = x[:, 1]
    R   = x[:, 2]
    #ax[i].plot(t,R1, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "m")
    #ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = next(colors))
    ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}, s = {round(S,3)}, r = {round(r,3)}, asymp = {round(R2[-1],3)}", color = next(colors))
    
    ax[i].set(ylabel="R")
    ax[i].set(xlabel="T")
    print(f's = {round(S,3)}, r = {round(r,3)}, asymp = {round(R2[-1],3)}')
    #print(f'max: {max(R1)}')
    #print(f'min: {min(R1)}')

fig.legend(loc='lower right')



#plt.show()