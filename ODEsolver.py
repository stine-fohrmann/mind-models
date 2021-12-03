# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t, s, a = 0.6, b = 0.1, freq = 1, mult = 1, add = 0):


 

    # constants for mutual inhibition
    K_0      = 0
    K_1      = 0.05
    K_2      = 0.1
    K_2prime = 0.5
    K_3      = 0.2
    K_4      = 1
    J_3      = 0.05
    J_4      = 0.05

    # activator inhibitor
    k_0      = 4
    k_1      = 1
    k_2      = 1
    k_2prime = 1
    k_3      = 1
    k_4      = 1
    k_5      = 0.1
    k_6      = 0.075
    j_3      = 0.3
    j_4      = 0.3
 
    S   = s + a*np.sin(b*t)
    # assign each ODE to a vector element
    R1 =  x[0]
    X   = x[1]
    R   = x[2]

    def E1(R1):
        return goldbeter_koshland(k_3*R1, k_4, j_3, j_4)
    
    def E(R):
        return goldbeter_koshland(K_3, K_4 * R, J_3, J_4)

    # activator inhibitor
    dR1dt = k_0*E1(R1) + k_1*S - k_2 * R1 -k_2prime * X * R1 
    dR1dt = dR1dt/freq # an attempt to manipulate the frequency of the oscillations
    dXdt  = k_5*R1     - k_6*X
    #dXdt  = dXdt *2

    # mutual inhibition with negative feedback (R_P) as signal
    #dRdt   = K_0 + K_1 * S - K_2 * R - K_2prime * E(R) * R
    dRdt   = K_0 + K_1 * (mult*R1+add) - K_2 * R - K_2prime * E(R) * R

    return [dR1dt, dXdt, dRdt]

def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G

''' creates some interesting results. for keeping
a_values = [0.2] # first value for amplibute to plot for
b_values = [0.2]

'''
""" For direct oscillatory signal for act-inhib
a_values = [0.6] # first value for amplibute to plot for
s_values = [1.3]    # for activ-inhib
b_values = [0.05]
"""

a_values = [0.21] # first value for amplibute to plot for
b_values = [0.18] # 0.18 works with 0.17 on a

a_values = [0.20] # first value for amplibute to plot for
b_values = [0.185] # 0.18 works with 0.17 on a

s_values = [0.2]    # for activ-inhib
r_values = [0.3]

# If we want to plot more things with different inputs
a_step = 0.01    # stepwise change of a
b_step = 0    # stepwise change of b 
s_step = 0    # stepwise change of signal
r_step = 0    # stepwise change of R_0 (ish)

plots = 1     # number of plots in the same window

try:
    for mul in range(plots-1):
        mul += 1           #accumulator
        # appending the different frequences and amplitudes
        a_values.append(a_values[0] + a_step * mul)
        b_values.append(b_values[0] + b_step * mul)
        s_values.append(s_values[0] + s_step * mul)
        r_values.append(r_values[0] + r_step * mul)
except:
    None

# Since actual max/min may be before the system has reached steady state, we have to take a later local extrema
def maxmin(lis):
    max =0; min = 0; counter = 1; maxVal,minVal = 0,0; peaks = []
    while max <  5 and min < 5:
        prev    = lis[counter-1]
        current = lis[counter]
        next    = lis[counter+1]
        if current > prev and current > next:
            maxVal = current
            max += 1
            peaks.append(counter)
        elif current < prev and current < next:
            minVal = current
            min += 1
        counter += 1
    return maxVal,minVal,peaks


t = np.linspace(0, 1000,500)

# initial condition
X_0  = 1
R1_0 = 1
init_cond = [R1_0,X_0, None] # R_0 is defined further up (r_values)

# creating the amount of plots given earlier

"""

"""
if plots != 1:
    fig,ax = plt.subplots(plots,1,figsize = (11,5))
else:
    fig,ax = plt.subplots(figsize = (11,5))
    ax = [ax]
""" TEMPORARY """
# fig,ax = plt.subplots(2,1,figsize = (11,5))

fig.suptitle("Response from mutual inhibition when \"negative feedback oscillator\" is the signal: for different freq and ampl")
colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(a_values))))

reqAmpl = 0.6 # The amplitude of oscillatory signal we know works 
reqCon  = 1.3 # the freq of the signal we know works
# reqAmpl = 0.7 # The amplitude of oscillatory signal we know works 
# reqCon  = 1.25 # the freq of the signal we know works
freq    = 7.4   # divide the R_P with some number, which will change frequencies, 8 works when there is no oscillatory signal in the activator inhibitor
for i in range(len(a_values)):
    # iterate through the different a and b values to plot them
    a = a_values[i]; b = b_values[i]; S = s_values[i]; r = r_values[i]
    init_cond[-1]=r
    signal = S + a*np.sin(b*t)
    # ax[i].plot(t,signal,  color = "g")
    
    ''' To extract information about R_P, so we can find out what values are needed to make is stay within the right values 
    E.g. we need to find out what to multiply the R_P signal with when used as input for dRdt'''
    x = odeint(odes, init_cond, t, args=(S,0,b, freq)) 
    R_P1  = x[:, 0]
    # ax[i].plot(t,R1, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "y") # plots the R_P, which is used as signal for dRdt
    maxR,minR,peaks = maxmin(R_P1)
    period = 1/(peaks[-1]-peaks[-2])
    print(period)
    amplR = (maxR-minR)/2                # The "amplitude" of the act-inhib
    midR  = maxR-amplR                   # The value which the act-inhib oscillates around
    mult  = reqAmpl/amplR                # what we want to multiply R_P with to stay within the right values
    add   = reqCon - midR*mult           # What we need to add to R_P to have the act-inhib at the right height
    R_P1  = [mult*x + add for x in R_P1]   # Making a list of the corrected act-inhib in order to plot it
    # ax[i].plot(t,R_P1, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "g") # plots the R_P, which is used as signal for dRdt
    
    x = odeint(odes, init_cond, t, args=(S,0,b,freq,mult,add )) # When dRdt get the right signal from R_P
    R1    = x[:, 2]
    ax[i].plot(t,R1, label = f"Mut-inhib + act-inhib without sinus", color = "b")
    # ax[i].plot(R_P1,R1, label = f"Mut-inhib + act-inhib without sinus", color = "m")

    x = odeint(odes, init_cond, t, args=(S,a,b,freq,mult,add )) # When dRdt get the right signal from R_P
    R_P = x[:, 0] # act-inhib including freq change and oscilatory signal
    X   = x[:, 1] 
    R   = x[:, 2] # mutual inhibition with act-inhib signal
    c = next(colors)
    R_P  = [mult*x + add for x in R_P]   # Making a list of the corrected act-inhib in order to plot it
    # ax[i].plot(R_P,R, label = f"Mut-inhib + act-inhib WITH sinus", color = c)
    ax[i].plot(t,  R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = c)
    # ax[i].plot(t,R_P,  "m") # plots the R_P, which is used as signal for dRdt
    

    # ax[i].plot(t,  R, label = f"Mut-inhib + act-inhib with sinus", color = "b")
    # for plotting in one figure
    # ax[1].plot(R_P,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = c)
    # ax[0].plot(t  ,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = c)
    
    # ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "g")
    # ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}, s = {round(S,3)}, r = {round(r,3)}, asymp = {round(R2[-1],3)}", color = next(colors))
    """
    print(f'max: {round(maxR,3)}')
    print(f'mid: {round(midR,3)}')
    print(f'min: {round(minR,3)}')
    print(f'amplitude: {round(amplR,3)}')
    """

    ax[i].set(ylabel="R")
    ax[i].set(xlabel="S")

fig.legend(loc='lower right')



plt.show()