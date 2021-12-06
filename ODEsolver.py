# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
print("Running...")
def odesR(x, t, s=0, amp=0, freq=0):
    # constants
    k_0 = 0
    k_1 = 0.05
    k_2 = 0.1
    k_2prime = 0.5
    k_3 = 0.2
    k_4 = 1
    J_3 = 0.05
    J_4 = 0.05
    S = s
    a = amp  # 0.6   # amplitude
    b = freq  # 0.02   # frequency
    S = S + a * np.sin(b * t)

    # assign each ODE to a vector element
    R = x[0]

    def E(R):
        return goldbeter_koshland(k_3, k_4 * R, J_3, J_4)

    # define each ODE: mutual inhibition model (f)
    dRdt = k_0 + k_1 * S - k_2 * R - k_2prime * E(R) * R

    return [dRdt]

def odes(x, t, s, a = 0, b = 1, freq = 1, mult = 1, add = 0, sig = False):

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
    if not sig:
        dRdt   = K_0 + K_1 * (mult*R1+add) - K_2 * R - K_2prime * E(R) * R
        dRdt = dRdt*0.8
    else:
        dRdt   = K_0 + K_1 * S - K_2 * R - K_2prime * E(R) * R

    return [dR1dt, dXdt, dRdt]

def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G

# Since actual max/min may be before the system has reached steady state, we have to take a later local extrema
def maxmin(lis):
    length = len(lis)
    max =0; min = 0; counter = 1; maxVal,minVal = 0,0; peaks = []
    while (max <  5 and min < 5) and length > counter+1:
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

def findScales(response, reqAmpl, reqCon):
    maxR,minR,peaks = maxmin(response)
    period = 1/(peaks[-1]-peaks[-2])
    amplR = (maxR-minR)/2                # The "amplitude" of the act-inhib
    midR  = maxR-amplR                   # The value which the act-inhib oscillates around
    mult  = reqAmpl/amplR                # what we want to multiply R_P with to stay within the right values
    add   = reqCon - midR*mult           # What we need to add to R_P to have the act-inhib at the right height
    return mult, add

def scaleValues(list, mult, add):
    #list = [mult*x + add for x in list]   # Making a list of the corrected act-inhib in order to plot it
    return [mult*x + add for x in list]   # Making a list of the corrected act-inhib in order to plot it

def run_model(R_0, signal=(0, 0, 1), color=None):
    S_init = signal[0]
    amp = signal[1]
    freq = signal[2]
    x = odeint(odesR, R_0, t, args=(S_init, amp, freq))
    R = x[:, 0]
    signal = S_init + amp * np.sin(freq * t)
    if not color:
        if R[-1] > R_0:
            color = 'lightskyblue'
        else:
            color = 'pink'
    #ax[0].plot(signal, R, color=color)
    ax[0].plot(S_init, R[-1], '.', color='black')

t = np.linspace(0, 1000,700)


fig,ax = plt.subplots(2,1,figsize = (9,5))
fig.suptitle("Response from mutual inhibition when \"negative feedback oscillator\" is the signal: for different freq and ampl")



"""
R_0_values = np.concatenate([np.linspace(0, 0.1, 1), np.linspace(0.1, 0.3, 20), np.linspace(0.3, 1, 2)])
S_init_values = np.linspace(0, 2, 100)
for R0 in R_0_values:
    for S_init in S_init_values:
        constant_signal = (S_init, 0, 0)
        run_model(R0, constant_signal)

"""    
reqAmpl = 0.65 # The amplitude of oscillatory signal we know works 
reqCon  = 1.2 # the freq of the signal we know works
freq    = 7.4   # divide the R_P with some number, which will change frequencies, 8 works when there is no oscillatory signal in the activator inhibitor
freq    = 7.8  # divide the R_P with some number, which will change frequencies, 8 works when there is no oscillatory signal in the activator inhibitor
a = 0.22; 
a = 0.31; 
b = 0.216; 
b = 0.212; 

S = 0.2; 
r = 0.3
X_0  = 1
R1_0 = 1

# iterate through the different a and b values to plot them
init_cond = [R1_0,X_0, r] # R_0 is defined further up (r_values)
signal = S + a*np.sin(b*t)
# ax[i].plot(t,signal,  color = "g")

''' To extract information about R_P, so we can find out what values are needed to make is stay within the right values 
E.g. we need to find out what to multiply the R_P signal with when used as input for dRdt'''
x = odeint(odes, init_cond, t, args=(S,0,b, freq)) 
R_P1  = x[:, 0]
mult,add  = findScales(R_P1, reqAmpl, reqCon) # find scaling values
R_P1      = scaleValues(R_P1, mult, add)
# ax[1].plot(t,R_P1, label = f"a = {round(a,3)}, b = {round(b,3)}", color = "g") # plots the R_P, which is used as signal for dRdt

x = odeint(odes, init_cond, t, args=(S,0,b,freq,mult,add )) # When dRdt get the right signal from R_P
R1    = x[:, 2]

# ax[1].plot(t,R1, label = f"Mut-inhib + act-inhib without sinus", color = "y")
# ax[0].plot(R_P1,R1, label = f"Mut-inhib + act-inhib without sinus", color = "m")

x = odeint(odes, init_cond, t, args=(S,a,b,freq,mult,add )) # When dRdt get the right signal from R_P
R_P = x[:, 0] # act-inhib including freq change and oscilatory signal
X   = x[:, 1] 
R   = x[:, 2] # mutual inhibition with act-inhib signal
R_P  = scaleValues(R_P,mult,add)   # Making a list of the corrected act-inhib in order to plot it

def removeItems(nparray, amount):
    if type(nparray) != list:
        lis = nparray.tolist()
    else:
        lis = nparray
    for element in range(amount):
        lis.pop(0)
    return lis

"""Only plot after the system is stable
"""
R_P = removeItems(R_P,200)
R   = removeItems(R,200)
t   = removeItems(t,200)

ax[0].plot(R_P,R, label = f"Phaseplane: Mut-inhib + act-inhib WITH sinus", color = "g")
ax[1].plot(t,  R, label = f"Full mutual inhibition response against time", color = "m")
# ax[1].plot(t,R_P, label = f"Full activator inhibitor response against ", color = "b") # plots the R_P, which is used as signal for dRdt
# ax[1].plot(signal,  R, label = f"Full mutual inhibition response against time", color = "m")




ax[0].set(ylabel="R")
ax[0].set(xlabel="S")

fig.legend(loc='lower right')



plt.show()