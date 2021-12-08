# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t, s=0, amp=0, freq=0):
    # constants
    k_0 = 0.01
    k_1 = 1
    k_2 = 1
    K_m1 = 0.05
    K_m2 = 0.05
    S = s
    amp = amp
    freq = freq
    S = S + amp * np.sin(freq * t)

    # assign each ODE to a vector element
    R_P = x[0]
    R_T = 1

    # define each ODE: sigmoidal model (c)
    dR_Pdt = (k_1 * S * (R_T - R_P)) / (K_m1 + R_T - R_P) - (k_2 * R_P) / (K_m2 + R_P)

    return [dR_Pdt]


# declare time vector
t = np.linspace(0, 50, 100)

S_values = np.linspace(0, 3, 100)
# S_values = [0, 0.5, 1, 1.5, 2]
# S_values = [1]
# initial condition
# R_P_values = [0.1]
# R_P_values = np.linspace(0, 1, 2)
R_P_values = np.linspace(0, 1, 10)
amps = [1]
# amps = np.linspace(0, 2, 5)
# freq = 0.1
# freqs = np.linspace(0, 1, 5)
freqs = [0]
colors = ['r', 'c', 'g', 'b', 'm']
counter2 = 0
steadystates = []
S_inits = []
for freq in freqs:
    for amp in amps:
        for R_P0 in R_P_values:
            counter = 0
            for S_init in S_values:
                x_0 = [R_P0]
                # x = odeint(odes, x_0, t, args=(S_init, amp, freq))
                x2 = odeint(odes, x_0, t, args=(S_init, 0, 0))
                # R_P = x[:, 0]
                R_P_const = x2[:, 0]
                # if counter < 5 and counter2 < 1:
                #     plt.plot(t, R_P, label=f'$S={S_init}$', color=colors[counter])
                # else:
                #     plt.plot(t, R_P, color=colors[counter])

                # if R_P[-1] > R_P0:
                #     color = 'lightskyblue'
                # else:
                #     color = 'pink'
                amp = 0
                freq = 0
                signal = S_init + amp * np.sin(freq * t)
                # plt.plot(signal, R_P)#, color=color)
                # plt.plot(t, R_P, label=f'$S_0 ={S_init} $', color=colors[counter])
                # plt.plot(t, R_P, label=f'$a ={amp} $')
                # plt.plot(t, R_P, label=f'$f ={freq} $')
                # plt.plot(t, R_P_const, '--', color=colors[counter])
                # plt.plot(S_init, R_P_const[-1], '.', color='black')
                S_inits.append(S_init)
                steadystates.append(R_P_const[-1])
                counter += 1
            counter2 += 1

plt.plot(S_inits[0:99], steadystates[0:99], '-', color='black')
# plt.axvline(x=1, ymin=0.1, ymax=0.9, color='b', linestyle='dashed')

# plt.title(f'$R_P(t)$ for oscillating $S = S_0 + {amp} * sin({freq}t)$ and $R_0={R_P_values[0]}$')
# plt.title(f'$R_P(t)$ for oscillating $S$ with varying frequency')
# plt.title('$R_P(t)$ for constant $S$ and different $R_{P0}$')
# plt.title('Steady state curve for sigmoidal model')
# plt.xlabel('$t$')
plt.xlabel('$S$')
plt.ylabel('$R_P$')
# plt.legend(loc='lower right')
# plt.legend()
plt.show()
