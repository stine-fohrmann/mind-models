# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t, s=0, amp=0, freq=0):
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


def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G


# declare time vector
t = np.linspace(0, 1000, 2000)

steadystates1 = []
steadystates2 = []
S_inits1 = []
S_inits2 = []


def run_model(R_0, signal=(0, 0, 0), color=None):
    S_init = signal[0]
    amp = signal[1]
    freq = signal[2]
    x = odeint(odes, R_0, t, args=(S_init, amp, freq))
    R = x[:, 0]
    signal = S_init + amp * np.sin(freq * t)
    # ax1.plot(t, signal)
    if not color:
        if R[-1] > R_0:
            color = 'lightskyblue'
        else:
            color = 'pink'
    # ax1.plot(signal, R, color=color)
    # ax1.plot(S_init, R[-1], '.', color='black')
    # ax2.plot(t, R, label=f'$S={round(S_init, 2)}$')
    # print(f'{S_init=}')
    # print(R[-1])

    if S_init < 1.8 and R[-1] < 0.2:
        S_inits1.append(S_init)
        steadystates1.append(R[-1])
    elif S_init > 0.7 and R[-1] > 0.3:
        S_inits2.append(S_init)
        steadystates2.append(R[-1])


fig, (ax1) = plt.subplots(1)

R_0_values = np.concatenate([np.linspace(0, 0.1, 1), np.linspace(0.1, 0.3, 20), np.linspace(0.3, 1, 2)])
S_init_values = np.linspace(0, 2, 100)
for R_0 in R_0_values:
    for S_init in S_init_values:
        constant_signal = (S_init, 0, 0)
        run_model(R_0, constant_signal)

ax1.plot(S_inits1[0:85], steadystates1[0:85], '-', color='black')
ax1.plot(S_inits1[84], steadystates1[84], '.', color='r', label='$S_{crit, 1}$')
i = 743+14
j = 783+31
ax1.plot(S_inits2[i:j], steadystates2[i:j], '-', color='black')
ax1.plot(S_inits2[i], steadystates2[i], '.', color='b', label='$S_{crit, 2}$')

for k in range(j-i):
    print(f'{S_inits2[k+i]=}, {steadystates2[i+k]}')

print(len(S_inits2))

R_0 = 0
S_init = 1.3
amp = 0.6
freq = 0.02
constant_signal = (S_init, 0, 0)
# run_model(R_0, constant_signal)
# run_model(R_0, (S_init, amp, freq), color='green')
# fig.suptitle(f'R(t) for S={constant_signal[0]} and S={S_init}+{amp}sin({freq}t)')
# fig.suptitle(f'S,R-diagram for S={S_init}+{amp}sin({freq}t)')
ax1.set_xlabel('$S$')
ax1.set_ylabel('$R$')
ax1.legend()

plt.show()
