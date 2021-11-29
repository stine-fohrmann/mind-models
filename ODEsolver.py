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

    # plt.figure(1)
    # plt.plot(t, S, '.', color='r')

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
t = np.linspace(0, 1000, 1000)

# Rt_plot = plt.figure(1)

# solve system of diff. eq.
S_values = [1.3]
# initial condition
R_0_values = [1]

# create figure with subplots

constant_signal = (1.3, 0, 0)


def run_model(R_0, signal=(0, 0, 0)):
    S_init = signal[0]
    amp = signal[1]
    freq = signal[2]
    x = odeint(odes, R_0, t, args=(S_init, amp, freq))
    R = x[:, 0]
    signal = S_init + amp * np.sin(freq * t)
    # ax1.plot(t, signal)
    ax1.plot(signal, R)
    ax2.plot(t, R, label=f'$S={round(S_init, 2)}$')


fig, (ax1, ax2) = plt.subplots(2)
S_init = 1.3
constant_signal = (S_init, 0, 0)
run_model(1, constant_signal)
amp = 0.6
freq = 0.02
osc_signal = (S_init, amp, freq)
run_model(1, osc_signal)
fig.suptitle(f'R(t) for S={constant_signal[0]} and S={S_init}+{amp}sin({freq}t)')

fig, (ax1, ax2) = plt.subplots(2)
S_init = 1
constant_signal = (S_init, 0, 0)
run_model(1, constant_signal)
amp = 0.72
freq = 0.02
run_model(1, (S_init, amp, freq))
fig.suptitle(f'R(t) for S={constant_signal[0]} and S={S_init}+{amp}sin({freq}t)')

fig, (ax1, ax2) = plt.subplots(2)
S_init = 1
constant_signal = (S_init, 0, 0)
run_model(1, constant_signal)
amp = 0.6
freq = 0.02
run_model(1, (S_init, amp, freq))
fig.suptitle(f'R(t) for S={constant_signal[0]} and S={S_init}+{amp}sin({freq}t)')

# plt.title(f'$R(t)$ for $S = {S_values[0]} + 0.6 sin(0.02 t)$ and $R_0={R_0_values[0]}$')
# plt.ylim(0, 2)
# plt.xlabel('$t$')
# plt.ylabel('$R(t)$')
# plt.legend(loc='lower right')

# RS_plot = plt.figure(2)
# plt.plot(S_values, R_asymptote, '.')
# plt.title('Stable values that $R$ approaches over time, depending on $S$')
# plt.xlabel('$S$')
# plt.ylabel(r'$R$ as $t \rightarrow \infty$')


plt.show()
