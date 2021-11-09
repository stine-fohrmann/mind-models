# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
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

    # assign each ODE to a vector element
    X = x[0]
    Y_P = x[1]
    R_P = x[2]

    # define each ODE: mutual inhibition model (f)
    dXdt = k_0 + k_1 * S - k_2 * X - k_2prime * R_P * X
    dY_pdt = (k_3 * X * (Y_T-Y_P))/(K_m3 + Y_T - Y_P) - (k_4 * Y_P)/(K_m4 + R_P)
    dR_Pdt = (k_5 * Y_P * (R_T - R_P))/(K_m5 + R_T - R_P) - (k_6 * R_P)/(K_m6 + R_P)

    return [dXdt, dY_pdt, dR_Pdt]

# declare time vector
t = np.linspace(0, 50, 100)

# S_values = np.linspace(0, 2, 10)
S_values = [2]

# initial condition
X_0 = 5
Y_P0 = 0.9
R_P0 = 0.1
init_cond = [X_0, Y_P0, R_P0]

# R_asymptote = []

for S in S_values:
    x = odeint(odes, init_cond, t, args=(S,))
    X = x[:, 0]
    Y_P = x[:, 1]
    R_P = x[:, 2]

fig, ax1 = plt.subplots()
plt.title('$X$, $Y_P$ and $R_P$ over time')

color = 'black'
ax1.set_xlabel('time')
ax1.set_ylabel('$X$', color=color)
plt.plot(t, X, label=f'$X$', color=color)
# plt.legend(loc='lower right')

ax2 = ax1.twinx()

color = 'r'
# ax2.set_xlabel('$time$')
ax2.set_ylabel('$Y_P$, $R_P$')#, color=color)
plt.plot(t, Y_P, label=f'$Y_P$', color=color)
# plt.legend(loc='lower right')
ax2.set_ylim([0, 1.1])

ax3 = ax1.twinx()

color = 'b'
# ax3.set_xlabel('time')
#ax3.set_ylabel('$R_P$', color=color)
plt.plot(t, R_P, label=f'$R_P$', color=color)
ax3.set_ylim([0, 1.1])

fig.tight_layout()
fig.legend(loc='lower right')
plt.show()
