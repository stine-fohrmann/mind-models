# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function that returns dy/dt
def odes(x, t, s):
    # constants
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
    S = s
    S = S + 2*np.sin(t*20)

    # assign each ODE to a vector element
    R = x[0]
    X = x[1]

    def E(R):
        return goldbeter_koshland(k_3*R, k_4, J_3, J_4)

    # define each ODE: mutual inhibition model (f)
    dRdt = k_0*E(R) + k_1*S - k_2 * R -k_2prime * X * R
    dXdt = k_5*R - k_6*X

    return [dRdt, dXdt]


def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G


# declare time vector
t = np.linspace(0, 100, 1000)

Rt_plot = plt.figure(1)

# solve system of diff. eq.
S_values = [0.2]
# initial condition
R_0_values = [1, 1]

R_asymptote = []

for S in S_values:
    x = odeint(odes, R_0_values, t, args=(S,))
    R = x[:, 0]
    X = x[:, 1]
    R_asymptote.append(R[-1])
    plt.plot(X, R, label=f'R against X')
    #plt.plot(X, -R, label=f'-R against X')


#plt.title('$R(t)$ for different $S$ and $R_0=0$')
plt.xlabel('$X$')
plt.ylabel('$R$')
plt.legend(loc='lower right')
#plt.xlim([0, 1.5])
#plt.ylim([0, 3])


RS_plot = plt.figure(2)
plt.plot(t, X)
plt.plot(t, R)
plt.title('Stable values that $R$ approaches over time, depending on $S$')
plt.xlabel('$S$')
plt.ylabel(r'$R$ as $t \rightarrow \infty$')

plt.show()