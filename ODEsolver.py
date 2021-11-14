# dependencies
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def odes(x, t, s, a = 0.6, b = 0.1):


    # constants
    k_0 = 0
    k_1 = 1
    k_1 = .9
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
    #a = 0.6   # amplitude
    #b = 0.1   # frequency
    S = S + a*np.sin(b*t)


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
    dXdt   = k_0 + k_1 * S - k_2 * X - k_2prime * R_P * X
    dY_pdt = (k_3 * X * (Y_T-Y_P))/(K_m3 + Y_T - Y_P) - (k_4 * Y_P)/(K_m4 + R_P)
    dR_Pdt = (k_5 * Y_P * (R_T - R_P))/(K_m5 + R_T - R_P) - (k_6 * R_P)/(K_m6 + R_P)
    dRdt   = K_0 + K_1 * R_P - K_2 * R - K_2prime * E(R) * R

    return [dXdt, dY_pdt, dR_Pdt, dRdt]


def goldbeter_koshland(u, v, J, K):
    G = (2 * u * K) / (v - u + v * J + u * K + np.sqrt((v - u + v * J + u * K) ** 2 - 4 * (v - u) * u * K))
    return G

t = np.linspace(0, 350, 300)
S_values = [1]
#a_values = [0.3]
#b_values = [0.0953]
a_values = [.8]
b_values = [0.065]
#amul = 1.6
#bmul = 0.8
amul = 0.2
bmul = 0.8
ainit,binit = a_values[0],b_values[0]
am, bm = a_values[0],b_values[0]
for mul in range(2):
    mul = mul + 1
    am = amul*mul + ainit
    bm = bmul* mul + binit
    print(f"a={round(am,4)}, b={round(bm,4)}")
    a_values.append(a_values[0]*am)
    b_values.append(b_values[0]*bm)


# initial condition
X_0 = 5
Y_P0 = 0.9
R_P0 = 0.1
R    = 0.01
init_cond = [X_0, Y_P0, R_P0, R]

#plt.figure(figsize=(14, 6), dpi=80)
fig,ax = plt.subplots(3,1,figsize = (11,5))
#fig.figure(figsize=(14, 6), dpi=80)

# R_asymptote = []
colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(a_values))))
S = 3
for i in range(len(a_values)):
#for a,b in a_values,b_values:
    a = a_values[i]
    b = b_values[i]
    x = odeint(odes, init_cond, t, args=(S,a,b))
    X   = x[:, 0]
    Y_P = x[:, 1]
    R_P = x[:, 2]
    R   = x[:, 3]
    c = next(colors)
    print(c)
    ax[i].plot(t,R, label = f"a = {round(a,3)}, b = {round(b,3)}", color = c)
#plt.plot(t,R,   color = 'm', label = "R")


# declare time vector


#RS_plot = plt.figure(1)

#plt.plot(t,R_P, color = 'r', label = "$R_P$")
#plt.plot(t,X,   color = 'b', label = "X")
#plt.plot(t,Y_P, color = 'g', label = "$Y_P$")
#plt.title('Stable values that $R$ approaches over time, depending on $S$')
#plt.xlabel('$S$')
#plt.ylabel(r'$R$ as $t \rightarrow \infty$')
#plt.legend(loc='lower right')
fig.legend(loc='lower right')
plt.show()