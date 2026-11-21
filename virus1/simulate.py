import numpy as np
import matplotlib.pyplot as plt


# 1. Define the derivatives 
def derivatives(y, p):

    S, I, D, A, R, T, H, E = y

    # Infection flux
    new_infections = S * (p["alpha"]*I + p["beta"]*D + p["gamma"]*A + p["delta"]*R)

    dS = -new_infections

    dI = new_infections - (p["epsilon"] + p["zeta"] + p["llambda"]) * I

    dD = p["epsilon"]*I - (p["eta"] + p["rho"]) * D

    dA = p["zeta"]*I - (p["theta"] + p["mu"] + p["kappa"]) * A

    dR = p["eta"]*D + p["theta"]*A - (p["nu"] + p["xi"]) * R

    dT = p["mu"]*A + p["nu"]*R - (p["sigma"] + p["tau"]) * T

    dH = p["llambda"]*I + p["rho"]*D + p["kappa"]*A + p["xi"]*R + p["sigma"]*T

    dE = p["tau"]*T

    
    return np.array([dS, dI, dD, dA, dR, dT, dH, dE])

# Parameters 

params = {
    "alpha": 0.570,
    "beta": 0.011,
    "gamma": 0.456,
    "delta": 0.011,

    "epsilon": 0.171,
    "theta": 0.371,

    "zeta": 0.125,
    "eta": 0.125,

    "mu": 0.017,
    "nu": 0.027,

    "tau": 0.01,

    "llambda": 0.034,
    "rho": 0.034,
    "kappa": 0.017,
    "xi": 0.017,
    "sigma": 0.017,
}

# Defining R0 reproduction number
'''
r1 = epsilon + zeta + llambda
r2 = eta + rho
r3 = theta + mu + kappa
r4 = nu + xi
r5 = sigma + tau

R0 = alpha/r1 + (beta*epsilon)/(r1*r2) + (gamma*zeta)/(r1*r3) + (delta*eta*epsilon)/(r1*r2*r4) + (delta*zeta*theta)/(r1*r3*r4)
'''

# Initial states of the population

S = 0.999   # mostly healthy and susceptible population
I = 0.001   # small infected population
D = 0.0
A = 0.0
R = 0.0
T = 0.0
H = 0.0
E = 0.0

y0 = np.array([S, I, D, A, R, T, H, E])


# 2. Setup time grid and step size
t_start, t_end = 0.0, 300.0
dt = 0.01  

num_steps = int((t_end - t_start) / dt)

# Create arrays to store the history
t = np.linspace(t_start, t_end, num_steps + 1)

y = np.zeros((num_steps + 1, 8))

y[0] = y0    # store initial condition


# 3. The Solver Loop (RK4)
for i in range(num_steps):
    # k1: slope at the beginning of the interval
    k1 = derivatives(y[i], params)
    
    # k2: slope at the midpoint using k1
    k2 = derivatives(y[i] + 0.5*dt*k1, params)
    
    # k3: slope at the midpoint using k2
    k3 = derivatives(y[i] + 0.5*dt*k2, params)
    
    # k4: slope at the end using k3
    k4 = derivatives(y[i] + dt*k3, params)
    
    # Weighted average to update x and y
    y[i + 1] = y[i] + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)



# 4. Plot all 8 variables in one figure
labels = [
    "S (Susceptible)",
    "I (Infected)",
    "D (Diagnosed)",
    "A (Ailing)",
    "R (Recognized)",
    "T (Threatened)",
    "H (Healed)",
    "E (Extinct)",
]

plt.figure(figsize=(10, 6))
for idx, label in enumerate(labels):
    plt.plot(t, y[:, idx], label=label, linewidth=2)

plt.title("SIDARTHE Model Simulation (RK4)", fontsize=14)
plt.xlabel("Time (days)", fontsize=12)
plt.ylabel("Population Fraction", fontsize=12)
plt.xlim(t_start, t_end)
plt.ylim(0, 1.0)
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()