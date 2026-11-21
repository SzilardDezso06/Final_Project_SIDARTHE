import numpy as np
import matplotlib.pyplot as plt

# parameters 

alpha = 0.570
beta = 0.011
gamma = 0.456
delta = 0.011

epsilon = 0.171
theta = 0.371

zeta = 0.125
eta = 0.125

mu = 0.017
nu = 0.027

tau = 0.01

llambda = 0.034
rho = 0.034
kappa = 0.017
xi = 0.017
sigma = 0.017

# Defining R0 reproduction number

r1 = epsilon + zeta + llambda
r2 = eta + rho
r3 = theta + mu + kappa
r4 = nu + xi
r5 = sigma + tau

R0 = alpha/r1 + (beta*epsilon)/(r1*r2) + (gamma*zeta)/(r1*r3) + (delta*eta*epsilon)/(r1*r2*r4) + (delta*zeta*theta)/(r1*r3*r4)


# initial states of the population

S = 0.999
I = 0.001
D = 0.0
A = 0.0
R = 0.0
T = 0.0
H = 0.0
E = 0.0

y = [S, I, D, A, R, T, H, E]