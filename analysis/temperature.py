import numpy as np

def compute_temperature(system, Nc=3, k_B=1.0):
    Ek = 0.5 * system.mass * np.sum(system.v ** 2)
    dof = 3 * system.N - Nc
    T = 2 * Ek / (dof * k_B)
    return T
