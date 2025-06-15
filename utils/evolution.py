import numpy as np
from utils.constants import amu, angstrom, eV, v_unit, a_unit

def velocity_verlet(positions, velocity, acceleration, forces, mass=12.0):
    """Velocity Verlet integration"""
    new_acceleration = forces / mass * a_unit # a(t+dt)
    positions += velocity + 0.5 * acceleration # r(t+dt)
    velocity += 0.5 * (acceleration + new_acceleration) # v(t+dt/2)
    return positions, velocity, new_acceleration
