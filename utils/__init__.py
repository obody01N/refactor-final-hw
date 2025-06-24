# import sys
# print(sys.path)
import numpy as np

def initialize_velocity_distribution(temperature, N_particles):
    """
    Initialize velocity distribution for a system at a given temperature.
    
    Parameters:
    - temperature: Temperature in Kelvin
    
    """
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    mass=12
    mass_kg = mass * 1.66053906660e-27  # Convert mass from amu to kg
    std_dev = np.sqrt(k_B * temperature / mass_kg)
    velocities = np.random.normal(0, std_dev, (N_particles, 3))  # 3D velocity vector
    mean_velocity = np.mean(velocities)
    # Adjust velocities to have zero mean (Dc=3)
    velocities -= mean_velocity
    return velocities
