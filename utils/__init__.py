from constants import amu, angstrom, eV, k_B, mass_kg
import numpy as np

def initialize_velocity_distribution(temperature):
    """
    Initialize velocity distribution for a system at a given temperature.
    
    Parameters:
    - temperature: Temperature in Kelvin
    
    """
    std_dev = np.sqrt(k_B * temperature / mass_kg)
    velocities = np.random.normal(0, std_dev, size=(3,))  # 3D velocity vector
    mean_velocity = np.mean(velocities)
    # Adjust velocities to have zero mean (Dc=3)
    velocities -= mean_velocity
    return velocities
