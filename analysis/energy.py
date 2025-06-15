import h5py
import numpy as np
from utils.constants import amu, angstrom, eV, v_unit, a_unit, dt

def calculate_kinetic_energy(positions, velocity):
    """Calculate kinetic, potential and total energies"""
    kinetic_energy = 0.5 * np.sum(velocity**2)  # amu*(angstrom / dt)**2
    kinetic_energy *= amu * angstrom**2 / (dt**2 * eV)  # Convert to eV
    return kinetic_energy

def plot_energy(kinetic_energy, potential_energy, total_energy):
    """Plot kinetic, potential and total energies"""
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(kinetic_energy, label='Kinetic Energy (eV)')
    plt.plot(potential_energy, label='Potential Energy (eV)')
    plt.plot(total_energy, label='Total Energy', linestyle='--')
    plt.xlabel('Time Step')
    plt.ylabel('Energy (eV)')
    plt.title('Energy vs Time Step')
    plt.legend()
    plt.grid()
    plt.show()