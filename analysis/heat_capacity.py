import numpy as np
from analysis.energy import calculate_energies

def calculate_heat_capacity(energy_list, temperature):
    """Calculate heat capacity from energy fluctuations"""
    pass

def run_temperature_sweep(T_start=1, T_end=11, T_number=5):
    """Run simulations at different temperatures"""
    pass

def plot_heat_capacity(heat_capacity_data):
    """Plot heat capacity as a function of temperature"""
    import matplotlib.pyplot as plt
    temperatures = heat_capacity_data['temperatures']
    heat_capacities = heat_capacity_data['heat_capacities']
    
    plt.figure(figsize=(8, 6))
    plt.plot(temperatures, heat_capacities, marker='o')
    plt.title('Heat Capacity vs Temperature')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat Capacity (J/K)')
    plt.grid()
    plt.show()