import numpy as np
# from analysis.energy import calculate_kinetic_energy
from utils.constants import k_B

def calculate_heat_capacity(energy_list, temperature, **kwargs):
    """根据能量涨落计算热容"""
    nx, ny, nz = kwargs.get('nx', 6), kwargs.get('ny', 6), kwargs.get('nz', 4)
    energies = np.array(energy_list)
    mean_E = np.mean(energies)
    mean_E2 = np.mean(energies**2)
    
    variance = mean_E2 - mean_E**2
    Cv = variance / (k_B * temperature**2)
    # 对粒子数取平均———得到单位原子热容
    N = nx*ny*nz
    Cv/=N
    # 单位原子热容换算成摩尔原子热容
    Cv_molar = Cv * 1.602176634e-19 * 6.02214076e23
    return Cv_molar

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