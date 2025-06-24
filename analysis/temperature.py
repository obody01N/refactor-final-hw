import numpy as np
import matplotlib.pyplot as plt

def plot_kinetic_temperature(kinetic_energy_list, N_particles, temperature, steps, eV=1.602176634e-19 ,k_B=1.380649e-23):
    temperature_list = kinetic_energy_list * eV * 2 / ((3 * N_particles - 3) * k_B)  # K
    plt.figure(figsize=(8, 5))
    plt.plot(temperature_list, label="Kinetic Energy (eV)")
    # 采用对数刻度
    # plt.yscale('log')
    plt.xlabel("time (10fs)")
    plt.ylabel("temperature (K)")
    plt.title(f"temperature(K) calculated based on the molecular kinetic theory at {temperature} K")
    plt.legend()
    plt.grid()
    plt.savefig(f"result\\{temperature}K\\Kinetic temperature VS time (10fs).png")
    plt.show()
