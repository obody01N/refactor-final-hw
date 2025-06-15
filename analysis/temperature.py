import numpy as np
import matplotlib.pyplot as plt

# def compute_temperature(system, Nc=3, k_B=1.0):
#     Ek = 0.5 * system.mass * np.sum(system.v ** 2)
#     dof = 3 * system.N - Nc
#     T = 2 * Ek / (dof * k_B)
#     return T

def plot_kinetic_temperature(kinetic_energy_list, N_particles, temperature, eV=1.602176634e-19 ,k_B=1.380649e-23):
    temperature_list = kinetic_energy_list * eV * 2 / ((3 * N_particles - 3) * k_B)  # K
    times = np.arange(len(temperature_list))  # 10fs time step

    plt.figure(figsize=(8, 5))
    plt.plot(times[4:], temperature_list[4:], label="Kinetic Energy (eV)")
    # 采用对数刻度
    # plt.yscale('log')
    plt.xlabel("time (10fs)")
    plt.ylabel("temperature (K)")
    plt.title(f"temperature(K) calculated based on the molecular kinetic theory at {temperature} K")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig("Kinetic temperature VS time (10fs).png")

# print((mass_kg * 9e16) / 1.380649e-23 / (3*N_particles-3)*N_particles) # 43373348753114.875