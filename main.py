import numpy as np
from utils.config import Config
from utils.lattice import generate_graphite, calculate_box_size, calculate_supercell_vectors
from utils.constants import v_unit, a_unit, amu, angstrom, eV
from utils.__init__ import initialize_velocity_distribution
from utils.evolution import velocity_verlet
from utils.boundary import apply_boundary_conditions, apply_boundary_conditions_2
from utils.potentials import compute_forces
from analysis.visualization import plot_initial_structure, animate_trajectory
from analysis.heat_capacity import calculate_heat_capacity, plot_heat_capacity
from analysis.energy import calculate_kinetic_energy, plot_energy
from analysis.msd import compute_diffusion_coefficient, compute_diffusion_coefficient_2, plot_msd, plot_diffusion_coefficient
from analysis.temperature import plot_kinetic_temperature
from analysis.rdf import plot_rdf
# import h5py
import os

def main():
    """Run a molecular dynamics simulation of graphite.

    positions - angstrom
    velocities - m/s = v_unit angstrom / dt
    dt - seconds
    forces - eV/angstrom
    acceleration - eV/(angstrom * amu) = a_unit angstrom / dt^2
    v_unit = 1/(angstrom / dt) = 1e-4
    a_unit = eV/(angstrom*amu) / (angstrom / dt^2) = 0.9648533215665326
    """
    config = Config()
    # Initialize system using config parameters
    # positions: N x 3 array of atom positions in Angstrom
    N_particles = config.nx * config.ny * config.nz * len(config.basis_atoms)
    animation_pos = np.zeros((config.steps, N_particles, 3))
    animation_vel = np.zeros((config.steps, N_particles, 3))
    animation_acc = np.zeros((config.steps, N_particles, 3))
    positions, boundary = generate_graphite(
        config.nx, config.ny, config.nz,
        basis_vectors=config.basis_vectors,
        basis_atoms=config.basis_atoms
    )
    box_size = calculate_box_size(config.nx, config.ny, config.nz)
    positions = apply_boundary_conditions(positions, box_size)  # Apply periodic boundary conditions
    animation_pos[0] = positions.copy()
    initial_positions = positions.copy()
    msd_list = np.zeros(config.steps)
    kinetic_energy = np.zeros(config.steps)
    potential_energy = np.zeros(config.steps)
    heat_capacity_list = []
    plot_initial_structure(positions)
    
    # simulation at different temperatures
    for temperature in config.temperature:
        print(f"Running simulation at temperature: {temperature} K")
        os.makedirs(f"result\\{temperature}K", exist_ok=True)
        # t = 0:
        initial_velocity = initialize_velocity_distribution(temperature, N_particles) # v(t=0)
        initial_velocity *= v_unit  # Convert to angstrom / dt
        force, acceleration, current_potential_energy = compute_forces(positions, box_size, two_potentials=config.two_potentials)
        msd_list[0] = 0
        animation_vel[0] = initial_velocity.copy()
        animation_acc[0] = acceleration.copy()
        potential_energy[0] = current_potential_energy
        kinetic_energy[0] = calculate_kinetic_energy(positions, initial_velocity)
        # t = 1:
        positions = positions + initial_velocity + 0.5 * acceleration
        wrapped_positions = apply_boundary_conditions(positions, box_size)
        force, acceleration, current_potential_energy = compute_forces(wrapped_positions, box_size, two_potentials=config.two_potentials)
        velocity = initial_velocity + 0.5 * acceleration  # v(t=0.5)
        animation_pos[1] = wrapped_positions.copy()
        animation_vel[1] = velocity.copy()
        animation_acc[1] = acceleration.copy()
        msd_list[1] = np.mean(np.sum((wrapped_positions - initial_positions) ** 2, axis=1))
        potential_energy[1] = current_potential_energy
        kinetic_energy[1] = calculate_kinetic_energy(wrapped_positions, velocity)
        # t >=2:
        for t in range(2,config.steps):
            # print(f"Step {t}/{config.steps-1} at temperature {temperature} K")
            force, acceleration, current_potential_energy = compute_forces(wrapped_positions, box_size, two_potentials=config.two_potentials)
            next_positions, velocity = velocity_verlet(positions, velocity, acceleration, force)
            # velocity是中间时刻(1.5)的，positions和acceleration是当前时刻(2)的
            wrapped_positions = apply_boundary_conditions(next_positions, box_size)
            animation_pos[t] = wrapped_positions.copy()
            animation_vel[t] = velocity.copy()
            animation_acc[t] = acceleration.copy()
            # 计算各参数
            kinetic_energy[t] = calculate_kinetic_energy(wrapped_positions, velocity)
            potential_energy[t] = current_potential_energy
            msd_list[t] = np.mean(np.sum((next_positions - initial_positions) ** 2, axis=1))
            positions = next_positions

        # 分子动力学可视化
        animate_trajectory(animation_pos, config.nx, config.ny, config.nz, N_particles, config.steps, temperature) # TODO:debug
        # # 轨迹数据存储
        # with h5py.File("result\\simulation_data.h5", "w") as h5file:
        #     h5file.create_dataset("positions", data=animation_pos)
        #     h5file.create_dataset("velocities", data=animation_vel)
        #     h5file.create_dataset("accelerations", data=animation_acc)

        # 数据分析
        total_energy = kinetic_energy + potential_energy
        plot_energy(kinetic_energy, potential_energy, total_energy, temperature)
        plot_msd(msd_list, config.steps, temperature)
        diffusion_coefficient = compute_diffusion_coefficient_2(msd_list)
        # diffusion_coefficient = compute_diffusion_coefficient(msd_list, config.dt, interval=10)
        plot_diffusion_coefficient(diffusion_coefficient, config.steps, temperature)
        heat_capacity = calculate_heat_capacity(total_energy, temperature, N_particles)
        heat_capacity_list.append(heat_capacity)
        plot_kinetic_temperature(kinetic_energy, N_particles, temperature, config.steps)
        plot_rdf(animation_pos, temperature)

    # 跨温度的综合分析
    plot_heat_capacity(np.array(heat_capacity_list), config.temperature)

if __name__ == "__main__":
    main()