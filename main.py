import numpy as np
from utils.config import Config
from utils.lattice import generate_graphite, calculate_box_size, calculate_supercell_vectors
from utils.constants import v_unit, a_unit, amu, angstrom, eV
from utils.__init__ import initialize_velocity_distribution
from utils.evolution import velocity_verlet
from utils.boundary import apply_boundary_conditions
from utils.potentials import compute_forces
from analysis.visualization import plot_initial_structure, animate_trajectory
from analysis.heat_capacity import calculate_heat_capacity, plot_heat_capacity
from analysis.energy import calculate_kinetic_energy, plot_energy
from analysis.msd import compute_diffusion_coefficient, compute_diffusion_coefficient_2, plot_msd, plot_diffusion_coefficient
from analysis.temperature import plot_kinetic_temperature
from analysis.rdf import plot_rdf

config = Config()

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

    # Initialize system using config parameters
    # positions: N x 3 array of atom positions in Angstrom
    N_particles = config.nx * config.ny * config.nz * len(config.basis_atoms)
    animation_pos = np.zeros((config.steps, N_particles, 3))
    positions, boundary = generate_graphite(
        config.nx, config.ny, config.nz,
        basis_vectors=config.basis_vectors,
        basis_atoms=config.basis_atoms
    )
    box_size = calculate_box_size(config.nx, config.ny, config.nz)
    supercell_vectors, inv_supercell = calculate_supercell_vectors(config.nx, config.ny, config.nz)
    positions = apply_boundary_conditions(positions, box_size, supercell_vectors, inv_supercell)
    animation_pos[0] = positions.copy()
    initial_positions = positions.copy()
    plot_initial_structure(positions)
    msd_list = np.zeros(config.steps)  # List to store mean square displacement at each time step
    kinetic_energy = np.zeros(config.steps)  # List to store kinetic energy at each time step
    potential_energy = np.zeros(config.steps)  # List to store potential energy at each time step
    
    # simulation at different temperatures
    for temperature in config.temperature:
        print(f"Running simulation at temperature: {temperature} K")
        # t = 0:
        initial_velocity = initialize_velocity_distribution(temperature) # v(t=0)
        initial_velocity *= v_unit  # Convert to angstrom / dt
        force, current_potential_energy = compute_forces(positions, box_size)
        acceleration = force / config.mass * a_unit
        msd_list[0] = 0
        potential_energy[0] = current_potential_energy
        # t = 1:
        positions = positions + initial_velocity + 0.5 * acceleration
        wrapped_positions = apply_boundary_conditions(positions, box_size)
        force, current_potential_energy = compute_forces(wrapped_positions, box_size)
        acceleration = force / config.mass * a_unit
        velocity += 0.5 * acceleration  # v(t=0.5)
        animation_pos[1] = wrapped_positions.copy()
        msd_list[1] = np.mean(np.sum((wrapped_positions - initial_positions) ** 2, axis=1))
        # t >=2:
        for t in range(2,config.steps):
            force, potential_energy = compute_forces(wrapped_positions, box_size)
            next_positions, velocity, acceleration = velocity_verlet(positions, velocity, acceleration, force)
            # velocity是中间时刻(1.5)的，positions和acceleration是当前时刻(2)的
            wrapped_positions = apply_boundary_conditions(next_positions, box_size)
            animation_pos[t] = wrapped_positions.copy()
            # 计算各参数
            kinetic_energy = calculate_kinetic_energy(wrapped_positions, velocity)
            msd_list[t] = np.mean(np.sum((next_positions - initial_positions) ** 2, axis=1))
            positions = next_positions

        # analysis
        total_energy = kinetic_energy + potential_energy
        plot_energy(kinetic_energy, potential_energy, total_energy)
        plot_msd(msd_list, config.steps)
        diffusion_coefficient = compute_diffusion_coefficient_2(msd_list, config.dt)
        # diffusion_coefficient = compute_diffusion_coefficient(msd_list, config.dt, interval=10)
        plot_diffusion_coefficient(diffusion_coefficient, config.steps)
        plot_heat_capacity(heat_capacity, config.temperature)
        plot_kinetic_temperature(kinetic_energy, N_particles, temperature)
        plot_rdf(animation_pos)

        animate_trajectory(positions, box_size, config.steps)

    # 跨温度的综合分析
    heat_capacity = calculate_heat_capacity(total_energy, config.temperature)

    # print("mean square distance each t: ", msd_list)
    # D, slope, intercept, r2 = compute_diffusion_coefficient(msd_list, config.dt, interval=10)
    # print(f"Diffusion coefficient D = {D:.4e} (fit R^2 = {r2:.4f})")

if __name__ == "__main__":
    main()