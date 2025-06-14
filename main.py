# main.py
import numpy as np
from utils.config import config
from utils.lattice import generate_graphite, calculate_box_size
from utils.constants import v_unit, a_unit
from utils.evolution import velocity_verlet
from utils.boundary import apply_boundary_conditions
from utils.potentials import compute_forces
from analysis.visualization import plot_initial_structure, animate_trajectory
from analysis.heat_capacity import calculate_heat_capacity
from analysis.energy import calculate_energies
from analysis.msd import compute_diffusion_coefficient

def main():
    # Initialize system using config parameters
    positions, boundary = generate_graphite(
        config.nx, config.ny, config.nz,
        basis_vectors=config.basis_vectors,
        basis_atoms=config.basis_atoms
    )
    box_size = calculate_box_size(config.nx, config.ny, config.nz)
    positions = apply_boundary_conditions(positions, box_size)
    initial_positions = positions.copy()
    plot_initial_structure(positions)
    
    # simulation at different temperatures
    for temperature in [11.0, 300.0]:
        print(f"Running simulation at temperature: {temperature} K")
        _, energy, msd = simulation(temperature=config.temperature, dt=config.dt,
               positions=positions, initial_positions=initial_positions, boundary=boundary, box_size=box_size)
    # analysis
    calculate_heat_capacity()
    



def simulation(temperature, dt, positions, initial_positions, boundary, box_size):
    """Run a molecular dynamics simulation of graphite at a given temperature(Kelvin).

    positions - angstrom
    velocities - m/s = v_unit angstrom / dt
    dt - seconds
    forces - eV/angstrom
    acceleration - eV/(angstrom * amu) = a_unit angstrom / dt^2
    v_unit = 1/(angstrom / dt) = 1e-4
    a_unit = eV/(angstrom*amu) / (angstrom / dt^2) = 0.9648533215665326
    """

    # Initialize system
    msd_list = []

    # Simulation loop
    for step in range(config.steps):
        forces = compute_forces(positions, boundary)
        positions = velocity_verlet(positions, forces, dt)
        energy = calculate_energies(positions, forces, boundary)
        msd_list.append(np.mean(np.sum((positions - initial_positions) ** 2, axis=1)))
        positions = apply_boundary_conditions(positions, box_size)

    animate_trajectory(positions, box_size, step)

    print("mean square distance each t: ", msd_list)
    D, slope, intercept, r2 = compute_diffusion_coefficient(msd_list, dt, interval=10)
    print(f"Diffusion coefficient D = {D:.4e} (fit R^2 = {r2:.4f})")

    return positions, energy

def main():


if __name__ == "__main__":
    main()