import numpy as np

def apply_boundary_conditions(positions, box_size):
    """
    Apply periodic boundary conditions for a cubic lattice.

    Parameters:
    positions : np.ndarray
        Array of shape (N_particles, 3) representing particle positions.
    box_size : np.ndarray
        Array of shape (3,) representing the dimensions of the cubic box [Lx, Ly, Lz].

    Returns:
    np.ndarray
        Adjusted positions within the periodic boundary conditions.
    """
    return positions - np.round(positions / box_size) * box_size

def apply_boundary_conditions_2(positions, supercell_vectors, inv_supercell):
    # 适用于任意晶格
    # rij = rij - np.round(rij / box_size) * box_size  # PBC处理
    return positions - np.round(positions @ inv_supercell) @ supercell_vectors
