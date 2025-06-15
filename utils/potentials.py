import numpy as np

def smooth_cutoff(r, ro, rc):
    """Cubic smooth cutoff function for LJ potential between ro and rc."""
    if r < ro:
        return 1.0
    elif r < rc:
        x = (rc - r) / (rc - ro)
        return 6 * x**5 - 15 * x**4 + 10 * x**3
    else:
        return 0.0

def lj_force_scalar(r, epsilon, sigma, ro, rc):
    """Lennard-Jones potential"""
    r6 = (sigma / r) ** 6
    r12 = r6 ** 2
    force_scalar = 24 * epsilon * (2*r12 - r6) / r**2 * smooth_cutoff(r, ro, rc)
    return force_scalar

def morse_force_scalar(r, D_e, a, r_e):
    """Morse potential"""
    exp_term = np.exp(-a * (r - r_e))
    force_scalar = 2 * a * D_e * (1 - exp_term) * exp_term
    return force_scalar

def lj_potential(r, epsilon, sigma, ro, rc):
    """Lennard-Jones potential with smooth cutoff."""
    if r < ro:
        return 0.0
    elif r < rc:
        r6 = (sigma / r) ** 6
        r12 = r6 ** 2
        potential = 4 * epsilon * (r12 - r6) * smooth_cutoff(r, ro, rc)
        return potential
    else:
        return 0.0
    
def morse_potential(r, D_e, a, r_e):
    """Morse potential without smooth cutoff."""
    exp_term = np.exp(-a * (r - r_e))
    potential = D_e * (1 - exp_term) ** 2
    return potential

def compute_forces(positions, box_size, **kwargs): 
    #  boundary,
    """Calculate forces using appropriate potentials"""

    # 提取参数
    cutoff = kwargs.get('cutoff', 3.0)  # 默认截断距离
    cutoff_face = kwargs.get('cutoff_face', 1.0)  # 面内截断距离
    two_potential = kwargs.get('two_potential', False)  # 是否使用两种势能
    epsilon = kwargs.get('epsilon', 0.0103)  # eV
    sigma = kwargs.get('sigma', 3.4)  # Angstrom
    a = kwargs.get('a', 1.0)  # Morse势参数
    D_e = kwargs.get('D_e', 0.1)  # eV
    r_e = kwargs.get('r_e', 1.0)  # 平衡距离
    ro = kwargs.get('ro', 0.5)  # 截断距离
    rc = kwargs.get('rc', 3.0)  # 截断距离
    mass = kwargs.get('mass', 12.0)  # amu, 碳原子质量

    # 计算加速度和力
    N = len(positions)
    accelerations = np.zeros((N, 3))
    forces = np.zeros((N, N, 3))  # 存储两两之间的力
    current_potential_energy = 0.0

    for i in range(N):
        for j in range(i + 1, N):
            rij = positions[i] - positions[j]
            rij = rij - np.round(rij / box_size) * box_size  # PBC处理
            r = np.linalg.norm(rij)
            if r < cutoff:
                if np.abs(rij[2]) < cutoff_face and two_potential == True: # 面内共价键用Morse势近似
                    force_scalar = morse_force_scalar(r, D_e, a, r_e)
                    potential = morse_potential(r, D_e, a, r_e)
                else: # 面间相互作用用LJ势描述
                    force_scalar = lj_force_scalar * smooth_cutoff(r, ro, rc)
                    potential = lj_potential(r, epsilon, sigma, ro, rc)
                current_potential_energy += potential
                force_vector = force_scalar * rij
                forces[i, j] = force_vector
                forces[j, i] = -force_vector
                accelerations[i] += force_vector / mass
                accelerations[j] -= force_vector / mass
    return forces, potential # eV/Å