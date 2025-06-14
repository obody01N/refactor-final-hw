import numpy as np
from utils.config import config

def generate_graphite(nx, ny, nz, basis_vectors, basis_atoms):
    positions = []
    boundary = []
    a1, a2, a3 = basis_vectors
    basis_frac = basis_atoms
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                shift = i*a1 + j*a2 + k*a3
                for l in range(len(basis_frac)):
                    pos = basis_frac[l][0] * a1 + basis_frac[l][1] * a2 + basis_frac[l][2] * a3
                    xyz = pos + shift
                    positions.append(np.append(xyz, 0.0))  # 最后一位是占位时间
                    # 检查是否为边界原子
                    is_boundary = False
                    # 检查x方向边界
                    if i == 0 or i == nx-1:
                        is_boundary = True
                    # 检查y方向边界
                    if j == 0 and (l == 0 or l == 1):
                        is_boundary = True
                    if j == ny-1 and (l == 2 or l == 3):
                        is_boundary = True
                    # 检查z方向边界
                    if k == 0 and (l == 0 or l == 2):
                        is_boundary = True
                    if k == nz-1 and (l == 1 or l == 3):
                        is_boundary = True
                    
                    if is_boundary:
                        boundary.append((atom_index,i,j,k,l))
                    
                    atom_index += 1
    
    return np.array(positions), boundary

def calculate_box_size(nx, ny, nz):
    return np.array([2.456, 2.456*np.sin(np.pi*2/3), 7.0]) * np.array([nx, ny, nz])

def supercell(nx, ny, nz, basis_vectors):
    supercell_vectors = basis_vectors * np.array([nx, ny, nz])
    inv_supercell = np.linalg.inv(supercell_vectors.T)
    return supercell_vectors, inv_supercell