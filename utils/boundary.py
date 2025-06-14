import numpy as np

def boundary_conditions(positions, box_size, supercell_vectors, inv_supercell):
    # 先转换为晶格坐标
    frac_coords = np.dot(positions[:, :3], inv_supercell)
    
    # 对每一维做模运算，实现周期性边界
    frac_coords = frac_coords % 1.0  # 自动限制在 [0, 1)
    
    # 转换回笛卡尔坐标
    wrapped_cartesian = np.dot(frac_coords, supercell_vectors.T)
    
    # 保留原来的时间维度（或其他额外列）
    wrapped_positions = np.copy(positions)
    wrapped_positions[:, :3] = wrapped_cartesian
    
    return wrapped_positions