import numpy as np

def apply_boundary_conditions(positions, supercell_vectors, inv_supercell):
    supercell_vectors = supercell_vectors.T
    inv_supercell = inv_supercell.T
    
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

# # pbc修正位置和位移
# def wrap(r,H,unwrap_s=None, s_now=None, s_prev=None):
    
#     s = (Hinv @ r.T).T      # 实坐标转分数坐标
#     s_wrapped = s % 1.0     # 所有分量进入 [0,1)
#     r_wrapped = (H @ s_wrapped.T).T  # 再转回实空间

#     # s_now, s_prev 是相邻两帧的分数坐标
#     delta_s = s_now - s_prev
#     delta_s -= np.round(delta_s)  # 减去整数晶胞跳跃数
#     unwrap_s += delta_s
#     unwrap_r = (H @ unwrap_s.T).T
