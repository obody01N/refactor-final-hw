import numpy as np
from scipy.stats import linregress

def msd(system, initial_positions):
    N = system.N
    r = system.r
    msd_value = 0.0
    for i in range(N):
        dr = r[i] - initial_positions[i]
        msd_value += np.dot(dr, dr)
    return msd_value / N

def compute_diffusion_coefficient(msd_list, dt, interval, fit_start_ratio=0.5):
    """
    使用线性拟合计算扩散系数 D
    msd_list: MSD 序列
    dt: 单步时间间隔
    interval: 每多少步存一次 msd
    fit_start_ratio: 从后半段开始拟合（避免初始非线性段）
    """
    times = np.array([i * dt * interval for i in range(len(msd_list))])
    msd_array = np.array(msd_list)

    start = int(len(times) * fit_start_ratio)  # 拟合起点
    slope, intercept, r_value, p_value, std_err = linregress(times[start:], msd_array[start:])

    D = slope / 6.0
    return D, slope, intercept, r_value**2