import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

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
    times = np.array([i * dt * interval for i in range(1,len(msd_list))])
    msd_array = np.array(msd_list)

    start = int(len(times) * fit_start_ratio)  # 拟合起点
    slope, intercept, r_value, p_value, std_err = linregress(times[start:], msd_array[start:])

    D = slope / 6.0
    return D, slope, intercept, r_value**2

def compute_diffusion_coefficient_2(msd_list):
    time_steps = np.arange(len(msd_list))
    diffusion_coefficient = msd_list[1:] / (6 * time_steps[1:])  # 6D for 3D space
    return diffusion_coefficient # MSD/10fs

def plot_msd(msd_list, steps, temperature):
    """绘制 MSD 曲线"""
    times = np.arange(0, steps, 1)
    # plt.figure(figsize=(10, 6))
    plt.plot(times, msd_list, label='MSD', color='blue')
    plt.xlabel('Time (10fs)')
    plt.ylabel('Mean Square Displacement (Å²)')
    plt.title('Mean Square Displacement vs Time')
    plt.legend()
    plt.grid()
    plt.savefig(f"result\\{temperature}K\\MSD VS time (10fs).png")
    plt.show()

def plot_diffusion_coefficient(diffusion_coefficient, steps, temperature):
    """绘制扩散系数随时间变化的图像"""
    times = np.arange(1, steps, 1)
    # plt.figure(figsize=(10, 6))
    plt.plot(times, diffusion_coefficient, label='Diffusion Coefficient', color='red')
    plt.xlabel('Time (10fs)')
    plt.ylabel('Diffusion Coefficient (0.1Å²/fs)')
    plt.title('Diffusion Coefficient vs Time')
    plt.legend()
    plt.grid()
    plt.savefig(f"result\\{temperature}K\\Diffusion Coefficient VS time (10fs).png")
    plt.show()