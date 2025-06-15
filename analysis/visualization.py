import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def plot_initial_structure(positions):
    """3D plot of initial structure"""

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='black', s=40)

    bond_threshold = 1.7
    N = len(positions)
    for i in range(N):
        for j in range(i + 1, N):
            dist = np.linalg.norm(positions[i, :3] - positions[j, :3])
            if dist < bond_threshold:
                xs = [positions[i, 0], positions[j, 0]]
                ys = [positions[i, 1], positions[j, 1]]
                zs = [positions[i, 2], positions[j, 2]]
                ax.plot(xs, ys, zs, color='gray', linewidth=1)

    ax.set_title("Initial Graphite Structure", fontsize=14)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    plt.tight_layout()
    plt.show()


def animate_trajectory(nx, ny, nz, postions):
    """Animate the trajectory"""

    particle_ids = positions[:, 0].astype(int)
    positions = positions[:, 1:4]
    times = positions[:, 4]

    def func(i, j, k):
        return i * (ny * nz) + j * (nz) + k

    # 只选择边界内的原子，否则边界处横跳干扰视线
    def select_in_boundary():
        in_boundary = []
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                for k in range(1, nz + 1):
                    in_boundary.append(func(i, j, k))
        particle_ids = particle_ids[in_boundary]
        positions = positions[in_boundary]
        times = times[in_boundary]
    #select_in_boundary()

    unique_times = np.unique(times)
    unique_ids = np.unique(particle_ids)
    N_particles = len(unique_ids)
    N_steps = len(unique_times)

    # 为每个时间步准备数据索引
    # 把数据重塑成 (steps, N_particles, 3)
    positions_reshaped = positions.reshape((N_steps, N_particles, 3))

    # 画图初始化
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # 给边界处原子染色，观察零温下原子运动是否符合微扰条件
    min_atom_indices = func(nx//2, ny//2, nz//2)
    max_atom_indices = nx - 1
    colors = np.full(N_particles, 'black')
    colors[min_atom_indices] = 'red'
    colors[max_atom_indices] = 'green'

    scat = ax.scatter(np.zeros(N_particles), np.zeros(N_particles), np.zeros(N_particles), s=30, c=colors)
    # scat = ax.scatter([], [], [], s=30)

    ax.set_xlim(np.min(positions[:,0]), np.max(positions[:,0]))
    ax.set_ylim(np.min(positions[:,1]), np.max(positions[:,1]))
    ax.set_zlim(np.min(positions[:,2]), np.max(positions[:,2]))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Graphite Simulation')

    def update(frame):
        pos = positions_reshaped[frame]
        scat._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
        ax.set_title(f"Time = {(unique_times[frame]*1e12):.3f} ps")
        return scat,

    ani = FuncAnimation(fig, update, frames=N_steps, interval=100, blit=False)

    plt.show()