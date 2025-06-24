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


def animate_trajectory(animation_pos, nx, ny, nz, N_particles, N_steps, temperature):
# def animate_trajectory(animation_pos, nx, ny, nz, N_particles, N_steps, unique_times):
    """Animate the trajectory with time labels"""

    def func(i, j, k):
        return i * (ny * nz) + j * nz + k

    # Initialize the plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Highlight specific atoms for visualization
    min_atom_index = func(nx // 2, ny // 2, nz // 2)
    max_atom_index = func(nx - 1, ny - 1, nz - 1)
    colors = np.full(N_particles, 'black')
    colors[min_atom_index] = 'red'
    colors[max_atom_index] = 'green'

    scat = ax.scatter(np.zeros(N_particles), np.zeros(N_particles), np.zeros(N_particles), s=30, c=colors)

    ax.set_xlim(np.min(animation_pos[:, :, 0]), np.max(animation_pos[:, :, 0]))
    ax.set_ylim(np.min(animation_pos[:, :, 1]), np.max(animation_pos[:, :, 1]))
    ax.set_zlim(np.min(animation_pos[:, :, 2]), np.max(animation_pos[:, :, 2]))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Graphite Simulation')

    def update(frame):
        pos = animation_pos[frame]
        scat._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        ax.set_title(f"Time = {frame * 1e-2:.2f} ps")  # Adjusted time calculation for simplicity
        return scat,

    ani = FuncAnimation(fig, update, frames=N_steps, interval=100, blit=False)

    ani.save(f"result\\{temperature}K\\trajectory_animation.gif", writer="pillow")  # Save the animation as a GIF file
    plt.show()

    # def update(frame):
    #     pos = positions_reshaped[frame]
    #     scat._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
    #     ax.set_title(f"Time = {(unique_times[frame] * 1e12):.3f} ps")
    #     return scat,

    # ani = FuncAnimation(fig, update, frames=N_steps, interval=100, blit=False)

    # plt.show()

    
# def animate_trajectory(animation_pos, nx, ny, nz, N_particles, N_steps):
#     """Animate the trajectory"""

#     # positions 是[id,x,y,z,t]的形式存储的
#     particle_ids = np.arange(N_particles)
#     # particle_ids = animation_pos

#     def func(i, j, k):
#         return i * (ny * nz) + j * (nz) + k

#     # positions_reshaped = positions.reshape((N_steps, N_particles, 3))
#     positions_reshaped = animation_pos

#     # 画图初始化
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     # 给边界处原子染色，观察零温下原子运动是否符合微扰条件
#     min_atom_indices = func(nx//2, ny//2, nz//2)
#     max_atom_indices = nx - 1
#     colors = np.full(N_particles, 'black')
#     colors[min_atom_indices] = 'red'
#     colors[max_atom_indices] = 'green'

#     scat = ax.scatter(np.zeros(N_particles), np.zeros(N_particles), np.zeros(N_particles), s=30, c=colors)
#     # scat = ax.scatter([], [], [], s=30)

#     ax.set_xlim(np.min(positions[:,0]), np.max(positions[:,0]))
#     ax.set_ylim(np.min(positions[:,1]), np.max(positions[:,1]))
#     ax.set_zlim(np.min(positions[:,2]), np.max(positions[:,2]))
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.set_title('Graphite Simulation')

#     def update(frame):
#         pos = positions_reshaped[frame]
#         scat._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
#         ax.set_title(f"Time = {(unique_times[frame]*1e12):.3f} ps")
#         return scat,

#     ani = FuncAnimation(fig, update, frames=N_steps, interval=100, blit=False)

#     plt.show()