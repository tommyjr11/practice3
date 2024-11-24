import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# 获取所有的CSV文件列表，按时间步排序
file_list = sorted(glob.glob('data/step_*.csv'))

# 初始化一个列表来存储所有时间步的rho数据
rho_data_list = []

# 检查是否找到数据文件
if not file_list:
    print("No data files found in the 'data' directory.")
    exit()

# 读取数据
for filename in file_list:
    with open(filename, 'r') as f:
        # 跳过时间信息行
        time_line = next(f)
        # 读取剩余的数据
        data = []
        for line in f:
            # 去除行末的换行符并分割
            values = line.strip().split(',')
            # 转换为浮点数
            try:
                float_values = [float(val) for val in values if val != '']
                data.append(float_values)
            except ValueError as e:
                print(f"Error parsing line in {filename}: {line}")
                continue
        # 将数据转换为numpy数组
        data = np.array(data)
        if data.size == 0:
            print(f"No data in file {filename}")
            continue
        # 计算每行的单元格数量
        num_cells_per_row = data.shape[1] // 4  # 每个单元格有4个值
        if data.shape[1] % 4 != 0:
            print(f"Data format error in file {filename}: columns not multiple of 4")
            continue
        # 重塑数据
        data = data.reshape((data.shape[0], num_cells_per_row, 4))
        # 提取rho（第一个维度是行，第二个是列，第三个是变量）
        rho = data[:, :, 0]
        v1 = data[:, :, 1]
        v2 = data[:, :, 2]
        p = data[:, :, 3]
        # 计算energy
        energy = 0.5*rho*(v1**2 + v2**2) + p/(1.4-1)
        rho_data_list.append(p)

# 检查是否成功读取rho数据
if not rho_data_list:
    print("No rho data extracted from files.")
    exit()

# 创建网格
ny, nx = rho_data_list[0].shape
X, Y = np.meshgrid(range(nx), range(ny))

# 创建图形和3D坐标轴
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 计算rho的全局最小值和最大值，用于设定Z轴范围
z_min = min(np.nanmin(rho) for rho in rho_data_list)
z_max = max(np.nanmax(rho) for rho in rho_data_list)

# 初始化绘图对象
surface = [ax.plot_surface(X, Y, rho_data_list[0], cmap='viridis')]

# 更新函数
def update_plot(frame):
    ax.clear()
    ax.set_zlim(z_min, z_max)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('rho')
    ax.set_title(f'Time Step: {frame}')
    # 绘制新的rho数据
    ax.plot_surface(X, Y, rho_data_list[frame], cmap='viridis')

# 创建动画
ani = animation.FuncAnimation(fig, update_plot, frames=len(rho_data_list), interval=200)

# 展示动画
plt.show()

# 绘制最后一个时间步的rho二维面图
rho_final = rho_data_list[-1]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, rho_final, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('rho')
ax.set_title('Density rho at the last time step')
plt.show()

# 选择中间的一行
mid_row = ny // 2
rho_line = rho_final[mid_row, :]

# 绘制折线图
plt.figure()
plt.plot(range(nx), rho_line, marker='o')
plt.title('Density rho at the last time step along the middle row')
plt.xlabel('X')
plt.ylabel('rho')
plt.grid(True)
plt.show()
