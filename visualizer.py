import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import sys

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    print("Reading", file_path, "...")
    
    time_steps = []
    trajectories = {}
    
    for line in lines:
        if line != "\n":
            elements = line.strip().split('\t')
            time = float(elements[1])
            time_steps.append(time)
            
            positions = elements[2:]
            for i, pos in enumerate(positions):
                x, y, z = map(float, pos.split(','))
                if i not in trajectories:
                    trajectories[i] = {'x': [], 'y': [], 'z': []}
                trajectories[i]['x'].append(x)
                trajectories[i]['y'].append(y)
                trajectories[i]['z'].append(z)
    
    return time_steps, trajectories

def animate_trajectories(time_steps, trajectories):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    scatter_plots = []
    line_plots = []
    
    for idx_i, i in enumerate(trajectories):
        scatter, = ax.plot([], [], [], 'o')
        if idx_i <= 7:
            line, = ax.plot([], [], [], '--')
        else:
            line, = ax.plot([], [], [], '-')
        scatter_plots.append(scatter)
        line_plots.append(line)
    
    def init():
        ax.set_xlim(-10e8, 10e8)
        ax.set_ylim(-10e8, 10e8)
        ax.set_zlim(-10e8, 10e8)
        for scatter, line in zip(scatter_plots, line_plots):
            scatter.set_data([], [])
            scatter.set_3d_properties([])
            line.set_data([], [])
            line.set_3d_properties([])
        return scatter_plots + line_plots
    
    def update(frame):
        time_index = frame % len(time_steps)
        for i, (scatter, line) in enumerate(zip(scatter_plots, line_plots)):
            x = trajectories[i]['x'][:time_index+1]
            y = trajectories[i]['y'][:time_index+1]
            z = trajectories[i]['z'][:time_index+1]
            scatter.set_data(x[-1:], y[-1:])
            scatter.set_3d_properties(z[-1:])
            line.set_data(x, y)
            line.set_3d_properties(z)
        return scatter_plots + line_plots
    
    ani = FuncAnimation(fig, update, frames=len(time_steps), init_func=init, blit=True, repeat=True, interval=0.001)
    plt.show()

def main(sys_args):
    print("MINORBIT2 Visualization Utility")
    
    if len(sys_args) > 1:
        file_path = sys_args[1]

    else:
        file_path = input("Name of the result file to animate: ")

    try:
        time_steps, trajectories = read_data(file_path)
    except FileNotFoundError:
        print("File not found:", file_path)
        if not file_path.endswith(".txt"):
            file_path = file_path + ".txt"
        else:
            print("")
            main([])

    try:
        time_steps, trajectories = read_data(file_path)
    except FileNotFoundError:
        print("File not found:", file_path)
        print("")
        main([])
        
    animate_trajectories(time_steps, trajectories)

if __name__ == "__main__":
    main(sys.argv)

