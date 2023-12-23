import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    return np.loadtxt(filename, delimiter=',')

def plot_system(S_i):
    """ Plot the state of the system """
    S_plot = S_i
    
    fig, axs = plt.subplots(2, 2, figsize=[10, 10])  # Creates a 2x2 grid

    # Density plot
    axs[0, 0].scatter(S_plot[:,0], S_plot[:,2])
    axs[0, 0].set_title('Density plot')
    axs[0, 0].set_xlabel('$Position ')
    axs[0, 0].set_xlim([0.2, 1.0])
    axs[0, 0].set_ylim([0, 1])
    axs[0, 0].set_ylabel('$Density ')

    # Velocity plot
    axs[0, 1].scatter(S_plot[:,0], S_plot[:,1])
    axs[0, 1].set_title('Velocity plot')
    axs[0, 1].set_xlabel('$Position ')
    axs[0, 0].set_xlim([0.2, 1.0])
    axs[0, 1].set_ylim([0, 1])
    axs[0, 1].set_ylabel('$Velocity ')

    # Energy plot
    axs[1, 0].scatter(S_plot[:,0], S_plot[:,3])
    axs[1, 0].set_title('Energy plot')
    axs[1, 0].set_xlabel('$Position ')
    axs[0, 0].set_xlim([0.2, 1.0])
    axs[1, 0].set_ylim([1.8, 2.6])
    axs[1, 0].set_ylabel('$Internal energy ')

    # Pressure plot
    axs[1, 1].scatter(S_plot[:,0], S_plot[:,4])
    axs[0, 0].set_xlim([0.2, 1.0])
    axs[1, 1].set_title('Pressure plot')
    axs[1, 1].set_xlabel('$Position $')
    axs[1, 1].set_ylim([0, 1.2])
    axs[1, 1].set_ylabel('$Pressure $')

    plt.savefig('c-plot.png', dpi=500)  # Guarda la gráfica como una imagen PNG con una resolución de 300 dpi
    plt.tight_layout()  # Adjusts the space between the plots for better layout
    plt.show()


S_i = load_data("output.csv")
plot_system(S_i)


