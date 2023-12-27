import matplotlib.pyplot as plt
import pandas as pd


def plot_system(S_i):
    """ Plot the state of the system """
    S_plot = S_i
    
    fig, axs = plt.subplots(2, 2, figsize=[10, 10])  # Creates a 2x2 grid
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    # Density plot
    axs[0, 0].scatter(S_plot['x'], S_plot['rho'], s=1)
    axs[0, 0].set_xlabel('x')
    #axs[0, 0].set_xlim([0.2, 1.0])
    axs[0, 0].set_ylim([0, 1.1])
    axs[0, 0].set_ylabel('Density')

    # Velocity plot
    axs[0, 1].scatter(S_plot['x'], S_plot['vx'], s=1)
    axs[0, 1].set_xlabel('x')
    #axs[0, 1].set_xlim([0.2, 1.0])
    axs[0, 1].set_ylim([0, 1])
    axs[0, 1].set_ylabel('Velocity ')

    # Energy plot
    axs[1, 0].scatter(S_plot['x'], S_plot['e'], s=1)
    axs[1, 0].set_xlabel('x ')
    #axs[1, 0].set_xlim([0.2, 1.0])
    axs[1, 0].set_ylim([1.8, 2.6])
    axs[1, 0].set_ylabel('Energy')

    # Pressure plot
    axs[1, 1].scatter(S_plot['x'], S_plot['P'], s=1)
    axs[1, 1].set_xlabel('x')
    #axs[1, 1].set_xlim([0.2, 1.1])
    axs[1, 1].set_ylim([0, 1.2])
    axs[1, 1].set_ylabel('Pressure')
    
    plt.savefig('c-plot.png', dpi=500)  # Guarda la gráfica como una imagen PNG con una resolución de 300 dpi
    plt.tight_layout()  # Adjusts the space between the plots for better layout
    plt.show()


mesh = pd.read_csv('/home/yo/Documents/Tesis/codes/SPH/shock-Tube/c++/mesh-1d.csv', delimiter=',')
plot_system(mesh)


