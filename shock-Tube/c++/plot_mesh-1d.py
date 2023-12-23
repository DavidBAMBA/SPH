import pandas as pd
import matplotlib.pyplot as plt

# Cargamos los datos
mesh = pd.read_csv('/home/yo/Documents/Tesis/codes/SPH/SPH/shock-Tube/c++/mesh-1d.csv', delimiter=',')

# Creamos una figura y un eje
fig, ax = plt.subplots()

# Graficamos la posición 'x' en función de la presión 'P'
ax.scatter(mesh['x'], mesh['P'], c='blue', label='Presión')

# Establecemos las etiquetas y el título
ax.set_xlabel('Posición x')
ax.set_ylabel('Presión')
ax.set_xlim(0.2, 1.0)
ax.set_ylim(0.0, 1.0)
ax.set_title('Posición x vs Presión')

# Mostramos la leyenda
ax.legend()

# Mostramos la gráfica
plt.savefig('plot-1d')
plt.show()
