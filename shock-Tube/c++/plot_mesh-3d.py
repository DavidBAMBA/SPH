import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

n = sys.argv[1]
mesh = pd.read_csv('/home/yo/Documents/Tesis/codes/SPH/Toy-star/mesh-3d.csv', delimiter=',')
print(mesh.head())

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(mesh['x'], mesh['y'], mesh['z'], c=mesh['rho'])  # Coloreado por densidad
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.set_title('Visualización de la Malla SPH')
plt.savefig(f'{n}.png')
#plt.show()

# Leer los datos del archivo CSV
df = mesh

# Graficar la coordenada x en función de la presión
#plt.figure(figsize=(10, 5))
#plt.scatter(df['x'], df['rho'])
#plt.title('Posición x vs Presión')
#plt.xlabel('Posición x')
#plt.ylabel('Presión')
#plt.show()

