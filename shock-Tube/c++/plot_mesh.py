import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

mesh = pd.read_csv('/home/yo/Documents/Tesis/codes/SPH/SPH/shock-Tube/c++/mesh.csv', delimiter=',')
print(mesh.head())

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(mesh['x'], mesh['y'], mesh['z'], c=mesh['rho'])  # Coloreado por densidad
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Visualización de la Malla SPH')
plt.show()


# Leer los datos del archivo CSV
df = pd.read_csv("mesh.csv")

# Graficar la coordenada x en función de la presión
plt.figure(figsize=(10, 5))
plt.scatter(df['x'], df['P'], s=1, alpha=0.5)
plt.title('Posición x vs Presión')
plt.xlabel('Posición x')
plt.ylabel('Presión')
plt.grid(True)
plt.show()

# Graficar la coordenada x en función de la energía
plt.figure(figsize=(10, 5))
plt.scatter(df['x'], df['e'], s=1, alpha=0.5)
plt.title('Posición x vs Energía')
plt.xlabel('Posición x')
plt.ylabel('Energía')
plt.grid(True)
plt.show()

