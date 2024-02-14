import pandas as pd
import matplotlib.pyplot as plt

# Cargar los datos del archivo CSV
mesh = pd.read_csv('/home/yo/Documents/Tesis/codes/SPH/blast-Wave/mesh-2d.csv', delimiter=',')
print(mesh.head())

# Graficar la coordenada x en función de la densidad
plt.figure(figsize=(10, 10))
plt.scatter(mesh['x'], mesh['y'], s=1)
#plt.xlim(-10, 10)
#plt.ylim(-5,5)
plt.xlabel(' x')
plt.ylabel('y')
plt.show()

# Graficar la coordenada x en función de la velocidad
plt.figure(figsize=(10, 5))
plt.scatter(mesh['x'], mesh['rho'], s=1)
plt.xlabel('Posición x')
plt.ylabel('')
plt.show()
