import math
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, x, y, rho, e, P):
        self.x = x
        self.y = y
        self.rho = rho
        self.e = e
        self.P = P

def Mesh(x1, x2, y1, y2, Nx_l, Nx_r, Ny):
    xdim = x2 - x1
    ydim = y2 - y1

    xstep_l = (xdim / 2.0) / Nx_l
    xstep_r = (xdim / 2.0) / Nx_r
    ystep = ydim / ((Ny - 1) * math.sqrt(3) / 2)

    mesh = []

    xOffset_l = xstep_l / 2
    xOffset_r = xstep_r / 2
    yOffset = ystep * math.sqrt(3) / 2

    for jj in range(1, Ny + 1):
        for ii in range(1, Nx_l + Nx_r + 1):
            y = y1 + (jj - 1) * yOffset

            if ii <= Nx_l: # Left side
                x = x1 + (ii - 1) * xstep_l + (jj % 2) * xOffset_l
                p = Particle(x, y, 1.0, 2.5, 1.0)
            else: # Right side
                ii_r = ii - Nx_l
                x = x1 + xdim / 2 + (ii_r - 1) * xstep_r + (jj % 2) * xOffset_r
                p = Particle(x, y, 0.2, 1.795, 0.1795)

            mesh.append(p)

    return mesh

# Ejemplo de uso
Nx_l = 20
Nx_r = 10
Ny = 10
mesh = Mesh(0, 10, 0, 3, Nx_l, Nx_r, Ny)

# Visualización
x_coords = [p.x for p in mesh]
y_coords = [p.y for p in mesh]

plt.scatter(x_coords, y_coords,s=2)
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Particle Mesh Distribution')
plt.show()
