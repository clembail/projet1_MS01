import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

a = 1.0
b = 1.0
Nx = 100
Ny = 98
data = "data_jacobi_para.csv"

Z = np.loadtxt(data, delimiter=",")
print(Z)

x = np.linspace(0,a,Nx+2)
y = np.linspace(0,b,Ny+2)
X,Y = np.meshgrid(x,y)


# Tracé
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("f(x,y)")

plt.show()