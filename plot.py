import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

a = 1.0
b = 1.0
Nx = 50
Ny = 100
data = "data_gauss-seidel_para.csv"

Z = np.loadtxt(data, delimiter=",")
print(Z)

y = np.linspace(0,b,Nx+2)
x = np.linspace(0,a,Ny+2)
X,Y = np.meshgrid(x,y)


# Tracé
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("f(x,y)")

plt.show()