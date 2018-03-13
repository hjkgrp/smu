from matplotlib import pyplot as plt
import numpy as np
import math
from numpy.linalg import inv

x = np.arange(0,1,0.01)
ye = np.zeros(len(x))
Xe = np.array([0, 0.4, 0.6, 1])
yc = np.zeros(len(x))
y = np.zeros(len(Xe))

A = 0.5
B = 10
C = -5

for i in range(0, len(x)):
	ye[i] = (6 * x[i] - 2) ** 2 * math.sin(12 * x[i] - 4)
	yc[i] = A * ye[i] + B * (x[i] - 0.5) - C
	

for i in range(0, len(Xe)):
	y[i] = (6 * Xe[i] - 2) ** 2 * math.sin(12 * Xe[i] - 4)

Psi = np.zeros((len(Xe), len(Xe)))

for i in range(0, len(Xe)):
	for j in range(0, len(Xe)):
		Psi[i][j] = np.exp(-(Xe[i]-Xe[j]) ** 2)

muhat = np.matmul(np.matmul(np.ones(len(Xe)).transpose(), inv(Psi)), y) / np.matmul(np.matmul(np.ones(len(Xe)).transpose(), inv(Psi)), np.ones(len(Xe)))

b = np.matmul(inv(Psi), (y - np.ones(len(Xe)) * muhat))


psi = np.zeros((len(Xe), len(x)))
for i in range(0, len(Xe)):
	for j in range(0, len(x)):
		psi[i][j] = np.exp(-(x[j]-Xe[i]) ** 2)


plt.plot(yc)
plt.plot(ye)
plt.plot(muhat + np.matmul(b, psi))


plt.show()
