from matplotlib import pyplot as plt
import numpy as np
import math
from numpy.linalg import inv

x = np.arange(0,1,0.01)
Xe = np.array([0, 0.4, 0.6, 1])
Xc = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
yc_cont = np.zeros(len(x)) #these 'cont' are just for visualization, not the actual data points 
ye_cont = np.zeros(len(x))
ye = np.zeros(len(Xe))
yc = np.zeros(len(Xc))
yc_atXe = np.zeros(len(Xe))
d = np.zeros(len(Xe))
rho = 1.87
theta = 1 
A = 0.5
B = 10
C = -5

def one(length, T):
    if T == 0:
        one  = np.ones(len(length))
        return np.reshape(one, (len(one),1))
    if T == 1:
        oneT  = np.ones(len(length)).transpose()
        return np.reshape(oneT, (1,len(oneT)))

for i in range(0, len(x)):
	ye_cont[i] = (6 * x[i] - 2) ** 2 * math.sin(12 * x[i] - 4)
	yc_cont[i] = A * ye_cont[i] + B * (ye_cont[i] - 0.5) - C
	
for i in range(0, len(Xe)):
	ye[i] = (6 * Xe[i] - 2) ** 2 * math.sin(12 * Xe[i] - 4)

ye = np.reshape(ye, (len(ye), 1))

for i in range(0, len(Xc)):
        yc[i] = A * (6 * Xc[i] - 2) ** 2 * math.sin(12 * Xc[i] - 4) + B * ((6 * Xc[i] - 2) ** 2 * math.sin(12 * Xc[i] - 4) - 0.5) - C

yc = np.reshape(yc, (len(yc), 1))

for i in range(0, len(Xe)): #for the yc at the positions of Xe
        yc_atXe[i] = A * ye[i] + B * (ye[i] - 0.5) - C

yc_atXe = np.reshape(yc_atXe, (len(yc_atXe), 1))

y = np.vstack((np.reshape(yc, (len(yc),1)), np.reshape(ye, (len(ye),1))))

d = np.reshape(ye - rho * yc_atXe,(len(ye),1)) #both are columns now

def Psi(X1, X2):
	Psi = np.zeros((len(X1), len(X2)))
	
	for i in range(0, len(X1)):
		for j in range(0, len(X2)):
			Psi[i][j] = np.exp(-theta*(X1[i]-X2[j]) ** 2)
	return Psi

muc = np.matmul(np.matmul(one(Xc,1), inv(Psi(Xc,Xc))), yc) / (np.matmul(np.matmul(one(Xc,1), inv(Psi(Xc, Xc))), one(Xc,0)))

mud = np.matmul(np.matmul(one(Xe,1), inv(Psi(Xe,Xe))), d) / (np.matmul(np.matmul(one(Xe,1), inv(Psi(Xe, Xe))), one(Xe,0)))

s2c = np.matmul(np.matmul((yc - one(Xc,0) * muc).transpose(), inv(Psi(Xc, Xc))), (yc - one(Xc,0) * muc)) / len(Xc)

s2d = np.matmul(np.matmul((d - one(Xe,0) * mud).transpose(), inv(Psi(Xe, Xe))), (d - one(Xe,0) * mud)) / len(Xe)

C = np.vstack((np.hstack((s2c * Psi(Xc,Xc), rho * s2c * Psi(Xc, Xe))), np.hstack((rho * s2c * Psi(Xe, Xc), rho**2 * s2c * Psi(Xe, Xe) + s2d * Psi(Xe, Xe)))))

c = np.vstack((rho * s2c * Psi(Xc, x), rho**2 * s2c * Psi(Xe, x) +  s2d * Psi(Xe, x))) 

mu = np.matmul(np.matmul(one(C,1), inv(C)), y) / (np.matmul(np.matmul(one(C,1), inv(C)), one(C,0)))

new = mu = np.matmul(np.matmul(c.transpose(), inv(C)), (y - one(y,0) * mu))

plt.plot(new)
plt.show()
