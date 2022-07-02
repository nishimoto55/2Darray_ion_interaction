# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 13:20:33 2022

@author: ni4mo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

 # -*constant for drawing graph *- #
AxesSize = 24
LabelSize = 24
FontSize = 28
MS = 5.0
LineWidth = 1.5
spines = 2.0
length = 3
width = 1.5
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'


# Charge
Q = 1.6 * 10 ** (-19) # [C]
# mass
m0 = 40 * 1.67 * 10 ** (-27) # [kg]
# Permittivity
epsilon0 = 8.85 * 10 ** (-12) # [F/m]
# pi
pi = np.pi

mu_meter = 10 ** (-6)
m_meter = 10 ** (-3)

# Number of iom at a string
num = 2

# distance between two ions
def L(x_1, x_2, z_1, z_2):
    return np.sqrt((x_1 - x_2)**2 + (z_1 - z_2)**2)

# Distance between center-ion (equibrium position)
def l(omega_z):
    return (1/2)**(2/3) * ((Q**2)/(4 * pi * epsilon0 * m0 * (2*pi*omega_z)**2))**(1/3)

omega_x = 0.946 * 10 ** 6 # [Hz]
omega_y = 1.03 * 10 ** 6 # [Hz]
omega_z = 0.19 * 10 ** 6 # [Hz]
z0 = l(omega_z)
x0 = 30 * mu_meter
y0 = 200 * mu_meter

# initialaize ion position
zpos = np.zeros(2*num).reshape(2*num)
xpos = np.zeros(2*num).reshape(2*num)
ypos = np.zeros(2*num).reshape(2*num)

#create ion position
for i in range(2*num):
    ypos[i] = y0

for i in range(num):
    xpos[i] = x0

for i in range(num):
    xpos[num+i] = -x0

j = 0
for i in range(num):
    zpos[i + j] = -z0
    zpos[i + j + 1] = z0
    j = j + 1
    

print("xpos = \n", xpos,"\nzpos = \n", zpos, "\n ypos\n", ypos)


#initialize K_xx
K_xx = np.zeros(16).reshape(4,4)
K_yy = np.zeros(16).reshape(4,4)
K_zz = np.zeros(16).reshape(4,4)

for i in range(2*num):
    S = 0
    for j in range(2*num):
        if i==j:
            t = 0
        else:
            t_1 = (-1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
            t_2 = (3*(xpos[i] - xpos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
            t_3 = (-1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
            t_4 = (3*(xpos[j] - xpos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
            t = t_1 + t_2 + t_3 + t_4
        S = S + t
    K_xx[i,i] = m0*(2 * pi * omega_x)**2 + (Q**2)/(8*pi*epsilon0)*S

for i in range(2*num):
    S = 0
    for j in range(i+1,2*num):
        t_1 = (1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
        t_2 = ((-3)*(xpos[i] - xpos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
        t_3 = (1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
        t_4 = ((-3)*(xpos[j] - xpos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
        S = t_1 + t_2 + t_3 + t_4
        K_xx[i,j] = (Q**2)/(8 * pi * epsilon0) * S
        K_xx[j,i] = K_xx[i,j]

print('K_xx = \n',K_xx)

for i in range(2*num):
    S = 0
    for j in range(2*num):
        if i==j:
            t = 0
        else:
            t_1 = (-1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
            t_2 = (3*(ypos[i] - ypos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
            t_3 = (-1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
            t_4 = (3*(ypos[j] - ypos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
            t = t_1 + t_2 + t_3 + t_4
        S = S + t
    K_yy[i,i] = m0*(2 * pi * omega_y)**2 + (Q**2)/(8*pi*epsilon0)*S

for i in range(2*num):
    S = 0
    for j in range(i+1,2*num):
        t_1 = (1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
        t_2 = ((-3)*(ypos[i] - ypos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
        t_3 = (1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
        t_4 = ((-3)*(ypos[j] - ypos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
        S = t_1 + t_2 + t_3 + t_4
        K_yy[i,j] = (Q**2)/(8 * pi * epsilon0) * S
        K_yy[j,i] = K_yy[i,j]

print('K_yy = \n',K_yy)

for i in range(2*num):
    S = 0
    for j in range(2*num):
        if i==j:
            t = 0
        else:
            t_1 = (-1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
            t_2 = (3*(zpos[i] - zpos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
            t_3 = (-1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
            t_4 = (3*(zpos[j] - zpos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
            t = t_1 + t_2 + t_3 + t_4
        S = S + t
    K_zz[i,i] = m0*(2 * pi * omega_z)**2 + (Q**2)/(8*pi*epsilon0)*S

for i in range(2*num):
    S = 0
    for j in range(i+1,2*num):
        t_1 = (1)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**3
        t_2 = ((-3)*(zpos[i] - zpos[j])**2)/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
        t_3 = (1)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**3    
        t_4 = ((-3)*(zpos[j] - zpos[i])**2)/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
        S = t_1 + t_2 + t_3 + t_4
        K_zz[i,j] = (Q**2)/(8 * pi * epsilon0) * S
        K_zz[j,i] = K_zz[i,j]

print('K_zz = \n',K_zz)

K_xy = np.zeros(16).reshape(4,4)
K_yx = np.zeros(16).reshape(4,4)
K_zy = np.zeros(16).reshape(4,4)
K_yz = np.zeros(16).reshape(4,4)

K_xz = np.zeros(16).reshape(4,4)
K_zx = np.zeros(16).reshape(4,4)

for i in range(2*num):
    S = 0
    for j in range(2*num):
        if j==i:
            t = 0
        else:
            t_1 = (3*(xpos[i] - xpos[j])*(zpos[i] - zpos[j]))/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
            t_2 = (3*(xpos[j] - xpos[i])*(zpos[j] - zpos[i]))/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
            t = t_1 + t_2
        S = S + t
    K_xz[i,i] = (Q**2)/(8*pi*epsilon0)*S

for i in range(2*num):
    S=0
    for j in range(i+1, 2*num):
        t_1 = (3*(xpos[i] - xpos[j])*(zpos[i] - zpos[j]))/(L(xpos[i], xpos[j], zpos[i], zpos[j]))**5
        t_2 = (3*(xpos[j] - xpos[i])*(zpos[j] - zpos[i]))/(L(xpos[j], xpos[i], zpos[j], zpos[i]))**5
        S = t_1 + t_2
        K_xz[i,j] = (Q**2)/(8 * pi * epsilon0) *S
        K_xz[j,i] = K_xz[i,j]
        
K_zx = K_xz.T

print('K_xz = \n', K_xz)
print('K_zx = \n', K_zx)

K = np.block([[K_xx,K_xy,K_xz],[K_yx,K_yy,K_yz], [K_zx, K_zy, K_zz]])

W = (1/m0) * K
print('K = \n', W)

# eigenvalue
Enf = sorted((LA.eig(W)[0] ** (1/2))/(2*pi))
normalized_Enf = sorted((LA.eig(W)[0] ** (1/2))/(2*pi * omega_z))

f = np.zeros(4).reshape(4)
for i in range(4):
    f[i] = normalized_Enf[i]

print('eigenvalue = \n',Enf)
print('eigenvalue',normalized_Enf)

#draw graph
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111)

ax.set_xlabel(r'$d$ ($\mu$m)',fontsize = FontSize)
ax.set_ylabel(r'$\omega_z/\omega_z0$',fontsize = FontSize)
ax.spines['top'].set_linewidth(spines)
ax.spines['bottom'].set_linewidth(spines)
ax.spines['left'].set_linewidth(spines)
ax.spines['right'].set_linewidth(spines)
ax.tick_params(direction = 'in'
               ,length = length
               ,width = width
               ,labelsize = LabelSize
               )
dis = x0*2
d = np.full(4, dis * 10 ** 6)
ax.scatter(d,f,color='black')
plt.show()


