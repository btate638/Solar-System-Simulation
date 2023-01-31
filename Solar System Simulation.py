# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 11:43:37 2022

@author: btate
"""


import numpy as np
import matplotlib.pyplot as plt
import timeit

start = timeit.default_timer()

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150

N = 11
G = 6.67408e-11

dt = 50000
t_max = 86400 * 365.25 * 1

names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Asteroid']
r = [0, 5.79e10, 1.082e11, 1.495978707e11, 2.280e11, 7.78340821e11, 1.4320e12, 2.8670e12, 4.5150e12, 5.9064e12, 4.9e11]
m = [1.98892e30, 3.30e23, 4.87e24, 5.9742e24, 6.42e23, 1.8986e27, 5.68e26, 8.68e25, 1.02e26, 1.30e22, 0.95e21]
x = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
y = [0, 5.79e10, 1.082e11, 1.495978707e11, 2.280e11, 7.78340821e11, 1.4320e12, 2.8670e12, 4.5150e12, 5.9064e12, 4.9e11]
vx = [0, 
      ((G*m[0])/r[1])**0.5, 
      ((G*m[0])/r[2])**0.5, 
      ((G*m[0])/r[3])**0.5, 
      ((G*m[0])/r[4])**0.5, 
      ((G*m[0])/r[5])**0.5, 
      ((G*m[0])/r[6])**0.5, 
      ((G*m[0])/r[7])**0.5, 
      ((G*m[0])/r[8])**0.5, 
      ((G*m[0])/r[9])**0.5, 
      ((G*m[0])/r[10])**0.5]
vy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

ax = np.zeros(N)
ay = np.zeros(N)
bx = np.zeros(N)
by = np.zeros(N)
cx = np.zeros(N)
cy = np.zeros(N)
dx = np.zeros(N)
dy = np.zeros(N)

avx = np.zeros(N)
avy = np.zeros(N)
bvx = np.zeros(N)
bvy = np.zeros(N)
cvx = np.zeros(N)
cvy = np.zeros(N)
dvx = np.zeros(N)
dvy = np.zeros(N)

#List of displacement values for each planet 
x_array = [[],[],[],[],[],[],[],[],[],[],[]]
y_array = [[],[],[],[],[],[],[],[],[],[],[]]
x_array_au = [[],[],[],[],[],[],[],[],[],[],[]]
y_array_au = [[],[],[],[],[],[],[],[],[],[],[]]

#Kinetic Energy 
KE = [[],[],[],[],[],[],[],[],[],[],[]]

#Gravitational Potential Energy
GPE = [[],[],[],[],[],[],[],[],[],[],[]]

#Total Energy
TE = [[],[],[],[],[],[],[],[],[],[],[]]

#Total Initial Linear Momentum of Planets
initial_mom = 0
for i in range(1, N):
        initial_mom += m[i] * (vx[i]**2 + vy[i]**2)**0.5 

#print(initial_mom)
vx[0] = initial_mom / m[0]
print(vx[0])

def accnx(mj, xj, yj, xi, yi):
    xdiff = xj - xi
    ydiff = yj - yi
    r = (xdiff**2 +ydiff**2)**0.5
    return ((G*mj*(xdiff)) / (r**3))
    
def accny(mj, xj, yj, xi, yi):
    xdiff = xj-xi
    ydiff = yj-yi
    r = (xdiff**2 +ydiff**2)**0.5
    return ((G*mj*(ydiff)) / (r**3))


for t in np.arange(dt, t_max, dt):
    #RK 1
    for i in range(0, N):
        ax[i] = vx[i]
        ay[i] = vy[i]
        avx[i] = 0
        avy[i] = 0
       
        for j in range(0, N):
            if i != j:  
                avx[i] += accnx(m[j], x[j], y[j], x[i], y[i])
                avy[i] += accny(m[j], x[j], y[j], x[i], y[i])
    #RK2        
    for i in range(0, N):
    
        bx[i] = vx[i] + (dt/2)*(avx[i]) 
        by[i] = vy[i] + (dt/2)*(avy[i]) 
        bvx[i] = 0
        bvy[i] = 0
    
        for j in range(0, N):
            if i != j:
                bvx[i] += accnx(m[j], x[j]+(dt*ax[j]/2), y[j]+(dt*ay[j]/2), x[i]+(dt*ax[i]/2), y[i]+(dt*ay[i]/2))
                bvy[i] += accny(m[j], x[j]+(dt*ax[j]/2), y[j]+(dt*ay[j]/2), x[i]+(dt*ax[i]/2), y[i]+(dt*ay[i]/2))
    #RK3
    for i in range(0, N):
        
        cx[i] = vx[i] + (dt/2)*(bvx[i]) 
        cy[i] = vy[i] + (dt/2)*(bvy[i]) 
        cvx[i] = 0
        cvy[i] = 0
        
        for j in range(0, N):
            if i != j:
                cvx[i] += accnx(m[j], x[j]+(dt*bx[j]/2), y[j]+(dt*by[j]/2), x[i]+(dt*bx[i]/2), y[i]+(dt*by[i]/2))
                cvy[i] += accny(m[j], x[j]+(dt*bx[j]/2), y[j]+(dt*by[j]/2), x[i]+(dt*bx[i]/2), y[i]+(dt*by[i]/2))
    #RK4            
    for i in range(0, N):
        
        dx[i] = vx[i] + (dt)*(cvx[i]) 
        dy[i] = vy[i] + (dt)*(cvy[i]) 
        dvx[i] = 0
        dvy[i] = 0
        
        for j in range(0, N):
            if i != j:
                dvx[i] += accnx(m[j], x[j]+(dt*cx[j]), y[j]+(dt*cy[j]), x[i]+(dt*cx[i]), y[i]+(dt*cy[i]))
                dvy[i] += accny(m[j], x[j]+(dt*cx[j]), y[j]+(dt*cy[j]), x[i]+(dt*cx[i]), y[i]+(dt*cy[i]))
                
    #if x[3] < 0 and x[3] + dt/6 * (ax[3] + 2*bx[3] + 2*cx[3] + dx[3]) > 0 :
        #print(t/86400)
        
    #RK Final  
    for i in range(0, N):
        x[i] = x[i] + dt/6 * (ax[i] + bx[i]*2 + cx[i]*2 + dx[i])
        y[i] = y[i] + dt/6 * (ay[i] + by[i]*2 + cy[i]*2 + dy[i])
        
        vx[i] = vx[i] + dt/6 * (avx[i] + bvx[i]*2 + cvx[i]*2 + dvx[i])
        vy[i] = vy[i] + dt/6 * (avy[i] + bvy[i]*2 + cvy[i]*2 + dvy[i])
         
        x_array[i].append(x[i])
        y_array[i].append(y[i])
        KE[i].append(0.5 * m[i] * (vx[i]**2 + vy[i]**2))
    
    for i in range(1, N):
        GPE_temp = 0
        for j in range(0, N):
            if i != j:
                
                xdiff = x[j]-x[i]
                ydiff = y[j]-y[i]
                r = (xdiff**2 +ydiff**2)**0.5
                
                GPE_temp += ((G * m[i] * m[j])/ r)
                
        GPE[i].append(GPE_temp)
        TE[i].append(GPE[i][-1] + KE[i][-1])
        
for i in range(0,N): 
    for j in range(len(x_array[i])):
        x_array_au[i].append(x_array[i][j]/1.495978707e11)
        y_array_au[i].append(y_array[i][j]/1.495978707e11)
    
#Plot
plt.plot(x_array_au[1], y_array_au[1],
         x_array_au[2], y_array_au[2],
         x_array_au[3], y_array_au[3],
         x_array_au[4], y_array_au[4],
         x_array_au[5], y_array_au[5],
         x_array_au[6], y_array_au[6],
         x_array_au[7], y_array_au[7],
         x_array_au[8], y_array_au[8],
         x_array_au[9], y_array_au[9],
         x_array_au[10], y_array_au[10]
         )
plt.legend([names[1],
           names[2],
           names[3],
           names[4],
           names[5],
           names[6],
           names[7],
           names[8],
           names[9],
           names[10]
           ],
           loc='upper right')
plt.axis('square')
plt.ylabel('AU')
plt.xlabel('AU')
plt.title("Plot of Solar System using Runge-Kutta method")
plt.show()

stop = timeit.default_timer()

print('Runtime:', '{:.2f}'.format(stop - start), 'seconds')
