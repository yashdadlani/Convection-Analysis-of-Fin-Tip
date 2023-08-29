import numpy as np
import matplotlib.pyplot as plt
from math import *
K = float(input("Enter a value of thermal conductivity : "))

Tb = 250 ; Ts = 25 ; h = 25 ; ht = 10
l = 0.1 ; w = 0.5 ; t = 0.002 ; Ac = w*t ; p = 2*(w+t); 
m = sqrt(h*p/(K*Ac))
hnd = ht/(m*K)

dx = float(input("Enter a value of step size : "))         
n = int(l/dx)
x = np.linspace(0,l,n+1) 
theta = np.empty(n+1,float) 
theta[0] = 1
a = -2 - ((m*dx)**2) 
b = 3 + (ht/K)*(2*dx) 
alpha = [0]*n
alpha[0] = 0
alpha[n-1] = -4-(1*a)
for i in range(1,n-1):
     alpha[i]= 1
beta = [0]*n
beta[n-1] = b-1
for i in range(0,n-1):
     beta[i]= a
delta = [1]*n
delta[n-1]=0
sol = [0]*n
sol[0] = -1
e = [0]*n
f = [0]*n
th = [0]*n
e[0] = delta[0]/beta[0]
f[0] = sol[0]/beta[0]
for i in range(1, n-1):
     d = (beta[i] - (alpha[i]*e[i-1]))
     e[i] = delta[i]/d
     f[i] = (sol[i] - alpha[i]*f[i-1])/d
th[n-1] = (sol[n-1] - (alpha[n-1]*f[n-2])) / (beta[n-1] - (alpha[n-1]*e[n-2]))
for j in range(n-2,-1,-1):
     th[j] = f[j] - e[j]*th[j+1]
 # print('th = ' , [round(i, 5) for i in th])
for i in range(1,n+1):
     theta[i] = th[i-1]
T = np.empty(n+1,float)
print("x\t\t T(numerical temperature)")
for i in range(0,n+1):
    T[i]= theta[i]*(Tb-Ts) + Ts
    print("%f \t %f" %(x[i],T[i]))

q_fof = -1*K*Ac*((T[1]-T[0])/(dx))
q_sof = -1*K*Ac*((4*T[1]-3*T[0]-T[2])/(2*dx))
if K==237:
 print("Numerical Heat Loss for Al from fof scheme : " + str(q_fof))
 print("Numerical Heat Loss for Al from sof scheme : " + str(q_sof))
 plt.plot(x,T,marker='o')
 plt.xlabel('Location in m')
 plt.ylabel('Dimensional temperature')
 plt.grid()
 plt.title("Temp. Distribution for Aluminum")
 plt.show()
if K==17:
 print("Numerical Heat Loss for SS from fof scheme : " + str(q_fof))
 print("Numerical Heat Loss for SS from sof scheme : " + str(q_sof))
 plt.plot(x,T,marker='o')
 plt.xlabel('Location in m')
 plt.ylabel('Dimensional temperature')
 plt.grid()
 plt.title("Temp. Distribution for SS")
 plt.show()
if K==0.8:
 print("Numerical Heat Loss for Glass from fof scheme : " + str(q_fof))
 print("Numerical Heat Loss for Glass from sof scheme : " + str(q_sof))
 plt.plot(x,T,marker='o')
 plt.xlabel('Location in m')
 plt.ylabel('Dimensional temperature')
 plt.grid()
 plt.title("Temp. Distribution for Glass")
 plt.show()



