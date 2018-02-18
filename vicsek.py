"""
==========================
Vicsek Model
==========================
"""

#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import cython
from matplotlib.widgets import Button
import sys

#Define Variables
v=0.03
eta=0.1
L=3.1
N=40

#Initial State

def initialize():
    global x,y,theta,theta_new,order,order_new,order_inst
    x=np.random.uniform(0,L,N)
    y=np.random.uniform(0,L,N)
    theta=np.random.uniform(-np.pi,np.pi,N)
    theta_new=np.random.uniform(-np.pi,np.pi,N)
    order=0
    order_new=1
    order_inst=[]

#Functions
def dist(i,j):
    return np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)

def update():
    for i in range(N):
        theta_near = []
        for j in range(N):
            if dist(i,j)<1:
                theta_near.append(theta[j])
        theta_new[i] = np.mean(theta_near)+eta*np.random.uniform(-np.pi,np.pi)
    for i in range(N):
        theta[i]=theta_new[i]
        x[i]=(x[i]+v*np.cos(theta[i]))%L
        y[i]=(y[i]+v*np.sin(theta[i]))%L

#Run Code

#time_limit=100
"""
order_arr=[]
eta_arr = np.linspace(0,1,20)
for eta in eta_arr:
    initialize()
    while abs(order-order_new)>0.001:
        order=order_new
        update()
        order_new=(np.sqrt(np.sum(np.cos(theta))**2+np.sum(np.sin(theta))**2))/N
        order_inst.append(order_new)
    order_arr.append(np.mean(order_inst))

plt.scatter(eta_arr,order_arr)
plt.show()
"""


#Plot

initialize()

def reset(event):
    sys.exit()

for t in range(100):
    axes = plt.gca()
    axes.set_xlim([0,L])
    axes.set_ylim([0,L])
    update()
    plt.quiver(x,y,np.cos(theta),np.sin(theta), units='width',scale=30,width=0.005)
    #plt.scatter(x, y)
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color='blue', hovercolor='0.975')
    button.on_clicked(reset)
    plt.pause(0.01)
    plt.gcf().clear()
