#!/usr/bin/python3

""" This program describes the motion of a droplet inside and outside a cloud until it reaches the ground, using a simplified approach. """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# -----------------------------------------------------------------------------------------------------
u = np.zeros(12000000)
r = np.zeros(12000000)
t = np.zeros(12000000)
z = np.zeros(12000000)

# constants

E_coalescence = 0.7
# E_collection = E_collision * E_coalescence

# collision efficiency
r_coll = [20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500, 600, 1000, 1400, 1800, 2400, 3000]
E_coll = [0.17, 0.37, 0.55, 0.58, 0.68, 0.76, 0.81, 0.83, 0.87, 0.87, 0.88, 0.88, 0.88, 0.88, 0.88, 0.86, 0.83, 0.81]

# data

M = 6*10**(-3)        # kg/m^3
rl = 1000             # kg/m^3
w = 5                 # m/s
k1 = 1.19*10**8       # m^(-1)*sec^(-1)
k2 = 8*10**3          # sec^(-1)
k3 = 2.01*10**2       # m^(1/2)*sec^(-1)
k4 = 2.11*10**2       # m^(1/2)*sec^(-1)
dz = 0.001            # m
zo = 250              # m
ro = 15*10**(-6)      # m
z_base = 2000         # m
depth = 6000          # m

u[0] = 0              # m/s
t[0] = 0              # s
z[0] = zo + z_base    # m
r[0] = ro             # m

# temperature at starting point
T = 10-0.005*250

i = 0
j = 0

####################################################################
# ascent

while u[i] <= w:

    i = i + 1

    # velocity
    if r[i-1] < 40*10**(-6):
        u[i] = k1*r[i-1]**2
    elif r[i-1] >= 40*10**(-6) and r[i-1] <= 600*10**(-6):
        u[i] = k2*r[i-1]
    elif r[i-1] > 600*10**(-6) and r[i-1] < 2*10**(-3):
        u[i] = k3*r[i-1]**(0.5)
    else:
        u[i] = k4*r[i-1]**(0.5)

    if u[i] >= w:
        i = i - 1
        break

    # collection efficiency
    if r[i-1] <= r_coll[j]*10**(-6):
        E = E_coll[j]*E_coalescence
        if j!=17:
            j = j + 1
           
    # radius, height, time   
    z[i] = z[i-1] + dz
    dr = E*M/(4*rl) * u[i]/(w-u[i]) * dz
    r[i] = r[i-1] + dr
    t[i] = t[i-1] + dz/(w-u[i])


t_asc = t[i]/60
r_max = r[i]*1000

if z[i]>z_base+depth :
    print('The droplet has exceeded the cloud!')
    print(z[i])
   
else:
    print('-----------------------------------------------------')
    print('Maximum height from starting point = ',z[i]-z[0], 'm')
    print('Maximum height from ground = ',z[i], 'm')
    print('Time of ascent = ',t_asc,'min')
    print('Size of droplet at maximum height = ',r_max, 'mm')
    print('-----------------------------------------------------')

####################################################################
# descent

while z[i] >= z_base :

    i = i + 1
    # height
    z[i] = z[i-1] - dz
    
    if z[i] < z_base :
        i = i - 1
        break

    # velocity
    if r[i-1] < 40*10**(-6):
        u[i] = k1*r[i-1]**2
    elif r[i-1] >= 40*10**(-6) and r[i-1] <= 600*10**(-6):
        u[i] = k2*r[i-1]
    elif r[i-1] > 600*10**(-6) and r[i-1] < 2*10**(-3):
        u[i] = k3*r[i-1]**(1/2)
    else:   
        u[i] = k4*r[i-1]**(1/2)

    # radius, time
    dr = E*M/(4*rl) * dz * u[i]/(u[i]-w)
    r[i] = r[i-1] + dr
    t[i] = t[i-1] + dz/(u[i]-w)
   
index_cloud = i
t_cloud = t[i]/60
t_desc = t_cloud-t_asc
r_min = r[i]*1000    

print('Time of descent = ', t_desc, 'min')
print('Size of droplet at cloud base = ', r_min, 'mm')
print('-----------------------------------------------------')

print('Total time inside cloud = ', t_cloud, 'min')
print('-----------------------------------------------------')    

t_fall = (z_base/u[i])/60
t_total = t_fall + t_cloud

print('Time of fall = ', t_fall, 'min')
print('Total time of rain process = ', t_total, 'min')
print('-----------------------------------------------------')    

####################################################################
# out of cloud

xi1 = 10**2 * 10**(-12)    # m^2/s
c = (0.7-1)*xi1             # m^2/s

ro = r[i]
dt = 0.001

while z[i] >= 0:

    i = i + 1
    # velocity, height, time
    u[i] = u[i-1]
    dz = u[i]*dt
    z[i] = z[i-1] - dz
    t[i] = t[i-1] + dt
    
    if z[i] < 0:
        i = i - 1
        break

    # radius
    r[i] = (r[i-1]**2+2*c*dt)**(1/2)


r_ground = r[i]*1000
diff = (r[i]-ro)*10**(6)

print('Size of droplet at ground = ', r_ground, 'mm')
print('Difference between the raindrop at cloud base and raindrop at ground = ', diff, 'microns')

########################################################
# plot 1

height = z[0:i+1]/1000
radius = r[0:i+1]*1000
time = t[0:i+1]/60
velocity = u[0:i+1]


height_cloud = z[0:index_cloud+1]/1000
time_cloud = t[0:index_cloud+1]/60
radius_cloud = r[0:index_cloud+1]*1000


fig1, ax1 = plt.subplots()
ax1.tick_params(axis="y",direction="in")
ax1.tick_params(axis="x",direction="in")
plt.plot(radius,height)
plt.xlabel('r (mm)', labelpad = 10)
plt.ylabel('H (km)', labelpad = 10)
plt.xlim(-0.1,7.,1)
plt.ylim(0,8,1)
plt.title('Height with radius of droplet', y = 1.02)            
plt.savefig('height-radius.jpg')

# plot 2
plt.clf()

fig2, ax2 = plt.subplots()
ax2.tick_params(axis="y",direction="in")
ax2.tick_params(axis="x",direction="in")
plt.plot(time, height, color='red')
plt.xlabel('t (min)', labelpad = 10)
plt.ylabel('H (km)' , labelpad = 10)
plt.ylim(0,8,0.5)
plt.xlim(0,30,2)
plt.title('Height with time', y = 1.02)            
plt.savefig('height-time.jpg')

# plot 3
plt.clf()

fig3, ax3 = plt.subplots()
ax3.tick_params(axis="y",direction="in")
ax3.tick_params(axis="x",direction="in")
plt.plot(radius_cloud, height_cloud, color='green')
plt.xlabel('R (mm)', labelpad = 10)
plt.ylabel('H (km)' , labelpad = 10)
plt.ylim(2,8,0.5)
plt.xlim(-0.1,7.,1)
plt.title('Height with radius inside the cloud', y = 1.02)        
plt.savefig('height-radius.jpg')

# plot 4
plt.clf()

fig4, ax4 = plt.subplots()
ax4.tick_params(axis="y",direction="in")
ax4.tick_params(axis="x",direction="in")
plt.plot(velocity, height, color='green')
plt.xlabel('u (m/s)', labelpad = 10)
plt.ylabel('H (km)' , labelpad = 10)
plt.ylim(0,8.,1)
plt.xlim(0,18,1)
plt.title('Height with velocity of droplet', y = 1.02)        
plt.savefig('height-velocity.jpg')

# plot 5
plt.clf()

fig5, ax5 = plt.subplots()
ax5.tick_params(axis="y",direction="in")
ax5.tick_params(axis="x",direction="in")
plt.plot( radius, velocity, color='orange')
plt.ylabel('u (m/s)', labelpad = 10)
plt.xlabel('r (mm)' , labelpad = 10)
plt.ylim(0,18.,1)
plt.xlim(0,7.,1)
plt.title('Velocity with radius of droplet', y = 1.02)        
plt.savefig('velocity-radius.jpg')
