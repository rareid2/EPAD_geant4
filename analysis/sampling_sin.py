import numpy as np
import matplotlib.pyplot as plt
import random

import seaborn as sns
sns.set_theme()

samples = []
for i in range(10000):
    u = random.random() 
    sin_dist = np.arccos(u*2 - 1)# - np.deg2rad(90)
    samples.append(np.rad2deg(sin_dist))

plt.hist(samples,bins=30)
plt.show()
plt.close()





# Get x values of the sine wave
time        = np.arange(0, 10, 0.1);


# Amplitude of the sine wave is sine of a variable like time
amplitude   = np.sin(time)

for i in range(10):
    amplitude += amplitude

# Plot a sine wave using time and amplitude obtained for the sine wave
plt.plot(time, amplitude)

# Give a title for the sine wave plot
plt.title('Sine wave')

# Give x axis label for the sine wave plot
plt.xlabel('Time')

 
# Give y axis label for the sine wave plot
plt.ylabel('Amplitude = sin(time)')


plt.grid(True, which='both')
plt.axhline(y=0, color='k')

plt.show()
plt.close()

npoints = 1000
fig = plt.figure()
ax = plt.axes(projection='3d')

phis = []
thetas = []
while len(phis) < npoints:
    u1 = random.random()*2-1 
    u2 = random.random()*2-1
    
    if u1**2 + u2**2 >= 1:
        pass
    else:
        phis.append(u1)
        thetas.append(u2)
phis = np.array(phis)
thetas = np.array(thetas)


r = 5
z = r * 2 * phis * np.sqrt(1- phis**2 - thetas**2)
x = r * 2 * thetas * np.sqrt(1- phis**2 - thetas**2)
y = r * (1 - 2 * (phis**2 + thetas**2))

# okay that works, let's give them a direction now

pitch_angle = [np.arccos(random.random()*2 - 1) for i in range(npoints)]
pitch_angle = np.array(pitch_angle)
#pitch_angles = np.arccos(pitch_angle*2 - 1)# - np.deg2rad(90)

gyros = [random.random()*2*np.pi for i in range(npoints)]
gyros = np.array(gyros)

# inforce inward pointing
ww = np.sin(pitch_angle)*np.cos(gyros)
uu = np.sin(pitch_angle)*np.sin(gyros)
vv = np.cos(pitch_angle)

u = -1*uu
v = -1*vv
w = -1*ww

ax.scatter(x,y,z)
ax.quiver(x, y, z, u, v, w, length=4, normalize=True)

plt.show()
plt.close()

zz = np.array([0,1,0])

# find angle
zz_norm = np.linalg.norm(zz)
# directions
uvw = [(u1,v1,w1) for (u1,v1,w1) in zip(u,v,w)]
uvw = np.array(uvw)

angs = []
for vec in uvw:
    u_norm = np.linalg.norm(vec)
    ang = np.arccos(np.dot(vec,zz) / (u_norm * zz_norm))
    ang = np.rad2deg(ang)
    angs.append(ang)

# okay the real issue is that we dont know how to define this vector
# for each random point on the sphere, the direction should be perpendicular to the z axis

plt.hist(angs,bins=30)
plt.show()
plt.close()

# okay lets work on the first part
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.quiver(0, 0, 0, 0, 0, 1, length=.5, normalize=True)

# vector with an arbitrary angle to z axis
# 45 deg angle
vec = []

test_angle = np.deg2rad(90)
az_angle = 30

xx = -1*np.sin(test_angle)*np.cos(az_angle)
yy = -1*np.sin(test_angle)*np.sin(az_angle)
zz = -1*np.cos(test_angle)

ax.quiver(1, 1, 0, xx, yy, zz, length=.5, normalize=True)


ax.set_zlim([-1,1])
#plt.show()
plt.close()