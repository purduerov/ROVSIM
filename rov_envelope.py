import numpy as np
import time
from random import random
from thrust_mapper import ThrustMapper
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

THRUST_MAX = 3.71
THRUST_MIN = -2.9
POINTS = 100000

tm = ThrustMapper()


def mag(vec):  # gets the magnitude of the vector really fast for small size vectors
    return vec.dot(vec)**0.5


def fibonacci_sphere(samples=1):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append((x, y, z))

    return points


start = time.time()

print('Generating points')
sphere_points = fibonacci_sphere(POINTS)
envelope_points = []
fine_points = []

print("Projecting")
for point_num, point in enumerate(sphere_points):
    thrust_vec = np.array(point)

    tm.fine = False
    # desired_thrust_final = np.array([*thrust_vec, 0, 0, 0], dtype=np.float)
    # desired_thrust_final = np.array([0, 0, 0, *thrust_vec], dtype=np.float)

    # thrust_output = tm.thruster_output(desired_thrust_final)
    # fine_points.append(thrust_vec * sum(-np.abs(np.array(thrust_output) / (3 * THRUST_MAX))))

    desired_force_final = np.array([*thrust_vec, 0, 0, 0], dtype=np.float)
    desired_torque_final = np.array([0, 0, 0, 0, 0, 0], dtype=np.float)
    thrust_output = tm.thruster_output2(desired_force_final, desired_torque_final, tm.pseudo_inverse)
    # thrust_output = tm.thruster_output(desired_force_final, desired_torque_final)

    thrust = np.array([0, 0, 0])
    for thruster in thrust_output:
        thrust = thrust + np.matmul(thrust_output, tm.thruster_output_map)[:3]

    envelope_points.append(thrust)

    if not point_num % (POINTS // 10):
        print(point_num)

print(f"Runtime: {time.time() - start}")

print(f"Inscribed Circle Radius: {min([mag(point) for point in envelope_points])}")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# sx, sy, sz = zip(*[list(point) for point in sphere_points])
px, py, pz = zip(*[list(point) for point in envelope_points])
# px1, py1, pz1 = zip(*[list(point) for point in fine_points])

# ax.scatter(sx, sy, sz, color=["red"])
ax.scatter(px, py, pz, color=["blue"])
# ax.scatter(px1, py1, pz1, color=["green"])

ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
ax.set_zlim(-100, 100)

plt.show(block=True)
