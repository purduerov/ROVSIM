import time
import numpy as np
import pybullet as p
from simulation import Simulation
import pygame

from thrust_mapper import ThrustMapper, mag


rate = 1/240.0
print(f"Ticks per Second: {1/rate}")

sim = Simulation(robot_path='/robots/rov/robot.urdf', dt=rate)
tm = ThrustMapper()


desired_thrust_final = np.array([0, 0, 0, 0, 0, 0], dtype=np.float)  # X Y Z R P Y
# direction_scale = np.array([3.71 / 0.266,  # this is the limits of the thrust envelope per direction
#                             3.71 / 0.731,
#                             3.71 / 0.250,
#                             3.71 / 1.733,
#                             3.71 / 2.218,
#                             3.71 / 1.043], dtype=np.float)
# direction_scale *= 0.50  # Split between Translation & Rotation

direction_scale = 1.0
# direction_scale = 3.71 / 4.229  # this is the minimum inscribed sphere
print(f"Direction Scale: {direction_scale}")

debug = False

object_vol = {  # m^3
    'ORIGIN': [0.0, [0, 0, 0]],
    'POWERBOX': [0.004324681831, [0, 0, 0]],
    'FRAME': [0.0008008324, [0, 0, 0]],
    'LOGIC': [0.000100863103, [0, 0, 0]],
    'AIRBOX': [0.000545305774, [0, 0, 0]],
    'THRUSTER': [0.0002243406, [0, 0, 0]],
    'FOOT': [0.0000211825744, [0, 0, 0]],
    'FOAM': [0.00370555386/2.55, [0, 0, 0]],  # Magic foam fudge factor
}


def cob_position():
    """Returns center of buoyancy of the robot
    Returns:
        pos -- (x, y, z) robot center of buoyancy
    """

    net_force = 0
    cob = np.array([0.0, 0.0, 0.0])

    for joint_num in range(p.getNumJoints(sim.robot)):
        res = p.getLinkState(sim.robot, joint_num)
        cob_pos = res[0]
        j_name = str(p.getJointInfo(sim.robot, joint_num)[1].decode()).split('_')[0]

        force = object_vol[j_name][0] * 1000 * 9.81
        cob += (np.array(cob_pos) + np.array(object_vol[j_name][1])) * force
        net_force += force

    return cob / net_force


joint_names = [str(p.getJointInfo(sim.robot, joint)[1]) for joint in range(p.getNumJoints(sim.robot))]
print('\n'.join(joint_names))
thrusters = [joint for joint in range(p.getNumJoints(sim.robot)) if "THRUSTER" in str(p.getJointInfo(sim.robot, joint)[1])]
print(thrusters)

line_ids = {}
for thruster in thrusters:
    line_id = p.addUserDebugLine(lineFromXYZ=[0, 0, 0], lineToXYZ=[0, 0, 0], lineColorRGB=[1, 0, 0], lineWidth=3, lifeTime=0)
    line_ids[thruster] = line_id
print(line_ids)

vel_line = p.addUserDebugLine(lineFromXYZ=[0, 0, 0], lineToXYZ=[0, 0, 0], lineColorRGB=[1, 0, 0], lineWidth=3, lifeTime=0)
com_line = p.addUserDebugLine(lineFromXYZ=[0, 0, 0], lineToXYZ=[0, 0, 0], lineColorRGB=[1, 0, 0], lineWidth=3, lifeTime=0)
bou_line = p.addUserDebugLine(lineFromXYZ=[0, 0, 0], lineToXYZ=[0, 0, 0], lineColorRGB=[1, 0, 0], lineWidth=3, lifeTime=0)


pygame.init()
joystick = None

# This probably only works for windows, so update as necessary
try:
    for i in range(pygame.joystick.get_count()):
        if pygame.joystick.Joystick(i).get_name() == "Xbox 360 Controller":
            joystick = pygame.joystick.Joystick(i)
            joystick.init()
            print(f"Connected to: {joystick.get_name()}")
            break
except pygame.error:
    print("No Gamepad Found")


joy_map_a = ["LX", "LY", "RX", "RY", "LT", "RT"]
joy_map_b = ["X", "C", "S", "T", "LT", "RT", "SE", "ST", "LS", "RS"]

old_axes = None
old_buttons = None

if joystick:
    old_axes = [joystick.get_axis(num) for num in range(joystick.get_numaxes())]
    old_buttons = [joystick.get_button(num) for num in range(joystick.get_numbuttons())]

key_map = {
    65297: np.array([1, 0, 0, 0, 0, 0], dtype=np.float64),  # arrow up +X
    65298: np.array([-1, 0, 0, 0, 0, 0], dtype=np.float64),  # arrow down -X
    65295: np.array([0, 1, 0, 0, 0, 0], dtype=np.float64),  # arrow right +Y
    65296: np.array([0, -1, 0, 0, 0, 0], dtype=np.float64),  # arrow left -Y
    32:    np.array([0, 0, 1, 0, 0, 0], dtype=np.float64),  # space bar +Z
    65306: np.array([0, 0, -1, 0, 0, 0], dtype=np.float64),  # left shift -Z

    107:   np.array([0, 0, 0, 1, 0, 0], dtype=np.float64),  # K +R
    104:   np.array([0, 0, 0, -1, 0, 0], dtype=np.float64),  # H -R
    106:   np.array([0, 0, 0, 0, 1, 0], dtype=np.float64),  # J +P
    110:   np.array([0, 0, 0, 0, -1, 0], dtype=np.float64),  # N -P
    98:   np.array([0, 0, 0, 0, 0, 1], dtype=np.float64),  # B +Y
    109:   np.array([0, 0, 0, 0, 0, -1], dtype=np.float64),  # M -Y
}

while True:
    if joystick:
        pygame.event.pump()

        axes = [joystick.get_axis(num) for num in range(joystick.get_numaxes())]
        buttons = [joystick.get_button(num) for num in range(joystick.get_numbuttons())]

        desired_thrust_final[0] = -axes[joy_map_a.index("LY")]
        desired_thrust_final[1] = -axes[joy_map_a.index("LX")]
        desired_thrust_final[2] = buttons[joy_map_b.index("RT")] - buttons[joy_map_b.index("LT")]

        desired_thrust_final[3] = max(axes[joy_map_a.index("RT")], 0) - max(axes[joy_map_b.index("LT")], 0)
        desired_thrust_final[4] = -axes[joy_map_a.index("RY")]
        desired_thrust_final[5] = -axes[joy_map_a.index("RX")]

        if buttons[joy_map_b.index("ST")] and not old_buttons[joy_map_b.index("ST")]:
            p.resetBasePositionAndOrientation(sim.robot, (0, 0, 0.75), p.getQuaternionFromEuler([0, 0, 0]))
            print("Reset Position")

        if buttons[joy_map_b.index("C")] and not old_buttons[joy_map_b.index("C")]:
            debug ^= True
            print(f"Debug Mode: {debug}")

        if buttons[joy_map_b.index("X")] and not old_buttons[joy_map_b.index("X")]:
            tm.fine ^= True
            print(f"Fine Control: {tm.fine}")

        old_axes = axes
        old_buttons = buttons

    events: dict = p.getKeyboardEvents()

    if events != {}:
        for event in events.items():
            if event[1] == 3:
                print(event)
                if not joystick and event[0] in key_map.keys():
                    print("FUCK")
                    desired_thrust_final += key_map[event[0]]

                if event[0] == 122:  # Z Debug
                    debug ^= True
                    print(f"Debug Mode: {debug}")
                if event[0] == 120:  # X Reset Vel
                    desired_thrust_final = np.array([0, 0, 0, 0, 0, 0], dtype=np.float)
                    print("Reset Velocity")
                if event[0] == 99:  # C Fine
                    tm.fine ^= True
                    print(f"Fine Control: {tm.fine}")
                if event[0] == 97:  # A Reset Pos
                    p.resetBasePositionAndOrientation(sim.robot, (0, 0, 0.75), p.getQuaternionFromEuler([0, 0, 0]))
                    print("Reset Position")

            if event[1] == 4:
                if not joystick and event[0] in key_map.keys():
                    desired_thrust_final -= key_map[event[0]]

    scaled_thrust_final = np.copy(desired_thrust_final)

    if mag(desired_thrust_final) > 1.0:
        scaled_thrust_final /= mag(scaled_thrust_final)

    desired_force = np.copy(scaled_thrust_final)
    desired_force[3:] = 0

    desired_torque = np.copy(scaled_thrust_final)
    desired_torque[:3] = 0

    thrust_output = tm.thruster_output(desired_force, desired_torque)
    # if debug:
    #     print(f"Thruster Outputs: {np.around(thrust_output, 4)}")

    robot_pos, robot_ori = p.getBasePositionAndOrientation(sim.robot)
    t_vel, r_vel = p.getBaseVelocity(sim.robot)
    t_vel, r_vel = np.array(t_vel), np.array(r_vel)

    t_drag = -100.0 * np.sqrt(t_vel.dot(t_vel)) * t_vel
    r_drag = -3.0 * np.array(r_vel)

    if debug:
        p.addUserDebugLine(robot_pos, robot_pos + t_vel * 2, [1, 0, 0], lifeTime=0, lineWidth=3, replaceItemUniqueId=vel_line)

    p.applyExternalTorque(sim.robot, -1, r_drag, flags=p.LINK_FRAME)
    p.applyExternalForce(sim.robot, -1, t_drag, robot_pos, flags=p.WORLD_FRAME)

    if debug:
        p.addUserDebugLine(sim.com_position(), sim.com_position() + np.array([0, 0, 0.2]), [1, 0, 0], lifeTime=0, lineWidth=2, replaceItemUniqueId=com_line)
        p.addUserDebugLine(cob_position(), cob_position() + np.array([0, 0, 0.2]), [0, 1, 0], lifeTime=0, lineWidth=2, replaceItemUniqueId=bou_line)
        # print(f"Distance Between COM and COB: {sim.com_position() - cob_position()}")

    for thrust_num, joint in enumerate(thrusters):
        state = p.getLinkState(sim.robot, joint, computeLinkVelocity=0)
        joint_pos = np.array(state[0])

        thrust_axis = np.matmul(np.array(p.getMatrixFromQuaternion(state[1])).reshape(3, 3), np.array([0, 0, 1]))
        thrust_vec = thrust_axis / mag(thrust_axis) * thrust_output[thrust_num]

        if mag(thrust_vec) > (3.71 * 1.01):
            print(f"THRUSTER: {thrust_num} | {mag(thrust_vec)}")
            thrust_vec *= 3.71 / mag(thrust_vec)

        if debug:
            p.addUserDebugLine(joint_pos, joint_pos + thrust_vec / 4, [0, 1, 0], lifeTime=0, lineWidth=3, replaceItemUniqueId=line_ids[joint])  # This is the thrust vectors

        p.applyExternalForce(sim.robot, joint, thrust_vec, joint_pos, flags=p.WORLD_FRAME)

    for joint in range(p.getNumJoints(sim.robot)):
        name = str(p.getJointInfo(sim.robot, joint)[1].decode())
        name = name.split('_')[0]

        pos = np.array(p.getLinkState(sim.robot, joint, computeLinkVelocity=0)[4])

        if pos[2] < 5:
            if name in object_vol.keys():
                volume = object_vol[name][0]
                buoyant_force = np.array([0, 0, 1]) * volume * 1000 * 9.81

                p.applyExternalForce(sim.robot, joint, buoyant_force, pos + np.array(object_vol[name][1]), flags=p.WORLD_FRAME)

    p.resetDebugVisualizerCamera(0.9, 30, -30, robot_pos)
    p.stepSimulation(sim.physics_client)
    time.sleep(rate)
