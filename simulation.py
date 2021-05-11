import math
import time
import numpy as np
import pybullet as p
import os


class Simulation:
    def __init__(self, robot_path, floor=True, fixed=False, transparent=False, gui=True, real_time=True, panels=False, use_urdf_inertia=True, dt=0.002):
        self.dir = os.path.dirname(os.path.abspath(__file__))
        self.gui = gui
        self.real_time = real_time
        self.t = 0
        self.start = time.time()
        self.dt = dt
        self.mass = None
        self.physics_client = None

        # Debug lines drawing
        self.lines = []
        self.current_line = 0
        self.last_lines_draw = 0
        self.line_colors = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]]

        # Instantiating bullet
        if gui:
            self.physics_client = p.connect(p.GUI)
        else:
            self.physics_client = p.connect(p.DIRECT)
        p.setGravity(0, 0, -9.81)

        # Light GUI
        if not panels:
            p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
            p.configureDebugVisualizer(p.COV_ENABLE_SEGMENTATION_MARK_PREVIEW, 0)
            p.configureDebugVisualizer(p.COV_ENABLE_DEPTH_BUFFER_PREVIEW, 0)
            p.configureDebugVisualizer(p.COV_ENABLE_RGB_BUFFER_PREVIEW, 0)

        p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 1)

        # Loading floor and/or plane ground
        if floor:
            self.floor = p.loadURDF(self.dir+'/bullet/plane.urdf')
        else:
            self.floor = None

        # Loading robot
        start_pos = [0, 0, 0]
        if not fixed:
            start_pos[2] = 1

        start_orientation = p.getQuaternionFromEuler([0, 0, 0])
        flags = p.URDF_USE_SELF_COLLISION_EXCLUDE_ALL_PARENTS
        if use_urdf_inertia:
            flags += p.URDF_USE_INERTIA_FROM_FILE

        self.robot = p.loadURDF(robot_path, start_pos, start_orientation, flags=flags, useFixedBase=fixed)

        # Engine parameters
        p.setPhysicsEngineParameter(fixedTimeStep=self.dt, maxNumCmdPer1ms=0)
        # p.setRealTimeSimulation(1)
        # p.setPhysicsEngineParameter(numSubSteps=1)

        # Retrieving joints and frames
        self.joints = {}
        self.joints_info = {}
        self.joints_indexes = {}
        self.frames = {}
        self.maxTorques = {}

        # Collecting the available joints
        n = 0
        for i in range(p.getNumJoints(self.robot)):
            joint_info = p.getJointInfo(self.robot, i)
            name = joint_info[1].decode('utf-8')
            if not name.endswith('_fixing') and not name.endswith('_passive'):
                if '_frame' in name:
                    self.frames[name] = i
                else:
                    self.joints_indexes[name] = n
                    n += 1
                    self.joints[name] = i
                    self.joints_info[name] = {
                        'type': joint_info[2]
                    }
                    if joint_info[8] < joint_info[9]:
                        self.joints_info[name]['lowerLimit'] = joint_info[8]
                        self.joints_info[name]['upperLimit'] = joint_info[9]

        # Changing robot opacity if transparent set to true
        if transparent:
            for i in range(p.getNumJoints(self.robot)):
                p.changeVisualShape(self.robot, i, rgbaColor=[0.3, 0.3, 0.3, 0.3])

        print('* Found '+str(len(self.joints))+' DOFs')
        print('* Found '+str(len(self.frames))+' frames')

    def camera_look_at(self, target):
        """Control the look of the visualizer camera
        Arguments:
            target {tuple} -- target as (x,y,z) tuple
        """
        if self.gui:
            params = p.getDebugVisualizerCamera()
            p.resetDebugVisualizerCamera(
                params[10], params[8], params[9], target)

    def get_robot_pose(self):
        """Gets the robot (origin) position
        Returns:
            (tuple(3), tuple(3)) -- (x,y,z), (roll, pitch, yaw)
        """
        pose = p.getBasePositionAndOrientation(self.robot)
        return pose[0], p.getEulerFromQuaternion(pose[1])

    def frame_to_world_matrix(self, frame):
        """Gets the given frame to world matrix transformation. can be a frame name
        from URDF/SDF or "origin" for the part origin
        Arguments:
            frame {str} -- frame name
        Returns:
            np.matrix -- a 4x4 matrix
        """

        if frame == 'origin':
            frame_in_world_frame = p.getBasePositionAndOrientation(self.robot)
        else:
            frame_in_world_frame = p.getLinkState(self.robot, self.frames[frame])

        return self.pose_to_matrix(frame_in_world_frame)

    def transformation(self, frame_a, frame_b):
        """Transformation matrix AtoB
        Arguments:
            frameA {str} -- frame A name
            frameB {str} -- frame B name
        Returns:
            np.matrix -- A 4x4 matrix
        """
        a_in_world = self.frame_to_world_matrix(frame_a)
        b_in_world = self.frame_to_world_matrix(frame_b)

        return np.linalg.inv(b_in_world) * a_in_world

    def set_robot_pose(self, pos, orn):
        """Sets the robot (origin) pose
        Arguments:
            pos {tuple} -- (x,y,z) position
            orn {tuple} -- (x,y,z,w) quaternions
        """
        p.resetBasePositionAndOrientation(self.robot, pos, orn)

    def reset(self, height=1.0, orientation='straight'):
        """Resets the robot for experiment (joints, robot position, simulator time)
        Keyword Arguments:
            height {float} -- height of the reset (m) (default: {0.55})
            orientation {str} -- orientation (straight, front or back) of the robot (default: {'straight'})
        """
        self.lines = []
        self.t = 0
        self.start = time.time()

        # Resets the robot position
        orn = [0, 0, 0]
        if orientation == 'front':
            orn = [0, math.pi/2, 0]
        elif orientation == 'back':
            orn = [0, -math.pi/2, 0]
        self.set_robot_pose([0, 0, height], p.getQuaternionFromEuler(orn))

        # Reset the joints to 0
        for entry in self.joints.values():
            p.resetJointState(self.robot, entry, 0)

    def frame(self, frame):
        """Gets the given frame
        Arguments:
            frame {str} -- frame name
        Returns:
            tuple -- (pos, orn), where pos is (x, y, z) and orn is quaternions (x, y, z, w)
        """
        joint_state = p.getLinkState(self.robot, self.frames[frame])
        return joint_state[0], joint_state[1]

    def frames(self):
        """Gets the available frames in the current robot model
        Returns:
            dict -- dict of str -> (pos, orientation)
        """
        frames = {}

        for name in self.frames.keys():
            joint_state = p.getLinkState(self.robot, self.frames[name])
            pos = joint_state[0]
            orientation = p.getEulerFromQuaternion(joint_state[1])
            frames[name] = [pos, orientation]

        return frames

    def reset_joints(self, joints):
        """Reset all the joints to a given position
        Arguments:
            joints {dict} -- dict of joint name -> angle (float, radian)
        """
        for name in joints:
            p.resetJointState(self.robot, self.joints[name], joints[name])

    def set_joints(self, joints):
        """Set joint targets for motor control in simulation
        Arguments:
            joints {dict} -- dict of joint name -> angle (float, radian)
        Raises:
            Exception: if a joint is not found, exception is raised
        Returns:
            applied {dict} -- dict of joint states (position, velocity, reaction forces, applied torque)
        """
        applied = {}

        for name in joints.keys():
            if name in self.joints:
                if name.endswith('_speed'):
                    p.setJointMotorControl2(
                        self.robot, self.joints[name], p.VELOCITY_CONTROL, targetVelocity=joints[name])
                else:
                    if name in self.maxTorques:
                        max_torque = self.maxTorques[name]
                        p.setJointMotorControl2(self.robot, self.joints[name], p.POSITION_CONTROL, joints[name], force=max_torque)
                    else:
                        p.setJointMotorControl2(self.robot, self.joints[name], p.POSITION_CONTROL, joints[name])

                applied[name] = p.getJointState(self.robot, self.joints[name])
            else:
                raise Exception("Can't find joint %s" % name)

        return applied

    def robot_mass(self):
        """Returns the robot mass
        Returns:
            float -- the robot mass (kg)
        """
        if self.mass is None:
            k = -1
            self.mass = 0
            while True:
                if k == -1 or p.getLinkState(self.robot, k) is not None:
                    d = p.getDynamicsInfo(self.robot, k)
                    self.mass += d[0]
                else:
                    break
                k += 1

        return self.mass

    def com_position(self):
        """Returns center of mass of the robot
        Returns:
            pos -- (x, y, z) robot center of mass
        """

        k = -1
        mass = 0
        com = np.array([0., 0., 0.])
        while True:
            if k == -1:
                pos, _ = p.getBasePositionAndOrientation(self.robot)
            else:
                res = p.getLinkState(self.robot, k)
                if res is None:
                    break
                pos = res[0]

            d = p.getDynamicsInfo(self.robot, k)
            m = d[0]
            com += np.array(pos) * m
            mass += m

            k += 1

        return com / mass

    def add_debug_position(self, position, color=None, duration=30):
        """Adds a debug position to be drawn as a line
        Arguments:
            position {tuple} -- (x,y,z) (m)
        Keyword Arguments:
            color {tuple} -- (r,g,b) (0->1) (default: {None})
            duration {float} -- line duration on screen before disapearing (default: {30})
        """
        if color is None:
            color = self.line_colors[self.current_line % len(self.line_colors)]

        if self.current_line >= len(self.lines):
            self.lines.append({})

        self.lines[self.current_line]['update'] = True
        self.lines[self.current_line]['to'] = position
        self.lines[self.current_line]['color'] = color
        self.lines[self.current_line]['duration'] = duration

        self.current_line += 1

    def draw_debug_lines(self):
        """Updates the drawing of debug lines"""
        self.current_line = 0
        if time.time() - self.last_lines_draw > 0.05:
            for line in self.lines:
                if 'from' in line:
                    if line['update']:
                        p.addUserDebugLine(line['from'], line['to'], line['color'], 2, line['duration'])
                        line['update'] = False
                    else:
                        del line['from']
                line['from'] = line['to']

            self.last_lines_draw = time.time()

    def contact_points(self):
        """Gets all contact points and forces
        Returns:
            list -- list of entries (link_name, position in m, normal force vector, force in N)
        """
        result = []
        contacts = p.getContactPoints(bodyA=self.floor, bodyB=self.robot)
        for contact in contacts:
            link_index = contact[4]
            if link_index >= 0:
                link_name = (p.getJointInfo(self.robot, link_index)[12]).decode()
            else:
                link_name = 'base'
            result.append((link_name, contact[6], contact[7], contact[9]))

        return result

    def auto_collisions(self):
        """Returns the total amount of N in auto-collisions (not with ground)
        Returns:
            float -- Newtons of collisions not with ground
        """
        total = 0
        for k in range(1, p.getNumJoints(self.robot)):
            contacts = p.getContactPoints(bodyA=k)
            for contact in contacts:
                if contact[2] != self.floor:
                    total += contact[9]
        return total

    def execute(self):
        """Executes the simulation infinitely (blocks)"""
        while True:
            self.tick()

    def tick(self):
        """Ticks one step of simulation. If real_time is True, sleeps to compensate real time"""
        self.t += self.dt
        self.draw_debug_lines()

        p.stepSimulation()
        delay = self.t - (time.time() - self.start)
        if delay > 0 and self.real_time:
            time.sleep(delay)

    @staticmethod
    def pose_to_matrix(pose):
        """Converts a pyBullet pose to a transformation matrix"""
        translation = pose[0]
        quaternion = pose[1]

        # NOTE: PyBullet quaternions are x, y, z, w
        rotation = p.getMatrixFromQuaternion([quaternion[3], quaternion[0],
                                              quaternion[1], quaternion[2]])

        m = np.identity(4)
        m[0:3, 0:3] = rotation
        m.T[3, 0:3] = translation

        return np.ndarray(m)

    @staticmethod
    def matrix_to_pose(matrix):
        """Converts a transformation matrix to a pyBullet pose"""
        arr = np.array(matrix)
        translation = list(arr.T[3, 0:3])
        quaternion = p.getEulerFromQuaternion(arr[0:3, 0:3])

        # NOTE: PyBullet quaternions are x, y, z, w
        quaternion = [quaternion[1], quaternion[2],
                      quaternion[3], quaternion[0]]

        return translation, quaternion
