import numpy as np
import numba as nb

SCALE = 0.0254

THRUST_MAX = 3.71
THRUST_MIN = -2.9

X_COMP = np.sin(7 * np.pi / 18)
Y_COMP = np.cos(7 * np.pi / 18)


@nb.njit
def min_scale(array):
    return min(abs(THRUST_MAX / max(array[array != 0])), abs(THRUST_MIN / min(array[array != 0])))


@nb.njit
def mag(vec):  # gets the magnitude of the vector really fast for small size vectors
    return np.dot(vec, vec) ** 0.5


class ThrustMapper:
    def __init__(self):
        self.com = np.array([0.0, 0.0, 1.4]) * SCALE  # This is the COM position
        self.thruster_locations = np.array([[+4.438,  5.679,  0.00],           # Thruster 1
                                            [-4.438,  5.679,  0.00],           # Thruster 2
                                            [-4.438, -5.679,  0.00],           # Thruster 3
                                            [+4.438, -5.679,  0.00],           # Thruster 4
                                            [+7.500,  7.313, -2.25],           # Thruster 5
                                            [-7.500,  7.313, -2.25],           # Thruster 6
                                            [-7.500, -7.313, -2.25],           # Thruster 7
                                            [+7.500, -7.313, -2.25]]) * SCALE  # Thruster 8

        self.thruster_directions = np.array([[0, 0, 1],  # Thruster 1
                                             [0, 0, 1],  # Thruster 2
                                             [0, 0, 1],  # Thruster 3
                                             [0, 0, 1],  # Thruster 4
                                             [+X_COMP, -Y_COMP, 0],   # Thruster 5
                                             [-X_COMP, -Y_COMP, 0],   # Thruster 6
                                             [-X_COMP,  Y_COMP, 0],   # Thruster 7
                                             [+X_COMP,  Y_COMP, 0]])  # Thruster 8

        self.origin = None
        self.torque_map = None
        self.thruster_output_map = None
        self.pseudo_inverse = None

        self.change_origin(np.array([0.0, 0.0, 1.1]))

        self.fine = True
        self.fine_scale = 1.5
        self.mix = 0.70

    def change_origin(self, offset):
        self.origin = self.thruster_locations - self.com + offset * SCALE

        self.torque_map = np.cross(self.origin, self.thruster_directions)
        self.thruster_output_map = np.concatenate((self.thruster_directions, self.torque_map), axis=1)
        self.pseudo_inverse = np.linalg.pinv(self.thruster_output_map)

    def thruster_output(self, commanded_force, commanded_torque):
        f_mag = mag(commanded_force)
        t_mag = mag(commanded_torque)

        if f_mag > 1.0:
            commanded_force = commanded_force / f_mag
        if t_mag > 1.0:
            commanded_torque = commanded_torque / t_mag

        force = np.matmul(commanded_force, self.pseudo_inverse)
        torque = np.matmul(commanded_torque, self.pseudo_inverse)

        if self.fine:
            force *= 3.7301 * self.fine_scale
            torque *= 1.031 * self.fine_scale  # Thrust envelop inscribed circle radius
        else:
            if force.any():
                force *= min_scale(force)
            if torque.any():
                torque *= min_scale(torque)

        result = force * f_mag * self.mix + torque * t_mag * (1 - self.mix)

        return result

    @staticmethod
    @nb.njit
    def thruster_output2(commanded_force, commanded_torque, pseudo_inverse, mix, fine_scale, fine=False):
        f_mag = mag(commanded_force)
        t_mag = mag(commanded_torque)

        if f_mag > 1.0:
            commanded_force = commanded_force / f_mag
        if t_mag > 1.0:
            commanded_torque = commanded_torque / t_mag

        force = np.dot(commanded_force, pseudo_inverse)
        torque = np.dot(commanded_torque, pseudo_inverse)

        if fine:
            force *= 3.7301 * fine_scale
            torque *= 1.031 * fine_scale  # Smallest sphere that fits in the thrust envelopes
        else:
            if force.any():
                force *= min_scale(force)
            if torque.any():
                torque *= min_scale(torque)

        result = force * f_mag * mix + torque * t_mag * (1 - mix)

        if not fine and result.any():
            result *= min_scale(result)

        return result

    @staticmethod
    @nb.njit
    def thrust_to_pwm(thrust_val):
        if thrust_val < -.04:
            pwm = 0.018 * (thrust_val ** 3) + 0.117 * (thrust_val ** 2) + 0.4981 * thrust_val - 0.09808
        elif thrust_val > .04:
            pwm = 0.0095 * (thrust_val ** 3) - 0.0783 * (thrust_val ** 2) + 0.4004 * thrust_val + 0.0986
        else:
            pwm = 0.0
        return pwm

    @staticmethod
    @nb.njit
    def pwm_to_thrust(pwm):
        if pwm < -.1:
            thrust_out = -.8944 * (pwm ** 3) - 2.971 * (pwm ** 2) + 0.9844 * pwm + .1005
        elif pwm > .1:
            thrust_out = -1.1095 * (pwm ** 3) + 3.9043 * (pwm ** 2) + 1.1101 * pwm - 0.113
        else:
            thrust_out = 0.0

        return thrust_out


if __name__ == '__main__':
    tm = ThrustMapper()

    force_input = np.array([1, 0, 0, 0, 0, 0], dtype=np.float)  # X Y Z
    torque_input = np.array([0, 0, 0, 0, 0, 1], dtype=np.float)  # Ro Pi Ya

    pwm_values = tm.thruster_output(force_input, torque_input)
    pwm_values2 = tm.thruster_output(force_input, torque_input)

    result1 = np.matmul(pwm_values, tm.thruster_output_map)
    result2 = np.matmul(pwm_values2, tm.thruster_output_map)

    print(list(np.around(np.array(pwm_values), 2)))
    print(list(np.around(np.array(pwm_values2), 2)))

    print(list(np.around(np.array(result1), 2)))
    print(list(np.around(np.array(result2), 2)))
