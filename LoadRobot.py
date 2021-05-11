import subprocess

robot = './robots/rov'
subprocess.call(['python', '../onshape_to_robot/onshape_to_robot.py', '../ROVSIM/robots/rov'])
