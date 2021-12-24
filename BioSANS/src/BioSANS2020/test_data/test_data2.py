"""

This module is for initial test

"""


from numpy import linspace, sin, cos, pi, vstack
label = ["Thank you for using BioSANS"]
theta = linspace(0, 2 * pi, 100)
x = 16 * (sin(theta) ** 3)
y = 13 * cos(theta)-5*cos(2*theta)-2*cos(3*theta)-cos(4*theta)
data = [vstack((x, y)).T]
