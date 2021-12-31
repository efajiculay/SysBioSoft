"""

This module is for initial test

"""


from numpy import linspace, sin, cos, pi, vstack, random
label = "Thank you for using BioSANS".split()
theta = linspace(0, 2 * pi, 100)
x = 16 * (sin(theta) ** 3)
tup = (x,)
for i in range(5):
    y = 13 * cos(theta)-5*cos(2*theta)-2*cos(3*theta)-cos(4*theta) \
        + 5*random.rand()
    tup = tup + (y,)
data = [vstack(tup).T]
