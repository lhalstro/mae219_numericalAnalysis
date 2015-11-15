
import numpy as np
import matplotlib.pyplot as plt
import os

n = 257
Re = 400
d = 1
U = 1.

nu = d*U/Re
print('nu', nu)
dx = d/n
print('dx', dx)
dt = 1*dx/U
print('dt', dt)