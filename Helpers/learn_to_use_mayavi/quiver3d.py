# -*- coding: utf-8 -*- 
# @Time : 2021/10/11 16:11 
# @Author : lepold
# @File : quiver3d.py

import numpy as np
from mayavi import mlab
from mayavi.mlab import *

x, y, z = np.mgrid[-2:3, -2:3, -2:3]
r = np.sqrt(x ** 2 + y ** 2 + z ** 4)
u = y * np.sin(r) / (r + 0.001)
v = -x * np.sin(r) / (r + 0.001)
w = np.zeros_like(z)
obj = quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
mlab.show()