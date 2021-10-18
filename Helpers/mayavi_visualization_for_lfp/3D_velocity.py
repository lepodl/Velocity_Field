# -*- coding: utf-8 -*- 
# @Time : 2021/10/16 22:32 
# @Author : lepold
# @File : 3D_velocity.py

import numpy as np
from scipy.io import loadmat
from mayavi import mlab
from mayavi.mlab import *

velocity = loadmat("../../Data/process_data/velocity_of_lfp_200hz.mat")
vx = velocity["Ux"]
vy = velocity["Uy"]
vz = velocity["Uz"]
x, y, z, t = vx.shape

seed_time = 12
vx = vx[:, :, :, seed_time]
vy = vy[:, :, :, seed_time]
vz = vz[:, :, :, seed_time]
xx, yy, zz = np.mgrid[0:x, 0:y, 0:z]
obj = quiver3d(xx, yy, zz, vx, vy, vz, line_width=3, scale_factor=1)
mlab.show()