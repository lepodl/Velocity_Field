# -*- coding: utf-8 -*- 
# @Time : 2021/10/18 11:44 
# @Author : lepold
# @File : 3D_brain.py

import h5py
import numpy as np
from mayavi import mlab

# file = h5py.File(r"E:\PycharmProjects\brain_information_flow\data\raw_data\100206.mat", "r")
# # file = loadmat(r"E:\PycharmProjects\brain_information_flow\data\raw_data\100206.mat")
# xyz = file["BrainImg"][:]
# xyz = xyz.T
# print("xyz shape", xyz.shape)
# x, y, z = np.nonzero(xyz[:,:,:,20])
# mlab.points3d(x, y, z, )
# mlab.show()

file = h5py.File("../../Data/raw_data/DTI_voxel_network.mat", "r")
xyz = file["dti_xyz2"][:]
print("xyz shape", xyz.shape)
x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
mlab.points3d(x, y, z, )
mlab.show()