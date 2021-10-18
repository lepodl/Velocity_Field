# -*- coding: utf-8 -*- 
# @Time : 2021/10/14 17:19 
# @Author : lepold
# @File : process_lfp.py

import numpy as np
import os
import h5py
import sparse
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt
import time


# Process lfp pipeline.

def np_move_avg(a, n=10, mode="same"):
    tmp = []
    for i in range(a.shape[1]):
        tmp.append(np.convolve(a[:, i], np.ones((n,)) / n, mode=mode))
    return np.stack(tmp, axis=1)


def show(data):
    seed = np.random.randint(low=0, high=data.shape[1], size=(4,))
    fig, ax = plt.subplots(4, 1, figsize=(8, 5), dpi=300)
    ax = ax.flatten()
    for i in range(4):
        ax[i].plot(data[:, seed[i]])
    fig.savefig("./show_lfp.png")
    return


# Step 1
start = time.time()
project_path = "../"
data_path = os.path.join(project_path, "Data/raw_data")
invert_index = np.load("../Data/raw_data/invert_index.npy")
lfp = np.load("../Data/raw_data/lfp_200hz.npy")
lfp = lfp[:, invert_index] # vmean shape: (1600, 22703)
# lfp = np_move_avg(lfp, n=3, mode="valid")
# show(lfp)
# lfp = lfp[700:1200, :]
# lfp = lfp[::3, :]
t, num = lfp.shape
bad_channel = np.where(lfp<-65.)[1]
bad_channel = np.unique(bad_channel)

print(f"lfp min: {np.sort(lfp.flatten())[:5]}")
print(f"lfp max: {np.sort(lfp.flatten())[-5:]}")
lfp = lfp.astype(np.float32)
print("flp shape", lfp.shape)

coordinate_path = os.path.join(data_path, "DTI_voxel_network.mat")
valid_data = loadmat(os.path.join(data_path, "All_gmwmi_dat.mat"))
valid_label = valid_data["label"][:].squeeze()
valid_coord = valid_data["xyz"][:]
valid_coord = valid_coord.T
xyz = h5py.File(coordinate_path, "r")["dti_xyz2"][:]
label = h5py.File(coordinate_path, "r")["dti_label_num"][:].squeeze()
temp = np.where(~np.isin(valid_label, label))[0]

coordinate_min = np.empty(3, dtype=np.int64)
coordinate_max = np.empty(3, dtype=np.int64)
for i in range(xyz.shape[1]):
    print(f"{i}th coordinate: min {valid_coord[:, i].min()}, max {valid_coord[:, i].max()}")
    coordinate_min[i] = valid_coord[:, i].min()
    coordinate_max[i] = valid_coord[:, i].max()

xyz = xyz - coordinate_min  # shape(22703, 3)
valid_coord = valid_coord - coordinate_min


for i in range(xyz.shape[1]):
    print(f"{i}th coordinate: min {xyz[:, i].min()}, max {xyz[:, i].max()}")
    coordinate_min[i] = xyz[:, i].min()
    coordinate_max[i] = xyz[:, i].max()
if xyz.shape[1] == 3:
    xyz = xyz.T

linear_interpolate_index = []
size = coordinate_max + 1
bad_xyz = xyz[:, bad_channel]
for i in temp:
    x, y, z = valid_coord[i]
    linear_interpolate_index.append(x + y*size[0] + z*size[0]*size[1]+1)
for i in range(bad_xyz.shape[1]):
    x, y, z = bad_xyz[:, i]
    linear_interpolate_index.append(x + y*size[0] + z*size[0]*size[1]+1)
linear_interpolate_index = np.array(linear_interpolate_index, dtype=np.int32)
linear_interpolate_index = np.unique(linear_interpolate_index)

shape = [coordinate_max[0] + 1, coordinate_max[1] + 1, coordinate_max[2] + 1, t]
print(f"Sparse matrix shape {shape}")
coords = np.empty([4, lfp.shape[1] * lfp.shape[0]], dtype=np.int64)
coords[(0, 1, 2), :] = np.broadcast_to(xyz[:, :, None], (xyz.shape[0], xyz.shape[1], t)).reshape(
    (3, -1))
coords[3, :] = np.broadcast_to(np.arange(t, dtype=np.int64), (lfp.shape[1], t)).reshape(-1)
data = lfp.T.reshape((-1))
# coords = tuple(x for x in coords)
coo = sparse.COO(coords=coords, data=data, shape=shape)
data = coo.todense()
savemat(os.path.join(data_path, "lfp_200hz.mat"), mdict={"brain_image": data, "interp_idnex":linear_interpolate_index})
print(f"Done! cost time {time.time() - start:.2f}")

