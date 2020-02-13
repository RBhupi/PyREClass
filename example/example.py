"""
Created on Thu Oct 19 14:40:07 2019
@author: Bhupendra Raut
@modified: Valentin Louf
@modification: Thu Feb 13 11:23:00 2020
"""
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from pyreclass import getWTClass

fname = "/home/bhupendra/projects/screim/data/testdata/CPOL_20161210_0730_GRIDS_2500m.nc"

with netCDF4.Dataset(fname, "r") as ncid:
    data2d = np.squeeze(ncid["corrected_reflectivity"][0, 6, :, :]) #data*
    data2d = data2d.filled(0)
    data2d[data2d < 0] = 0

wt_class = getWTClass(data2d, res_km=2.5, conv_scale_km=20)

fig, ax = pl.subplots(1, 2, figsize=(12, 5))
ax = ax.ravel()

ax[0].pcolormesh(data2d)
ax[1].pcolormesh(wt_class)

pl.show()