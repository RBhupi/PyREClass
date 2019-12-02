#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:40:07 2019
@author: Bhupendra Raut
"""

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plot
from wtclass import getWTClass

import matplotlib.cm as cm



fname = "/home/bhupendra/projects/screim/data/testdata/CPOL_20161210_0730_GRIDS_2500m.nc"

nc = Dataset(fname, "r")
data2d = np.squeeze(nc.variables["corrected_reflectivity"][0, 6, :, :]) #data*

plot.matshow(data2d)


data2d=data2d.filled(0)
data2d[data2d<0]=0

wt_class = getWTClass(data2d, res_km=2.5, conv_scale_km=20)

plot.matshow(wt_class)
plot.colorbar()
