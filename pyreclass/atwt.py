#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:34:17 2019
@author: Bhupendra Raut
"""


import numpy as np
from math import log, floor


def atwt2d(data2d, max_scale=-1):
    """Computes A trous wavelet transform (ATWT)
    Computes ATWT of the 2d array up to \code{max_scale}.
    If \code{max_scale} is outside the boundaries, number of scales will be reduced.
    Data is mirrored at the boundaries.'Negative WT are removed. Not tested for non-square data.
    @param data2d 2d image as array or matrix.
    @param max_scale computes wavelets up to \code{max_scale}. Leave blank for maximum possible scales.
    @return array containing ATWT of input image with added 3rd dimention for scales.
    @note Need to break this into smaller functions.
    @author Bhupendra Raut and Dileep M. Puranik
    @seealso Press et al. (1992) Numerical Recipes in C."""

    dims = data2d.shape
    max_possible_scales = getMaxScales(dims)
    if(max_scale<0 or max_possible_scales < max_scale):
        max_scale = max_possible_scales

    ny = dims[0]
    nx = dims[1]


    wt = np.zeros((max_scale, ny, nx))
    temp1 = np.zeros(dims)
    temp2 = np.zeros(dims)

    sf=(0.0625, 0.25, 0.375) #scaling function

    #start wavelet loop
    for scale in range(1, max_scale+1):
        print(scale)
        x1 = 2**(scale-1)
        x2 = 2 * x1

        #Row-wise smoothing
        for i in range(0, nx):
            #find the indices for prev and next points on the line
            prev2 = abs(i-x2)
            prev1 = abs(i-x1)
            next1 = (i+x1)
            next2 = (i+x2)

            #If these indices are outside the image, "mirror" them
            #Sometime this causes issues at higher scales.
            if next1 > nx-1:
                next1 = 2*(nx-1) - next1

            if next2 > nx-1:
                next2 = 2*(nx-1) - next2

            if prev1<0 or prev2 <0 :
                prev1 = next1
                prev2 = next2

            for j in range(0, ny) :
                left2  =  data2d[j, prev2]
                left1  =  data2d[j, prev1]
                right1  =  data2d[j, next1]
                right2  =  data2d[j, next2]
                temp1[j, i]  =  sf[0] * (left2+right2) + sf[1] * (left1 + right1) + sf[2] * data2d[j, i]

        #Column-wise smoothing
        for i in range(0, ny):

            prev2 = abs(i-x2)
            prev1 = abs(i-x1)
            next1 = (i+x1)
            next2 = (i+x2)

            #If these indices are outside the image use next values
            if next1 > ny-1 :
                next1 = 2*(ny-1) - next1

            if next2 > ny-1 :
                next2 = 2*(ny-1) - next2

            if prev1<0 or prev2 <0 :
                prev1 = next1
                prev2 = next2
            #print ("scale "+ str(scale) + " i=" + str(i), " prev2=" +str(prev2)+ " prev1=" + str(prev1) + " next2=" + str(next2))
            for j in range(0, nx) :
                top2  =  temp1[prev2, j]
                top1  =  temp1[prev1, j]
                bottom1  =  temp1[next1, j]
                bottom2  =  temp1[next2, j]
                temp2[i, j]  =  sf[0] * (top2+bottom2) + sf[1] * (top1 + bottom1) + sf[2] * temp1[i, j]

        wt[scale-1, :, :] = data2d - temp2
        data2d[:] = temp2

    return(wt)




def getMaxScales(dims):
    """Calculate the mximum possible scale of ATWT for given dimensions.
    @param data_dim output of the \code{dim(data2d)} for given matrix or array.
    @return integer value of the maximum scale.
    """
    min_dims = min(dims)
    max_scale = log(min_dims)/log(2)
    return(floor(max_scale))
