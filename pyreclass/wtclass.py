#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 23:12:19 2019
@author: Bhupendra Raut
"""

import numpy as np
from math import log
from atwt import atwt2d



def getWTClass(dbz_data, res_km, conv_scale_km=20):
    """ Compute scan-by-scan ATWT of radar volume data and classify echoes.
    
     Converts dBZ to rain rates using standard Z-R relationship.
     This is to transform the normally distributed dBZ to gamma-like distribution.
     @param \code{vol_data} 3D array containing radar data. Last dimension should be levels.
     @param \code{res_km} resolution of the radar data in km.
     @param \code{conv_scale} approximate scale break (in km) between convective and stratiform scales.
     @return Sum of wavelets upto \code{conv_scale} for each scan.
     @export
     @seealso \code{\link{getWTSum}}"""
    
    #save the mask for missing data.
    #data_mask <- ifelse(test = is.na(dbz_data), yes = NA, no = 0)
    dbz_data[dbz_data==np.nan]=0
    scale_break = getScaleBreak(res_km, conv_scale_km)
    dbz_data_t = dbz2rr(dbz_data) #transform the dbz data
    wt_sum = getWTSum(dbz_data_t, scale_break)
    wt_class = labelClasses(wt_sum, dbz_data)
    return(wt_class)




def labelClasses(wt_sum, vol_data) :
    """ Lables 1. stratiform, 2. intense/heavy convective  and 
    3. moderate+transitional convective regions using given thersholds."""
     
    conv_wt_threshold = 5 #WT value more than this is strong convection
    tran_wt_threshold = 2 #WT value for moderate convection
    min_dbz_threshold = 10 #pixels below this value are not classified
    conv_dbz_threshold = 30 # pixel below this value are not convective

    #I first used negative numbers to annotate the categories. Then multiply it by -1.
    wt_class = np.where((wt_sum>=conv_wt_threshold) & (vol_data>=conv_dbz_threshold), -2, 0)
    wt_class = np.where((wt_sum<conv_wt_threshold) & (wt_sum>=tran_wt_threshold) &\
                        (vol_data>=conv_dbz_threshold), -3, wt_class)
    wt_class = np.where((wt_class==0) & (vol_data>=min_dbz_threshold), -1, wt_class)
    
    wt_class = wt_class*-1
    wt_class = np.where((wt_class==0), np.nan, wt_class)

    return(wt_class)




def dbz2rr(dbz, ZRA=200, ZRB=1.6):
    """Uses standard values for ZRA=200 and ZRB=1.6.
    @param dbz array, vector or matrix of reflectivity in dBZ
    @return rr rain rate in \code{mm/hr}"""
    rr = ((10.0**(dbz/10.0))/ZRA)**(1.0/ZRB)
    return(rr)

def getScaleBreak(res_km, conv_scale_km):
    """#' compute scale break for convection and stratiform regions.
    WT will be computed upto this scale and features will be designated as convection.
    @param res_km resolution of the image.
    @param conv_scale_km expected size of spatial variations due to convection.
    @return dyadic integer scale break in pixels. """
    scale_break = log((conv_scale_km/res_km))/log(2)+1
    return(round(scale_break))





def getWTSum(vol_data, conv_scale):
    """returns sum of WT upto given scale.
    
    Works with both 2d scans and Volume data."""
    dims = vol_data.shape

    #if data is 2d
    if(len(dims)==2):
        wt = atwt2d(vol_data, max_scale = conv_scale)
        wt_sum = np.sum(wt, axis=(0))
    else:			#else for volume data
    	num_levels = dims[3]
    	wt_sum = np.zeros(dims)

    	for lev in range(num_levels):
            if vol_data[:, :, lev].max < 1:
                next() #this needs reviewing
            wt=atwt2d(vol_data[lev, :, :], max_scale = conv_scale)

            #sum all the WT scales.
            wt_sum[lev, :, :] = np.sum(wt, axis=(0))


    #negative WT do not corresponds to convection in radar
    #wt_sum[wt_sum<0] = 0 #check calling function
    return(wt_sum)

