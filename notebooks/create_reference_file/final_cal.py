"""Create the final post-flash reference file from the stacked fits files. 


Authors
-------
    Catherine Martlin, 2023
    Heather Kurtz, 2016
Use
---
    This module is intended to be called by ``.py``
    as part of the UVIS Post-flash reference file creation.
"""

import csv
import glob
import os

from astropy.io import fits
import numpy as np

def final_cal(filename, prior_pf_file, error_file, shutter_flag, current):
    
    # Read in the two fits files
    pf_hdulist= fits.open(filename)
    data1=pf_hdulist[0].data
    data4=pf_hdulist[1].data
    
    pf_hdulist.close
    
    pf_hdulist_error= fits.open(error_file)
    error1=pf_hdulist_error[0].data
    error4=pf_hdulist_error[1].data
    
    pf_hdulist_error.close
    
    # read in the header of the current reference file
    
    prior_pf_file = fits.open(prior_pf_file)
    pf_hdr = prior_pf_file[0].header
    pf_hdr1 = prior_pf_file[1].header
    pf_hdr2 = prior_pf_file[2].header
    pf_hdr3 = prior_pf_file[3].header
    pf_hdr4 = prior_pf_file[4].header
    pf_hdr5 = prior_pf_file[5].header
    pf_hdr6 = prior_pf_file[6].header
    
    # Divide by the exposure time, multiply by the scale factor and by the gain:
    if shutter_flag == 'B':
        if current == 'low':
            data1_scaled=((data1/100.00)*0.03623)*1.56
            data4_scaled=((data4/100.00)*0.03623)*1.56

            error1_cor=(((error1*0.03623))/100.0)
            error4_cor=(((error4*0.03623))/100.0)
        
        elif current == 'med':
            data1_scaled=((data1/100.00)*0.03623*28.96)*1.56
            data4_scaled=((data4/100.00)*0.03623*28.96)*1.56

            error1_cor=(((error1*0.03623)*28.96)/100.0)
            error4_cor=(((error4*0.03623)*28.96)/100.0) 

    elif shutter_flag == 'A':
        if current == 'low':
            data1_scaled=((data1/100.00)*0.03639)*1.56
            data4_scaled=((data4/100.00)*0.03639)*1.56

            error1_cor=((error1*0.03639)/100.0)
            error4_cor=((error4*0.03639)/100.0)
            
        elif current == 'med':
            data1_scaled=((data1/100.00)*0.03639*28.96)*1.56
            data4_scaled=((data4/100.00)*0.03639*28.96)*1.56

            error1_cor=(((error1*0.03639)*28.96)/100.0)
            error4_cor=(((error4*0.03639)*28.96)/100.0)
    
    # Creates other extentions needed by the reference file
    dq1=np.zeros((2070,4206))
    dq2=np.zeros((2070,4206))
    
    # Add in rows of zeros for proper padding
    c=np.zeros((19, 4096))
    data1_19=np.concatenate((data1_scaled,c), axis=0)
    
    data4_19=np.concatenate((c,data4_scaled), axis=0)
    
    error1_19=np.concatenate((error1_cor,c), axis=0)
    
    error4_19=np.concatenate((c,error4_cor), axis=0)
    
    # Add 25 rows of zero to the edges
    d=np.zeros((2070,25 ))
    
    data4_f25=np.concatenate((data4_19,d), axis=1)
    data1_f25=np.concatenate((data1_19,d), axis=1)
    
    data4_e25=np.concatenate((d,data4_f25), axis=1)
    data1_e25=np.concatenate((d,data1_f25), axis=1)
    
    error4_f25=np.concatenate((error4_19,d), axis=1)
    error1_f25=np.concatenate((error1_19,d), axis=1)
    
    error4_e25=np.concatenate((d,error4_f25), axis=1)
    error1_e25=np.concatenate((d,error1_f25), axis=1)
    
    # Add in 60 columns of zeros in the middle
    idx=[]
    for i in range(60):
        idx.append(2073)

    data4_all=np.insert(data4_e25, idx, 0, axis=1)
    data1_all=np.insert(data1_e25, idx, 0, axis=1)
    
    error4_all=np.insert(error4_e25, idx, 0, axis=1)
    error1_all=np.insert(error1_e25, idx, 0, axis=1)
      
    print('data_4')
    print('Min:', np.min(data4_all))
    print('Max:', np.max(data4_all))
    print('Mean', np.mean(data4_all))
    print('Stdev:', np.std(data4_all))      
        
    print('data_1')
    print('Min:', np.min(data1_all))
    print('Max:', np.max(data1_all))
    print('Mean', np.mean(data1_all))
    print('Stdev:', np.std(data1_all))    
    
    # All data together with the 7 extensions and the header
    cal_out_file = filename[:-5] + '_final_cal.fits'
    
    prihdu = fits.PrimaryHDU(header=pf_hdr)
    
    single_extension1 = fits.ImageHDU(data = data1_all.astype(np.float32), header = pf_hdr1)
    single_extension2 = fits.ImageHDU(data = error1_all.astype(np.float32), header = pf_hdr2)
    single_extension3 = fits.ImageHDU(data = dq1.astype(np.float32), header = pf_hdr3)
    
    single_extension4 = fits.ImageHDU(data = data4_all.astype(np.float32), header = pf_hdr4)
    single_extension5 = fits.ImageHDU(data = error4_all.astype(np.float32), header = pf_hdr5)
    single_extension6 = fits.ImageHDU(data = dq2.astype(np.float32), header = pf_hdr6)
    
    all_extensions = [prihdu, single_extension1, single_extension2, 
                      single_extension3, single_extension4, single_extension5, 
                      single_extension6]
    
    myhdulist = fits.HDUList(all_extensions)
    
    myhdulist[0]._bitpix = 16
    
    myhdulist.writeto(cal_out_file, overwrite=True)
    
    prior_pf_file.close