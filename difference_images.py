#Filename: diference_images.py
#Description: This makes a difference image by subtracting two of the images from one another.
#Date: July 28, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

current = os.getcwd()
#gets current directory

base_path = '/Users/hkurtz/Desktop/post_flash/combine_im/med_test'
ori = os.path.join(base_path, 'ori_flt')
os.chdir(base_path)

list_of_files= glob.glob('*flt.fits')
#print(list_of_files)
#list_of_files=['cal_image_A.fits','cal_image_B.fits']
#cop_list=['wc52031oi_fls.fits','wc52031pi_fls.fits']
#list_of_files=['medium_cal_image_A.fits','medium_cal_image_B.fits']
#cop_list=['x4f1554li_fls.fits','x4f1554mi_fls.fits']

#hdulist_1 = fits.open(list_of_files[0])
#data1=hdulist_1[1].data
#data4=hdulist_1[4].data

for i in range(len(list_of_files)):
    #print(i)
#this gets the images to subtract and subtracts each individually from the first image read in.
    os.chdir(base_path)
    hdulist = fits.open(list_of_files[i])
    dataB4=hdulist[4].data
    data=hdulist[1].data
    all_data=np.concatenate((data,dataB4),axis=0)

    os.chdir(ori)
    hdulist_1 = fits.open(list_of_files[i])
    data1=hdulist_1[1].data
    data4=hdulist_1[4].data
    all_data1=np.concatenate((data1,data4),axis=0)
    
    print(np.shape(all_data1))

    #this subtracts the fits and writes a file of the difference
    datafinal=data1-data
    outfile= 'med_first_diff-' + list_of_files[i]
    fits.writeto(outfile, datafinal,clobber=True)
    
    datafinal4=data4-dataB4
    outfile= 'med_second_diff-' + list_of_files[i]
    fits.writeto(outfile, datafinal4,clobber=True)
    
    #this displays the difference file
    #if i!= 0:
        #this plots a histogram of the difference
    n, b, histogram=plt.hist(all_data1.ravel(), bins=1000,range=(-20, 20))
    plt.show()
        
        #this plots a histogram of the difference
    n, b, histogram=plt.hist(all_data.ravel(), bins=1000,range=(-20, 20))
    plt.show()
    #else:
     #   continue
    

    #print the stats of the image
    print('Min:', np.min(datafinal))
    print('Max:', np.max(datafinal))
    print('Mean', np.mean(datafinal))
    print('Stdev:', np.std(datafinal))    
    print('Min:', np.min(datafinal4))
    print('Max:', np.max(datafinal4))
    print('Mean', np.mean(datafinal4))
    print('Stdev:', np.std(datafinal4)) 
    
    
    

    
    
    hdulist.close()
    hdulist_1.close()