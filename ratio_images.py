#Filename: diference_images.py
#Description: This reads the header of a file to find what shutter was used. Then it moves the files to the correct folder based on the shutter.
#Date: August 3, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os


current = os.getcwd()
#gets current directory

base_path = '/Users/hkurtz/Desktop/post_flash/combine_im'
#shutter_b = os.path.join(base_path, 'isr_images')
os.chdir(base_path)


#gets the images to divide
#list_of_files= glob.glob('b*.fits')
filea='cal_image_A.fits'
fileb='cal_image_B.fits'
hdulist_1 = fits.open(filea)
data1=hdulist_1[1].data
data4=hdulist_1[4].data

#for i in range(len(list_of_files)):
#divides each individually from the first image read in.

hdulist = fits.open(fileb)
dataB4=hdulist[4].data
data=hdulist[1].data
 
all_dataA=np.concatenate((data4,data1),axis=0)
all_dataB=np.concatenate((dataB4,data),axis=0)
    #takes ratio and writes to file
datafinal=all_dataA/all_dataB
outfile= 'ratio_a_to_b_all.fits'
fits.writeto(outfile, datafinal,clobber=True)
    
    #datafinal4=data4/dataB4
    #outfile= 'ratio_' + list_of_files[0] + list_of_files[i]
    #isrfits.writeto(outfile, datafinal4,clobber=True)
    
#    if i!= 0:
        #plot the ratio
        
        #plot histogram of the ratio
 #       histogram=plt.hist(datafinal.ravel(), bins=1000,range=(0, 2))
  #      plt.show()
        
        #plot the ratio

        #plot histogram of the ratio
   #     histogram=plt.hist(datafinal4.ravel(), bins=1000,range=(0, 2))
    #    plt.show()
    #else:
     #   continue
        

    #stats
    #print(list_of_files[i])
    #print('mean 2012',np.mean(data1))
    #print('mean 201?',np.mean(data))
    #print('median',np.median(datafinal))
    #print ('mean',np.mean(datafinal))
    #print('std',np.std(datafinal))
    #print('Min:', np.min(datafinal))
    #print('Max:', np.max(datafinal))
       
    #print('median',np.median(datafinal4))
    #print ('mean',np.mean(datafinal4))
    #print('std',np.std(datafinal4))
    #print('Min:', np.min(datafinal4))
    #print('Max:', np.max(datafinal4))    
    #if np.mean(datafinal) < 0.97:
     #   print(list_of_files[i])
      #  print ('mean',np.mean(datafinal))
       # print ('mean',np.mean(datafinal4))
    #else:
     #   continue

    
   
hdulist.close()

hdulist_1.close()

os.chdir(current)