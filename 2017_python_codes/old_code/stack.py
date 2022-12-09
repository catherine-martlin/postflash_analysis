#Filename: mask.py
#Description: This creates a mask of all the cosmic rays in the DQ and applies it to the image 
#Date: August 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

#gets current directory
current = os.getcwd()
base_path = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis/shutter_a/med_100/'
#shange directory to the base path
os.chdir(base_path)


#list of files
list_of_files= glob.glob('*flt.fits')

#gets the file size
hdr = fits.getheader(list_of_files[0], 1)
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']
nf = len(list_of_files)
set_data=fits.getdata(list_of_files[0], 1)
#makes empty array to be filled with data
data_array_1 = np.empty((nf, ny, nx), dtype=float)
data_array_2 = np.empty((nf, ny, nx), dtype=float)
error_array_1 = np.empty((nf, ny, nx), dtype=float)
error_array_2 = np.empty((nf, ny, nx), dtype=float)

#makes array for rms
rms_1 = np.zeros(len(list_of_files), dtype=float)
rms_2 = np.zeros(len(list_of_files), dtype=float)
total_error_1 = np.zeros_like(set_data, dtype=float)
total_error_2 = np.zeros_like(set_data, dtype=float)

#read in the data and the DQs from the .fits for both extensions
for i , f in enumerate(list_of_files):

    
    data_1=fits.getdata(f, 1)
    data_2=fits.getdata(f, 4)
    error_1=fits.getdata(f, 2)
    error_2=fits.getdata(f, 5)
    DQ_1=fits.getdata(f, 3)
    DQ_2=fits.getdata(f, 6)
   # print(error_1)

#set the mask to a boolean array of the same size and shape as the data array    
    mask_1=np.zeros_like(data_1, dtype=bool)
    mask_2=np.zeros_like(data_2, dtype=bool)
    #I set the DQ to true because in the masking step it is understood that 1 will be masked and zeros are fine so the logic in this step must be reversed.
#create the mask for the data by setting any pixel with the flag value to true.
    mask_1[(0<DQ_1)]= True
    mask_2[(0<DQ_2)]= True
    error_1[(0<DQ_1 & (DQ_1<2**13))]=0.00001
    error_2[(0<DQ_2& (DQ_2<2**13))]=0.00001
    error_1_sq=error_1**2
    error_2_sq=error_2**2
#mask the data
    masked_data_1= ma.array(data=data_1, mask=mask_1)
    masked_data_2= ma.array(data=data_2, mask=mask_2)
    
    
#resets the data to array for stacking
    data_array_1[i, :, :] = masked_data_1
    rms_1[i] = masked_data_1.std()
    data_array_2[i, :, :] = masked_data_2
    rms_2[i] = masked_data_2.std()


    total_error_1=total_error_1+(error_1_sq)
    total_error_2=total_error_2+(error_2_sq)


print(total_error_1)
sr_total_error_1=np.sqrt(total_error_1)
sr_total_error_2=np.sqrt(total_error_2)
fin_error_1=(sr_total_error_1/(float(len(list_of_files))))
fin_error_2=(sr_total_error_2/(float(len(list_of_files))))

#create the median image
image_median_1 = np.median(data_array_1, axis=0)
image_median_2 = np.median(data_array_2, axis=0)


#create mean image
image_mean_1= np.mean(data_array_1, axis=0)
image_mean_2= np.mean(data_array_2, axis=0)


#create wighted mean image
#get wights
weights_1 = 1./rms_1**2
weights_1 /= weights_1.sum()

weights_2 = 1./rms_2**2
weights_2 /= weights_2.sum()


weights_expand_1 = np.tile(weights_1[..., np.newaxis, np.newaxis], (1, ny, nx))
weighted_mean_1 = np.sum(weights_expand_1 * data_array_1, axis=0)

weights_expand_2 = np.tile(weights_2[..., np.newaxis, np.newaxis], (1, ny, nx))
weighted_mean_2 = np.sum(weights_expand_2 * data_array_2, axis=0)



#prints stats
print('image one', image_median_1.std(), image_mean_1.std(), weighted_mean_1.std())

print('image two', image_median_2.std(), image_mean_2.std(), weighted_mean_2.std())




#writes to the new file   

#median
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU(image_median_1))
new_hdul.append(fits.ImageHDU(image_median_2))
new_hdul.writeto('A_median_pf.fits',clobber=True)

#error
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU(fin_error_1))
new_hdul.append(fits.ImageHDU(fin_error_2))
new_hdul.writeto('A_error_pf.fits',clobber=True)
    
#mean
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU(image_mean_1))
new_hdul.append(fits.ImageHDU(image_mean_2))
new_hdul.writeto('A_mean_pf.fits',clobber=True)


#weighted mean
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU(weighted_mean_1))
new_hdul.append(fits.ImageHDU(weighted_mean_2))
new_hdul.writeto('A_weighted_mean_pf.fits',clobber=True)



#changes bac to current directory
os.chdir(current)
