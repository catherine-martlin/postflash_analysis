#file: isr_image
#description: makes images for isr
#Date: January 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
import os
import numpy as np
import glob

current = os.getcwd()
base_path = '/Users/hkurtz/Desktop/post_flash/combine_im'
os.chdir(base_path)

for file in glob.glob('A_median_pf_100.fits'):

	hdulist_a1 = fits.open(file)
	dataA1=hdulist_a1[0].data
	dataA4=hdulist_a1[1].data

	all_data=np.concatenate((dataA1,dataA4), axis=0)
	all_data_norm=all_data/np.mean(all_data)
	out_file='full'+file
	fits.writeto(out_file, all_data_norm,clobber=True)
	
	print(file)
	print((np.min(all_data))*1.5)
	print((np.max(all_data))*1.5)
	print((np.mean(all_data))*1.5)
	print((np.std(all_data))*1.5)

os.chdir(current)