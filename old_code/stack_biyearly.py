#Filename: stack_year.py
#Description: This creates a mask of all the DQ and applies it to the image 
#Date: August 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

y2020=[]
y2018_19=[]
y2016_17=[]
y2014_15=[]
y2012_13=[]

current = os.getcwd()
base_path = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis_2/shutter_a/med_100/'
os.chdir(base_path)

files= glob.glob('*flt.fits')
print(len(files))


for i , f in enumerate(files):
    obs_date=fits.getheader(f)['DATE-OBS']
    #print(obs_date)
    if '2016' in obs_date:
        y2016_17.append(f)
    elif '2015' in obs_date:
        y2014_15.append(f)
    elif '2014' in obs_date:
        y2014_15.append(f)
    elif '2013' in obs_date:
        y2012_13.append(f)
    elif '2012' in obs_date:
        y2012_13.append(f)
    elif '2017' in obs_date:
        y2016_17.append(f)
    elif '2018' in obs_date:
        y2018_19.append(f)
    elif '2019' in obs_date:
        y2018_19.append(f)
    elif '2020' in obs_date:
        y2020.append(f)
    else:
        print('file',f,'has no year')

#print('2020',len(y2020))
print('2018 19',len(y2018_19))
print('2016 17',len(y2016_17))
print('2014 15',len(y2014_15))
print('2012 13',len(y2012_13))#, y2012)



def stack(list_of_files,out_file,error_file):
    hdr = fits.getheader(list_of_files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf = len(list_of_files)
    set_data=fits.getdata(list_of_files[0], 1)
    data_array_1 = np.empty((nf, ny, nx), dtype=float)
    data_array_2 = np.empty((nf, ny, nx), dtype=float)
    rms_1 = np.zeros(len(list_of_files), dtype=float)
    rms_2 = np.zeros(len(list_of_files), dtype=float)
    error_array_1 = np.empty((nf, ny, nx), dtype=float)
    error_array_2 = np.empty((nf, ny, nx), dtype=float)
    total_error_1 = np.zeros_like(set_data, dtype=float)
    total_error_2 = np.zeros_like(set_data, dtype=float)
    for i , f in enumerate(list_of_files):
        data_1=fits.getdata(f, 1)
        data_2=fits.getdata(f, 4)
        error_1=fits.getdata(f, 2)
        error_2=fits.getdata(f, 5)
        DQ_1=fits.getdata(f, 3)
        DQ_2=fits.getdata(f, 6)
        #data_1=data1
        #data_2=data2
        mask_1=np.zeros_like(data_1, dtype=bool)
        mask_2=np.zeros_like(data_2, dtype=bool)
        mask_1[DQ_1>=2**13]= True
        mask_2[DQ_2>=2**13]= True
        error_1[(0<DQ_1 & (DQ_1<2**13))]=0.00001
        error_2[(0<DQ_2& (DQ_2<2**13))]=0.00001
        error_1_sq=error_1**2
        error_2_sq=error_2**2
        masked_data_1= ma.array(data=data_1, mask=mask_1)
        masked_data_2= ma.array(data=data_2, mask=mask_2)
        data_array_1[i, :, :] = masked_data_1
        rms_1[i] = masked_data_1.std()
        data_array_2[i, :, :] = masked_data_2
        rms_2[i] = masked_data_2.std()
        total_error_1=total_error_1+(error_1_sq)
        total_error_2=total_error_2+(error_2_sq)
    sr_total_error_1=np.sqrt(total_error_1)
    sr_total_error_2=np.sqrt(total_error_2)
    fin_error_1=(sr_total_error_1/(float(len(list_of_files))))
    fin_error_2=(sr_total_error_2/(float(len(list_of_files))))
    image_median_1 = np.median(data_array_1, axis=0)
    image_median_2 = np.median(data_array_2, axis=0)
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(image_median_1))
    new_hdul.append(fits.ImageHDU(image_median_2))
    new_hdul.writeto(out_file, overwrite=True)
    #error
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(fin_error_1))
    new_hdul.append(fits.ImageHDU(fin_error_2))
    new_hdul.writeto(error_file,overwrite=True)

def main():
    stack(y2012_13,'y2012_13_median.fits', 'y2012_13_error.fits')
    stack(y2014_15,'y2014_15_median.fits', 'y2014_15_error.fits')
    stack(y2016_17,'y2016_17_median.fits', 'y2016_17_error.fits')
    stack(y2018_19,'y2018_19_median.fits', 'y2018_19_error.fits')
    #stack(y2020,'y2020_median.fits', 'y2020_error.fits')

main()
