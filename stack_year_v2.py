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
y2019=[]
y2018=[]
y2017=[]
y2016=[]
y2015=[]
y2014=[]
y2013=[]
y2012=[]

current = os.getcwd()
base_path = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis_2/shutter_b/med_100/'
os.chdir(base_path)

files= glob.glob('*flt.fits')
print(len(files))


for i , f in enumerate(files):
    obs_date=fits.getheader(f)['DATE-OBS']
    #print(obs_date)
    if '2016' in obs_date:
        y2016.append(f)
    elif '2015' in obs_date:
        y2015.append(f)
    elif '2014' in obs_date:
        y2014.append(f)
    elif '2013' in obs_date:
        y2013.append(f)
    elif '2012' in obs_date:
        y2012.append(f)
    elif '2017' in obs_date:
        y2017.append(f)
    elif '2018' in obs_date:
        y2018.append(f)
    elif '2019' in obs_date:
        y2019.append(f)
    elif '2020' in obs_date:
        y2020.append(f)
    else:
        print('file',f,'has no year')

print('2020',len(y2020))
print('2019',len(y2019))
print('2018',len(y2018))
print('2017',len(y2017))
print('2016',len(y2016))
print('2015',len(y2015))
print('2014',len(y2014))
print('2013',len(y2013))
print('2012',len(y2012))#, y2012)



def stack(list_of_files,out_file,error_file):
    hdr = fits.getheader(list_of_files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf = len(list_of_files)
    data_array_1 = np.empty((nf, ny, nx), dtype=float)
    data_array_2 = np.empty((nf, ny, nx), dtype=float)
    set_data=fits.getdata(list_of_files[0], 1)
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
    stack(y2012,'y2012_median_test.fits', 'y2012_error_test.fits')
    stack(y2013,'y2013_median_test.fits', 'y2013_error_test.fits')
    stack(y2014,'y2014_median_test.fits', 'y2014_error_test.fits')
    stack(y2015,'y2015_median_test.fits', 'y2015_error_test.fits')
    stack(y2016,'y2016_median_test.fits', 'y2016_error_test.fits')
    stack(y2017,'y2017_median_test.fits', 'y2017_error_test.fits')
    stack(y2018,'y2018_median_test.fits', 'y2018_error_test.fits')
    stack(y2019,'y2019_median_test.fits', 'y2019_error_test.fits')
    stack(y2020,'y2020_median_test.fits', 'y2020_error_test.fits')

main()
