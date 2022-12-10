from astropy.io import fits
import numpy as np
import numpy.ma as ma

def stack(list_of_files,outfile,error_file):
    '''This function will stack a set of FITS images and create a masked
    median file along with an error file.
    
    Parameters
    ----------
    list_of_files: str
        This is the list you create from the path column of the
        subset of the pandas database.

    outfile: str
        This is the path and filename you want for your outfile.

    error_file: str
    This is the path and filename you want for your
        outfile of the calcuated error.
    
    '''
    #gets the file size
    hdr = fits.getheader(list_of_files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf = len(list_of_files)
    # Setting up the empty data, rms, and error arrays and getting the data
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
        #read in the data and the DQs from the .fits for both extensions
        data_1=fits.getdata(f, 1)
        data_2=fits.getdata(f, 4)
        error_1=fits.getdata(f, 2)
        error_2=fits.getdata(f, 5)
        DQ_1=fits.getdata(f, 3)
        DQ_2=fits.getdata(f, 6)
        #set the mask to a boolean array of the same size and shape as the data array
        mask_1=np.zeros_like(data_1, dtype=bool)
        mask_2=np.zeros_like(data_2, dtype=bool)
        #I set the DQ to true because in the masking step it is understood that 1 will
        #be masked and zeros are fine so the logic in this step must be reversed.
        #create the mask for the data by setting any pixel with the flag value to true.
        mask_1[DQ_1>=2**13]= True
        mask_2[DQ_2>=2**13]= True
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
    sr_total_error_1=np.sqrt(total_error_1)
    sr_total_error_2=np.sqrt(total_error_2)
    fin_error_1=(sr_total_error_1/(float(len(list_of_files))))
    fin_error_2=(sr_total_error_2/(float(len(list_of_files))))
    #create the median image
    image_median_1 = np.median(data_array_1, axis=0)
    image_median_2 = np.median(data_array_2, axis=0)
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(image_median_1))
    new_hdul.append(fits.ImageHDU(image_median_2))
    new_hdul.writeto(outfile, overwrite=True)
    #error
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(fin_error_1))
    new_hdul.append(fits.ImageHDU(fin_error_2))
    new_hdul.writeto(error_file,overwrite=True)
