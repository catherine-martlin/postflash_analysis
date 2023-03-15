#! /usr/bin/env python

""" Generates WFC3 UVIS post-flash reference files.

This script serves as a pipeline to create WFC3 UVIS post-flash reference
files and is a wrapper around several modules that perform subtasks of
the postflash reference file creation algorithm.

Authors
-------
    Catherine A. Martlin, 2023
    Heather Kurtz, 2016

Use
---

    This script is intended to be executed via the command line as such:
    ::
        python cal_uvis_make_postflash.py
"""

import argparse
import csv
import glob
import logging
import os

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import pandas as pd

from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail


@log_fail
@log_info
def cal_uvis_make_postflash_main():
    """The wrapper of the postflash creation.

    Parameters
    ----------
    """

    logging.info('')
    logging.info('')
    logging.info('')
    #years = [2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022]
    cadence = 1
    years=[2022]

    for curr in ['low','med']:
        for shutter in ['A','B']:
            for y in years:
                #A_shutter_paths_year, A_shutter_outfile_year, A_shutter_error_outfile_year, A_shutter_fullframe_pf_year = create_reference_file(y, working_directory, today, cadence, postflash_data, shutter=shutter)
                #B_shutter_paths_year, B_shutter_outfile_year, B_shutter_error_outfile_year, B_shutter_fullframe_pf_year = create_reference_file(y, working_directory, today, cadence, postflash_data, shutter=shutter)
                #stack(A_shutter_paths_year, A_shutter_outfile_year, A_shutter_error_outfile_year)
                #stack(B_shutter_paths_year, B_shutter_outfile_year, B_shutter_error_outfile_year)
                #filename = working_directory + '2022_fullframe_B_flc_stack_2023-03-13_med.fits'
                #prior_pf_file = working_directory + '6c82014gi_fls.fits' # 2021 low B
                #error_file = working_directory + '2022_fullframe_B_flc_error_stack_2023-03-13_med.fits'
                print(curr)
                print(shutter)
                print(y)
                #final_cal(filename, prior_pf_file, error_file, shutter_flag=shutter, current=curr)
    
    change_permissions(working_directory)
        
        
def stack(list_of_files, outfile, error_file):
    """This function will stack a set of FITS images and create a masked
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
    
    """
    # Gets the file size
    hdr = fits.getheader(list_of_files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf = len(list_of_files)
    
    # Setting up the empty data, rms, and error arrays
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
        #Read in the data and the DQs from the FITS for both extensions
        data_1 = fits.getdata(f, 1)
        data_2 = fits.getdata(f, 4)
        error_1 = fits.getdata(f, 2)
        error_2 = fits.getdata(f, 5)
        DQ_1 = fits.getdata(f, 3)
        DQ_2 = fits.getdata(f, 6)
        
        #Set the mask to a boolean array of the same size and shape as the data array
        mask_1 = np.zeros_like(data_1, dtype=bool)
        mask_2 = np.zeros_like(data_2, dtype=bool)
        
        #DQ = true because in the masking step a 1 will be masked and 0 allowed,
        #so the logic in this step must be reversed.
        #Create a mask for the data by setting any pixel with the flag value to True.
        mask_1[DQ_1>=2**13] = True
        mask_2[DQ_2>=2**13] = True
        error_1[(0<DQ_1 & (DQ_1<2**13))] = 0.00001
        error_2[(0<DQ_2& (DQ_2<2**13))] = 0.00001
        error_1_sq = error_1**2
        error_2_sq = error_2**2
        
        #Mask the data
        masked_data_1 = ma.array(data=data_1, mask=mask_1)
        masked_data_2 = ma.array(data=data_2, mask=mask_2)
        
        #Resets the data to an array for stacking
        data_array_1[i, :, :] = masked_data_1
        rms_1[i] = masked_data_1.std()
        data_array_2[i, :, :] = masked_data_2
        rms_2[i] = masked_data_2.std()
        total_error_1 = total_error_1+(error_1_sq)
        total_error_2 = total_error_2+(error_2_sq)
        
    sr_total_error_1 = np.sqrt(total_error_1)
    sr_total_error_2 = np.sqrt(total_error_2)
    fin_error_1 = (sr_total_error_1/(float(len(list_of_files))))
    fin_error_2 = (sr_total_error_2/(float(len(list_of_files))))
    
    #Create the median image
    image_median_1 = np.median(data_array_1, axis=0)
    image_median_2 = np.median(data_array_2, axis=0)
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(image_median_1))
    new_hdul.append(fits.ImageHDU(image_median_2))
    new_hdul.writeto(outfile, overwrite=True)
    
    #Error
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(fin_error_1))
    new_hdul.append(fits.ImageHDU(fin_error_2))
    new_hdul.writeto(error_file,overwrite=True)

def create_reference_file(year, working_directory, today, cadence, postflash_data, shutter):
    """This function will use input information to run the stack()
    function and create specific file names for the outputs.
    
    Parameters
    ----------
    year: int
        The year as an integer input.
    
    working_directory: str
        The path the files will be saved to. Needs trailing '\'
    
    today: str
        Format isn't important, needed to add date to filenames.
    
    cadence: int
        Number of years you want stacked. 1 = 1 year, 2 = biyearly, etc.
    
    postflash_data: pandas dataframe
        Dataframe to define the data you are using.
    
    Returns
    -------
    paths_year: str
        List of all FITS file needed to stack.
    
    outfile_year: str
        Filename for the outfile specified by year.
    
    error_outfile_year: str
        Filename for the error outfile specified by year.
    
    """
    fullframe_pf = postflash_data.loc[(postflash_data['subarray'] == False) & (postflash_data['shutter'] == '{}'.format(shutter)) & (postflash_data['flash_cur'] == 'MED') & (postflash_data['flash_dur'] == 100.0)]
    if cadence == 1:
        # Allow 2012 reference file to be created with some 2013 data
        if year == 2012:
            fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-11-14 00:00:00'.format(str(year+1)))]
        else:
        # All single years besides 2012
            fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-01-01 00:00:00'.format(str(year+1)))]
            paths_year = fullframe_pf_year.path.tolist()
            print(len(paths_year))
            outfile_year = '{}{}_fullframe_{}_flc_stack_{}.fits'.format(working_directory,str(year), shutter, today)
            error_outfile_year = '{}{}_fullframe_{}_flc_error_stack_{}.fits'.format(working_directory,str(year), shutter, today)
            print(outfile_year)
    else:
        # Create reference files of different year cadences
        fullframe_pf = postflash_data.loc[(postflash_data['subarray'] == False) & (postflash_data['shutter'] == '{}'.format(shutter)) & (postflash_data['flash_cur'] == 'MED') & (postflash_data['flash_dur'] == 100.0)]
        fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-01-01 00:00:00'.format(str(year+cadence)))]
        paths_year = fullframe_pf_year.path.tolist()
        print(len(paths_year))
        outfile_year = '{}{}_cadence{}_fullframe_{}_flc_stack_{}.fits'.format(working_directory,str(year), str(cadence), shutter, today)
        error_outfile_year = '{}{}_cadence{}_fullframe_{}_flc_error_stack_{}.fits'.format(working_directory,str(year), str(cadence), shutter, today)
        print(outfile_year)

    return paths_year, outfile_year, error_outfile_year, fullframe_pf_year

def final_cal(filename, prior_pf_file, error_file, shutter_flag, current):
    """This function will take the individual data and error fits files and
    create a single output FITS file. The header of that file will need to be
    updated and then the file can be delivered to CRDS.
    
    Parameters
    ----------
    filename: str
        
    
    prior_pf_file: str
        
    
    error_file: str
        
    
    shutter_flag: str
        A or B
    
    current: str
        'low' or 'med'
        
    """
    # Read in the data and error files
    pf_hdulist = fits.open(filename)
    data1 = pf_hdulist[0].data
    data4 = pf_hdulist[1].data
    pf_hdulist.close
    
    pf_hdulist_error = fits.open(error_file)
    error1 = pf_hdulist_error[0].data
    error4 = pf_hdulist_error[1].data
    pf_hdulist_error.close
    
    # Read in the header of the current reference file
    prior_pf_file = fits.open(prior_pf_file)
    pf_hdr = prior_pf_file[0].header
    pf_hdr1 = prior_pf_file[1].header
    pf_hdr2 = prior_pf_file[2].header
    pf_hdr3 = prior_pf_file[3].header
    pf_hdr4 = prior_pf_file[4].header
    pf_hdr5 = prior_pf_file[5].header
    pf_hdr6 = prior_pf_file[6].header
    
    # Divide by the exposure time, multiply by the scale factor and
    # multiply by the gain:
    if shutter_flag == 'B':
        if current == 'low':
            data1_scaled=((data1/100.00)*0.03623)*1.56
            data4_scaled=((data4/100.00)*0.03623)*1.56

            error1_cor = (((error1*0.03623))/100.0)
            error4_cor = (((error4*0.03623))/100.0)
        
        elif current == 'med':
            data1_scaled = ((data1/100.00)*0.03623*28.96)*1.56
            data4_scaled = ((data4/100.00)*0.03623*28.96)*1.56

            error1_cor = (((error1*0.03623)*28.96)/100.0)
            error4_cor = (((error4*0.03623)*28.96)/100.0)

    elif shutter_flag == 'A':
        if current == 'low':
            data1_scaled = ((data1/100.00)*0.03639)*1.56
            data4_scaled = ((data4/100.00)*0.03639)*1.56

            error1_cor = ((error1*0.03639)/100.0)
            error4_cor = ((error4*0.03639)/100.0)
            
        elif current == 'med':
            data1_scaled = ((data1/100.00)*0.03639*28.96)*1.56
            data4_scaled = ((data4/100.00)*0.03639*28.96)*1.56

            error1_cor = (((error1*0.03639)*28.96)/100.0)
            error4_cor = (((error4*0.03639)*28.96)/100.0)
    
    # Creates other extentions needed by the reference file
    dq1 = np.zeros((2070,4206))
    dq2 = np.zeros((2070,4206))
    
    # Add in rows of zeros for proper padding
    c = np.zeros((19, 4096))
    data1_19 = np.concatenate((data1_scaled,c), axis=0)
    data4_19 = np.concatenate((c,data4_scaled), axis=0)
    
    error1_19 = np.concatenate((error1_cor,c), axis=0)
    error4_19 = np.concatenate((c,error4_cor), axis=0)
    
    # Add 25 rows of zero to the edges
    d=np.zeros((2070,25 ))
    
    data4_f25 = np.concatenate((data4_19,d), axis=1)
    data1_f25 = np.concatenate((data1_19,d), axis=1)
    
    data4_e25 = np.concatenate((d,data4_f25), axis=1)
    data1_e25 = np.concatenate((d,data1_f25), axis=1)
    
    error4_f25 = np.concatenate((error4_19,d), axis=1)
    error1_f25 = np.concatenate((error1_19,d), axis=1)
    
    error4_e25 = np.concatenate((d,error4_f25), axis=1)
    error1_e25 = np.concatenate((d,error1_f25), axis=1)
    
    # Add in 60 columns of zeros in the middle
    idx=[]
    for i in range(60):
        idx.append(2073)

    data4_all = np.insert(data4_e25, idx, 0, axis=1)
    data1_all = np.insert(data1_e25, idx, 0, axis=1)
    
    error4_all = np.insert(error4_e25, idx, 0, axis=1)
    error1_all = np.insert(error1_e25, idx, 0, axis=1)
      
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

if __name__ == '__main__':

    module = os.path.basename(__file__).replace('.py', '')
    configure_logging(module)
    logging.getLogger(module)

    cal_uvis_make_postflash_main()
