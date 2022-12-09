#Filename: sort_files.py
#Description: This reads the header of a file to find what shutter was used. Then it moves the files to the correct folder based on the shutter.
#Date: July 28, 2016
#Author: Heather Kurtz
from astropy.io import fits
import shutil
#to move fiels
import os
#to move files
import glob
#gets lots of files and gives their files names
import pdb
#this is the python debugging tool it is really cool and very helpful


current = os.getcwd()
#gets current directory

base_path = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis_2'
shutter_a = os.path.join(base_path, 'shutter_a')
shutter_b = os.path.join(base_path, 'shutter_b')

os.chdir(base_path)
#takes you to the directory where the files are

#gets the files
for file in glob.glob('*fits'):
    #pdb.set_trace()
    value = fits.getheader(file)['SHUTRPOS']
    value = value.lower().strip()
    #print('{} -> {}'.format(file,value))

    if value.lower().strip() in 'a':
    	fileName= os.path.join(shutter_a, file)
    	curentName= os.path.join(base_path, file)
    	os.rename(curentName, fileName)
    	print('Moving {} to shutter_a/'.format(file))
    	#shutil.move('/grp/hst/wfc3b/post_flash_cal/shutter_a/'+file, 
    	 #           '/grp/hst/wfc3b/post_flash_cal/'+file)
    elif value.lower().strip() in 'b':
    	curentName= os.path.join(base_path, file)
    	fileName=os.path.join(shutter_b, file)
    	os.rename(curentName, fileName)
    	print('Moving {} to shutter_b/'.format(file))
    	#shutil.move(/grp/hst/wfc3b/post_flash_cal/shutter_b/fitsName, /grp/hst/wfc3b/post_flash_cal/fitsName)

    else:
        print('Shutter is not specified for {}'.format(file))
    


os.chdir(current)
#returns you to your origanal directory

