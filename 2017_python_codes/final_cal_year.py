#Filename: final_cal.py
#Description: This makes the final adjustments to the pf file to make it a referance file.
#Date: August 22, 2016
#Author: Heather Kurtz



from astropy.io import fits
import numpy as np
import csv
import os
import glob



def final_cal(shutterB,currentfileB,errorB,a_flag):


	
	stats=[]
	
	# reads in the fits files B
	hdulistB= fits.open(shutterB)
	dataB1=hdulistB[0].data
	dataB4=hdulistB[1].data
	
	hdulistB.close
	
	hdulistBe= fits.open(errorB)
	errorB1=hdulistBe[0].data
	errorB4=hdulistBe[1].data
	
	hdulistBe.close
	
	# read in the header of the current reference file
	
	hdulistc= fits.open(currentfileB)
	Bhdr=hdulistc[0].header
	Bhdr1=hdulistc[1].header
	Bhdr2=hdulistc[2].header
	Bhdr3=hdulistc[3].header
	Bhdr4=hdulistc[4].header
	Bhdr5=hdulistc[5].header
	Bhdr6=hdulistc[6].header
	
	#divids that data by the exposure time, multiplies by the scale factor and multiplies by the gain
	if a_flag == 'B':
	
		dataB1_low=((dataB1/100.00)*0.03623)*1.56
		dataB4_low=((dataB4/100.00)*0.03623)*1.56
		
		errorB1_cor=((errorB1*0.03623)/100.0)
		errorB4_cor=((errorB4*0.03623)/100.0)

	elif a_flag == 'A':
		dataB1_low=((dataB1/100.00)*0.03639)*1.56
		dataB4_low=((dataB4/100.00)*0.03639)*1.56

		errorB1_cor=((errorB1*0.03639)/100.0)
		errorB4_cor=((errorB4*0.03639)/100.0)
	
	#creates the other extentions needed for the fits referance file
	dq1=np.zeros((2070,4206))
	dq2=np.zeros((2070,4206))
	
	#adds zeros in rows to data to obtain the correct form
	c=np.zeros((19, 4096))
	dataB1_19=np.concatenate((dataB1_low,c), axis=0)
	
	dataB4_19=np.concatenate((c,dataB4_low), axis=0)
	
	errorB1_19=np.concatenate((errorB1_cor,c), axis=0)
	
	errorB4_19=np.concatenate((c,errorB4_cor), axis=0)
	
	#this makes the edges have 25 coloms of zeros
	d=np.zeros((2070,25 ))
	
	dataB4_f25=np.concatenate((dataB4_19,d), axis=1)
	dataB1_f25=np.concatenate((dataB1_19,d), axis=1)
	
	dataB4_e25=np.concatenate((d,dataB4_f25), axis=1)
	dataB1_e25=np.concatenate((d,dataB1_f25), axis=1)
	
	errorB4_f25=np.concatenate((errorB4_19,d), axis=1)
	errorB1_f25=np.concatenate((errorB1_19,d), axis=1)
	
	errorB4_e25=np.concatenate((d,errorB4_f25), axis=1)
	errorB1_e25=np.concatenate((d,errorB1_f25), axis=1)
	
	#add 60 coloms of zeros in the middle
	idx=[]
	for i in range(60):
	    idx.append(2073)
	
	#idx= list(range(2073,2133))

	dataB4_all=np.insert(dataB4_e25, idx, 0, axis=1)
	dataB1_all=np.insert(dataB1_e25, idx, 0, axis=1)
	
	errorB4_all=np.insert(errorB4_e25, idx, 0, axis=1)
	errorB1_all=np.insert(errorB1_e25, idx, 0, axis=1)
	
	    
	print('dataB4')
	print('Min:', np.min(dataB4_all))
	print('Max:', np.max(dataB4_all))
	print('Mean', np.mean(dataB4_all))
	print('Stdev:', np.std(dataB4_all))      
	    
	print('dataB1')
	print('Min:', np.min(dataB1_all))
	print('Max:', np.max(dataB1_all))
	print('Mean', np.mean(dataB1_all))
	print('Stdev:', np.std(dataB1_all)) 
	
	#write out the new files!!!!!
	#each individual extention
	#fits.writeto(outfile1,dataA1_all )
	outfile1= 'B4m.fits'
	#fits.writeto(outfile1,dataB4_all )
	outfile1= 'B1m.fits'
	#fits.writeto(outfile1,dataB1_all )
	
	
	#shutter B all together with all 7 extensions including the header
	cal_out_file = shutterB[:-5] + 'final_cal.fits' 
	prihdu = fits.PrimaryHDU(header=Bhdr)
	
	single_extension1=fits.ImageHDU(data=dataB1_all.astype(np.float32),header=Bhdr1)
	single_extension2=fits.ImageHDU(data=errorB1_all.astype(np.float32),header=Bhdr2)
	single_extension3=fits.ImageHDU(data=dq1.astype(np.float32),header=Bhdr3)
	
	single_extension4=fits.ImageHDU(data=dataB4_all.astype(np.float32),header=Bhdr4)
	single_extension5=fits.ImageHDU(data=errorB4_all.astype(np.float32),header=Bhdr5)
	single_extension6=fits.ImageHDU(data=dq2.astype(np.float32),header=Bhdr6)
	
	all_extensions=[prihdu,single_extension1, single_extension2, single_extension3,single_extension4,single_extension5,single_extension6]
	myhdulist=fits.HDUList(all_extensions)
	myhdulist[0]._bitpix = 16
	myhdulist.writeto(cal_out_file, overwrite=True)
	
	hdulistc.close


def main():
	
	#gets current directory
	current = os.getcwd()
	
	base_path = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis_2/shutter_b/med_100/'
	os.chdir(base_path)
	
	files= glob.glob('*median.fits')
	print(len(files))
	currentfileB='wc52031pi_fls.fits'
	flag = 'B'

	for f in files:
		errorB= f[:-11] + 'error.fits'
		final_cal(f,currentfileB,errorB,flag)


	base_path2 = '/grp/hst/wfc3v/hkurtz/pf_2020/analysis_2/shutter_a/med_100/'
	os.chdir(base_path2)
	
	files= glob.glob('*median.fits')
	print(len(files))
	currentfileA='wc52031oi_fls.fits'
	flag = 'A'

	for fl in files:
		errorA= fl[:-11] + 'error.fits'
		final_cal(fl,currentfileA,errorA,flag)

	os.chdir(current)
	
main()


















