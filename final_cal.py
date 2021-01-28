#Filename: final_cal.py
#Description: This makes the final adjustments to the pf file to make it a referance file.
#Date: August 22, 2016
#Author: Heather Kurtz



from astropy.io import fits
import numpy as np
import csv
import os

#gets current directory
current = os.getcwd()

shutterA='A_median_pf_100.fits'
shutterB='B_median_pf_100.fits'
currentfileB='wc52031pi_fls.fits'
currentfileA='wc52031oi_fls.fits'
errorA='A_error_pf_100.fits'
errorB='B_error_pf_100.fits'

stats=[]

# reads in the fits files A
hdulistA= fits.open(shutterA)
dataA1=hdulistA[0].data
dataA4=hdulistA[1].data

hdulistA.close

hdulistAe= fits.open(errorA)
errorA1=hdulistAe[0].data
errorA4=hdulistAe[1].data

hdulistAe.close

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
hdulistc= fits.open(currentfileA)
Ahdr=hdulistc[0].header
Ahdr1=hdulistc[1].header
Ahdr2=hdulistc[2].header
Ahdr3=hdulistc[3].header
Ahdr4=hdulistc[4].header
Ahdr5=hdulistc[5].header
Ahdr6=hdulistc[6].header

hdulistc= fits.open(currentfileB)
Bhdr=hdulistc[0].header
Bhdr1=hdulistc[1].header
Bhdr2=hdulistc[2].header
Bhdr3=hdulistc[3].header
Bhdr4=hdulistc[4].header
Bhdr5=hdulistc[5].header
Bhdr6=hdulistc[6].header

#divids that data by the exposure time, multiplies by the scale factor and multiplies by the gain
dataA1_low=((dataA1/100.00)*0.03639)*1.56
dataA4_low=((dataA4/100.00)*0.03639)*1.56

dataB1_low=((dataB1/100.00)*0.03639)*1.56
dataB4_low=((dataB4/100.00)*0.03639)*1.56

errorA1_cor=((errorA1*0.03639)/100.0)
errorA4_cor=((errorA4*0.03639)/100.0)
errorB1_cor=((errorB1*0.03639)/100.0)
errorB4_cor=((errorB4*0.03639)/100.0)

#creates the other extentions needed for the fits referance file
dq1=np.zeros((2070,4206))
dq2=np.zeros((2070,4206))

#adds zeros in rows to data to obtain the correct form
c=np.zeros((19, 4096))
dataA1_19=np.concatenate((dataA1_low, c), axis=0)
dataB1_19=np.concatenate((dataB1_low,c), axis=0)

dataA4_19=np.concatenate((c,dataA4_low), axis=0)
dataB4_19=np.concatenate((c,dataB4_low), axis=0)

errorA1_19=np.concatenate((errorA1_cor, c), axis=0)
errorB1_19=np.concatenate((errorB1_cor,c), axis=0)

errorA4_19=np.concatenate((c,errorA4_cor), axis=0)
errorB4_19=np.concatenate((c,errorB4_cor), axis=0)

#this makes the edges have 25 coloms of zeros
d=np.zeros((2070,25 ))

dataA4_f25=np.concatenate((dataA4_19,d), axis=1)
dataB4_f25=np.concatenate((dataB4_19,d), axis=1)
dataA1_f25=np.concatenate((dataA1_19,d), axis=1)
dataB1_f25=np.concatenate((dataB1_19,d), axis=1)

dataA4_e25=np.concatenate((d,dataA4_f25), axis=1)
dataB4_e25=np.concatenate((d,dataB4_f25), axis=1)
dataA1_e25=np.concatenate((d,dataA1_f25), axis=1)
dataB1_e25=np.concatenate((d,dataB1_f25), axis=1)

errorA4_f25=np.concatenate((errorA4_19,d), axis=1)
errorB4_f25=np.concatenate((errorB4_19,d), axis=1)
errorA1_f25=np.concatenate((errorA1_19,d), axis=1)
errorB1_f25=np.concatenate((errorB1_19,d), axis=1)

errorA4_e25=np.concatenate((d,errorA4_f25), axis=1)
errorB4_e25=np.concatenate((d,errorB4_f25), axis=1)
errorA1_e25=np.concatenate((d,errorA1_f25), axis=1)
errorB1_e25=np.concatenate((d,errorB1_f25), axis=1)

#add 60 coloms of zeros in the middle
idx=[]
for i in range(60):
    idx.append(2073)

#idx= list(range(2073,2133))

dataA4_all=np.insert(dataA4_e25, idx, 0, axis=1)
dataA1_all=np.insert(dataA1_e25, idx, 0, axis=1)
dataB4_all=np.insert(dataB4_e25, idx, 0, axis=1)
dataB1_all=np.insert(dataB1_e25, idx, 0, axis=1)

errorA4_all=np.insert(errorA4_e25, idx, 0, axis=1)
errorA1_all=np.insert(errorA1_e25, idx, 0, axis=1)
errorB4_all=np.insert(errorB4_e25, idx, 0, axis=1)
errorB1_all=np.insert(errorB1_e25, idx, 0, axis=1)

print('dataA4')
print('Min:', np.min(dataA4_all))
print('Max:', np.max(dataA4_all))
print('Mean', np.mean(dataA4_all))
print('Stdev:', np.std(dataA4_all))      
    
print('dataA1')
print('Min:', np.min(dataA1_all))
print('Max:', np.max(dataA1_all))
print('Mean', np.mean(dataA1_all))
print('Stdev:', np.std(dataA1_all)) 
    
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
outfile1= 'A4m.fits'
#fits.writeto(outfile1,dataA4_all )
outfile1= 'A1m.fits'
#fits.writeto(outfile1,dataA1_all )
outfile1= 'B4m.fits'
#fits.writeto(outfile1,dataB4_all )
outfile1= 'B1m.fits'
#fits.writeto(outfile1,dataB1_all )

#shutter A all together with all 7 extensions including the header
prihdu = fits.PrimaryHDU(header=Ahdr)

single_extension1=fits.ImageHDU(data=dataA1_all.astype(np.float32),header=Ahdr1)
single_extension2=fits.ImageHDU(data=errorA1_all.astype(np.float32),header=Ahdr2)
single_extension3=fits.ImageHDU(data=dq1.astype(np.float32),header=Ahdr3)

single_extension4=fits.ImageHDU(data=dataA4_all.astype(np.float32),header=Ahdr4)
single_extension5=fits.ImageHDU(data=errorA4_all.astype(np.float32),header=Ahdr5)
single_extension6=fits.ImageHDU(data=dq2.astype(np.float32),header=Ahdr6)

all_extensions=[prihdu,single_extension1, single_extension2, single_extension3,single_extension4,single_extension5,single_extension6]
myhdulist=fits.HDUList(all_extensions)
myhdulist[0]._bitpix = 16
myhdulist.writeto('cal_image_A.fits', clobber=True)

#shutter B all together with all 7 extensions including the header
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
myhdulist.writeto('cal_image_B.fits', clobber=True)

hdulistc.close