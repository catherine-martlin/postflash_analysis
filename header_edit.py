#Filename: header_edit.py
#Description: This makes the final adjustments PF header
#Author: Heather Kurtz
#Date: August 23, 2016


from astropy.io import fits
import os

current = os.getcwd()
base_path = '/Users/hkurtz/Desktop/post_flash/combine_im'
os.chdir(base_path)

shutter_a='cal_image_A.fits'
shutter_b='cal_image_B.fits'

#open files in update mode
hdulistA= fits.open(shutter_a,mode='update')
hdulistB= fits.open(shutter_b,mode='update')


hdr_a=hdulistA[0].header
hdr_b=hdulistB[0].header

#change keyword values for the new file
hdr_a["DATE"]='2017-03-03T00:00:00'
hdr_b["DATE"]='2017-03-03T00:00:00'

hdr_a['FILENAME']='cal_image_A.fits'
hdr_b['FILENAME']='cal_image_B.fits'

hdr_a['ROOTNAME']='cal_image_A'
hdr_b['ROOTNAME']='cal_image_B'

hdr_a['CAL_VER']='calwf 3.3'
hdr_b['CAL_VER']='calwf 3.3'

hdr_a['BIASCORR']='COMPLETE'
hdr_b['BIASCORR']='COMPLETE'

hdr_a['EXPSCORR']='COMPLETE'
hdr_b['EXPSCORR']='COMPLETE'

hdr_a['CRCORR']='COMPLETE'
hdr_b['CRCORR']='COMPLETE'

hdr_a['DARKCORR']='COMPLETE'
hdr_b['DARKCORR']='COMPLETE'

hdr_a['BPIXTAB']='iref$u5d2012li_bpx.fits'
hdr_b['BPIXTAB']='iref$u5d2012li_bpx.fits'

hdr_a['ATODTAB']='N/A'
hdr_b['ATODTAB']='N/A'

hdr_a['DARKFILE']='iref$zcv15508i_drk.fits'
hdr_b['DARKFILE']='iref$zcv15508i_drk.fits'

hdr_a['PFLTFILE']='N/A'
hdr_b['PFLTFILE']='N/A'

hdr_a['GRAPHTAB']='mtab$yah1742qm_tmg.fits'
hdr_b['GRAPHTAB']='mtab$yah1742qm_tmg.fits'

hdr_a['COMPTAB']='mtab$yah1742rm_tmc.fits'
hdr_b['COMPTAB']='mtab$yah1742rm_tmc.fits'

hdr_a['IDCTAB']='iref$yas1621ai_idc.fits'
hdr_b['IDCTAB']='iref$yas1621ai_idc.fits'

hdr_a['SCALENSE']=0.0
hdr_b['SCALENSE']=0.0

hdr_a['INITGUES']='none'
hdr_b['INITGUES']='none'

hdr_a['CRSIGMAS']=0.0
hdr_b['CRSIGMAS']=0.0

hdr_a['CRRADIUS']=0.0
hdr_b['CRRADIUS']=0.0

hdr_a['CRTHRESH']=0.0
hdr_b['CRTHRESH']=0.0

hdr_a['BADINPDQ']=0
hdr_b['BADINPDQ']=0

hdr_a['REJ_RATE']=0.0
hdr_b['REJ_RATE']=0.0

hdr_a['PEDIGREE']='INFLIGHT 27/08/2012 19/07/2016'
hdr_b['PEDIGREE']='INFLIGHT 27/08/2012 19/07/2016'

#hdr_a['USEAFTER']='Jan 01 2012 00:00:00'
#hdr_b['USEAFTER']='Jan 01 2012 00:00:00'

#delets old comments
del hdr_a['COMMENT'] 
del hdr_b['COMMENT']

#adds my new comments
hdr_a['COMMENT']='Reference file generated by H.Kurtz. March 3, 2017'
hdr_b['COMMENT']='Reference file generated by H.Kurtz. March 3, 2017'

hdr_a['SHUTRPOS']='A'
hdr_b['SHUTRPOS']='B'




#deleat old history so as not to confuse the two files
del hdr_a['HISTORY'] 
del hdr_b['HISTORY']





#creat the new history for my data based on the information for the old history file
hdr_a['HISTORY']='This is the version 2.0 post-flash reference for low current,'
hdr_a['HISTORY']='shutter A, and unbinned pixels.'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='It was generated as follows:'

#hdr_a['HISTORY']='First a high signal-to-noise ratio 100 sec medium current shutter A'
#hdr_a['HISTORY']='flash image was generated using 111 images from proposals 13560, 14372,' 
#hdr_a['HISTORY']='14006, 13069, 13078:'
#hdr_a['HISTORY']=''
#hdr_a['HISTORY']='These images were combined to make a low current referance file.'
#hdr_a['HISTORY']='Then using the anayisis done for the previouse referance file,'
#hdr_a['HISTORY']='a scale factor of 28.96 was used to convert from the low current'
#hdr_a['HISTORY']='to the medium current.'
#hdr_a['HISTORY']=''
#hdr_a['HISTORY']='Finally the reference files for LOW current were multiplied by this'
#hdr_a['HISTORY']='scaling factor and re-delivered as refernence files for MEDIUM'
#hdr_a['HISTORY']='current.'hdr_a['HISTORY']=''
#hdr_a['HISTORY']=''
hdr_a['HISTORY']='First a high signal-to-noise ratio 100 sec medium current shutter A'
hdr_a['HISTORY']='flash image was generated using 111 images from proposals 13560, 14372,' 
hdr_a['HISTORY']='14006, 13069, 13078:'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='The images were calibrated in CALWF3 version 3.3. Populating the'
hdr_a['HISTORY']='DQ array. Which were then masked before median combining the '
hdr_a['HISTORY']='images.'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='Finally the 100 sec medium current flash image was scaled to the'
hdr_a['HISTORY']='appropriate brightness for a 1 sec low current flash image.  Using'
hdr_a['HISTORY']='the scaling factor from the previous post flash reference file of'
hdr_a['HISTORY']='0.03639'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='The effective flash in 1 sec was computed for the medium current'
hdr_a['HISTORY']='images by dividing by the appropriate flash duration 100 sec.'
hdr_a['HISTORY']='the 1 sec medium flash image was then multiplied by 0.03639 to'
hdr_a['HISTORY']='get the low flash 1 sec low current reference file.'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='CALWF3 expects the post flash image to be in units of electrons,'
hdr_a['HISTORY']='so the reference file values have been scaled by 1.56'
hdr_a['HISTORY']=''
hdr_a['HISTORY']='H. Kurtz -- 2017march3'

#hdr_a['HISTORY']='Files used to make this:'



hdr_b['HISTORY']='This is the version 2.0 post-flash reference for low current,'
hdr_b['HISTORY']='shutter B, and unbinned pixels.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='It was generated as follows:'
#hdr_b['HISTORY']=''
#hdr_b['HISTORY']='First a high signal-to-noise ratio 100 sec medium current shutter B'
#hdr_b['HISTORY']='flash image was generated using 112 images from proposals: 13568, 14372,' 
#hdr_b['HISTORY']='13069, 14006, 13078, 13560.'
#hdr_b['HISTORY']=''
#hdr_b['HISTORY']='These images were combined to make a low current referance file.'
#hdr_b['HISTORY']='Then using the anayisis done for the previouse referance file,'
#hdr_b['HISTORY']='a scale factor of 28.96 was used to convert from the low current'
#hdr_b['HISTORY']='to the medium current.'
#hdr_b['HISTORY']=''
#hdr_b['HISTORY']='Finally the reference files for LOW current were multiplied by this'  
#hdr_b['HISTORY']='scaling factor and re-delivered as refernence files for MEDIUM'        
#hdr_b['HISTORY']='current.'
#hdr_b['HISTORY']=''
hdr_b['HISTORY']='First a high signal-to-noise ratio 100 sec medium current shutter B'
hdr_b['HISTORY']='flash image was generated using 112 images from proposals: 13568, 14372,' 
hdr_b['HISTORY']='13069, 14006, 13078, 13560.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='The images were recalibrated in CALWF3 version 3.3. Populating the'
hdr_b['HISTORY']='DQ array with cosmic ray positions. Which were then masked'
hdr_b['HISTORY']='before median combining the images.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='Finally the 100 sec medium current flash image was scaled to the'
hdr_b['HISTORY']='appropriate brightness for a 1 sec low current flash image.  Using'
hdr_b['HISTORY']='the scaling factor from the previous post flash reference file of'
hdr_b['HISTORY']='0.03639'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='The effective flash in 1 sec was computed for the medium current'
hdr_b['HISTORY']='images by dividing by the appropriate flash duration 100 sec.'
hdr_b['HISTORY']='the 1 sec medium flash image was then multiplied by 0.03639 to'
hdr_b['HISTORY']='get the low flash 1 sec low current reference file.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='CALWF3 expects the post flash image to be in units of electons,'
hdr_b['HISTORY']=' so the reference file values have been scaled by 1.56'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='112 images from proposals 13568, 13560, 14372, 14006, 13069, 13078: were'
hdr_b['HISTORY']='recalibrated in CALWF3 version 3.3, the median image created, scaled by'
hdr_b['HISTORY']='0.03639, and multiplied by the gain.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='This file is different from the previous post flash file because it'
hdr_b['HISTORY']='uses newer data and more files are used in the creation of this file.'
hdr_b['HISTORY']='Statistically the two are similar. This file is on average 0.017 fainter'
hdr_b['HISTORY']='than the previouse file.'
hdr_b['HISTORY']=''
hdr_b['HISTORY']='H. Kurtz -- 2017march3'
#hdr_b['HISTORY']='Files used to make this:'







hdulistA.close()
hdulistB.close()








