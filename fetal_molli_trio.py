# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 15:19:28 2016

@author: sportnoy
"""

import numpy as np
import blood_tools
import os
import pylab as plt
from uncertainties import ufloat
from misc import fitting

parent_folder= os.getcwd()

scans_to_process=['/3T_data/10mar2016_trio/se_mc_comp_full_spoil 48ms ETL8_68']


#scans_to_process=['27feb2016_trio/T1map_TRUFI_40DEG_74',
#                  '10mar2016_trio/T1map_TRUFI_40DEG shorter_69',
#                  '30mar2016_trio/T1map_TRUFI_40DEG shorter_78',
#                  '23apr2016_trio/T1map_TRUFI_40DEG shorter_82',
#                  '27apr2016_trio/T1map_TRUFI_40DEG shorter_82']
                  
MCHCs=np.array([33.1,34.1,33.5,33.0,33.8,34.4])      

TI_full=np.arange(0,20000,1)                  

T1_dictionary={}
for jj, scan in enumerate(scans_to_process): 
    plt.figure()
    folder = parent_folder + scan
    images,TI=blood_tools.read_dicoms(folder,['InversionTime'])
    ti=np.array([element['InversionTime'] for element in TI])        
    
    rois=blood_tools.load_ROIs(folder + '/rois')
    bga_data=blood_tools.read_bga(os.path.split(folder)[0]+'/sO2_Hct')  
    sO2s=bga_data[0]
    ctHbs=bga_data[1]
    metHbs=bga_data[2]
    Hcts=ctHbs/MCHCs[jj]     
    patient_number=jj                       
    
    for mm, roi in enumerate(rois):  
        sig=blood_tools.calc_ROI_mean(roi,images)
        inversion_recovery=fitting.model('abs(M0*(1-2*aa*exp(-x/T1)))',
                                     {'M0':sig.max(),'aa':1,'T1':1000})  
        #Initial T1 guess based on null time    
        T1_guess=-ti[np.where(sig==sig.min())]/np.log(0.5) 
        if(np.size(T1_guess)>1):
            inversion_recovery['T1']=T1_guess[0]
        else:
            inversion_recovery['T1']=T1_guess  
       # inversion_recovery['T1']=2000                             
        inversion_recovery.fit(ti,sig)
        
        print inversion_recovery['aa']
        
        plt.subplot(4,4,mm+1)
        plt.plot(ti, sig,'o')
        plt.plot(TI_full,inversion_recovery(TI_full))          
        T1_corr=inversion_recovery['T1'].value*(2*inversion_recovery['aa'].value-1)  
        

        T1_full=ufloat(T1_corr,10)
        
        key=(sO2s[mm], Hcts[mm], metHbs[mm], patient_number)
        
        T1_dictionary[key]=T1_full 
        
plt.figure()        
R1s=[]
T1s=[]
T1_errs=[]
R1_errs=[]
sO2s=[]
Hcts=[] 
 
for key in T1_dictionary:
    if (key[3]==5):
        T1s.append(T1_dictionary[key].nominal_value)
        T1_errs.append(T1_dictionary[key].std_dev)
        Hcts.append(key[1])
        R1=1000/T1_dictionary[key]
        R1s.append(R1.nominal_value)
        R1_errs.append(R1.std_dev)
        sO2s.append(key[0])
#        
#plt.plot(Hcts,R1s,'o')

import matplotlib
plt.scatter(Hcts, T1s, c=sO2s, s=125, vmin=0,vmax=1.0, marker='h',edgecolor='black',cmap=matplotlib.cm.cool) 



#Check chi-square surface
#pt=0
#specimen=6
#
#folder = parent_folder + scans_to_process[pt]
#images,TI=blood_tools.read_dicoms(folder,['InversionTime'])
#ti=np.array([element['InversionTime'] for element in TI])        
#rois=blood_tools.load_ROIs(folder + '/rois')
#roi=rois[specimen]
#
#sig=blood_tools.calc_ROI_mean(roi,images)
#inversion_recovery=fitting.model('abs(M0*(1-2*aa*exp(-x/T1)))',
#                             {'M0':sig.max(),'aa':1,'T1':1000}) 
#                             
#
##Initial T1 guess based on null time    
##T1_guess=-ti[np.where(sig==sig.min())]/np.log(0.5) 
##if(np.size(T1_guess)>1):
##    inversion_recovery['T1']=T1_guess[0]
##else:
##    inversion_recovery['T1']=T1_guess  
#inversion_recovery['T1']=1400.0
##inversion_recovery['T1'].freeze()
#inversion_recovery.fit(ti,sig)
#chi_min=inversion_recovery.chi2()
#plt.figure()
#plt.plot(ti, sig,'o')
#plt.plot(TI_full,inversion_recovery(TI_full))    
#
#T1_corr=inversion_recovery['T1'].value*(2*inversion_recovery['aa'].value-1)  
#T1s=np.arange(800,2000,1)
##aas=np.arange(0.1,1.9,0.01)
#M0s=np.arange(100,600,2)
##chi_surf=inversion_recovery.contours(inversion_recovery['T1'],inversion_recovery['aa'],T1s,aas) 
#chi_surf=inversion_recovery.contours(inversion_recovery['T1'],inversion_recovery['M0'],T1s,M0s)       