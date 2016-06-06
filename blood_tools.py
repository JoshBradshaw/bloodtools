# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 14:35:19 2015

@author: sportnoy
Sharon's Tools for blood data analysis
"""
import numpy as np
import cPickle
import os
import fitting
import pylab as plt
import dicom
import sys
import matplotlib
import xlwt
import itertools
import re

def calc_ROI_mean(roi,images):
    """Calculate the mean of an ROI across a list of images"""
    if isinstance(images,list):
        mean=np.zeros(len(images))
        for jj in np.arange(len(images)):
            mean[jj]=images[jj][roi.get_indices()].mean()
    else:
        mean=images[roi.get_indices()].mean()
    return mean

def calc_ROI_serr(ROI,images):
    """Calculate the mean of an ROI across a list of images"""  
    if isinstance(images,list):
        serr=np.zeros(len(images))
        for jj in np.arange(len(images)):
            serr[jj]=images[jj][ROI.get_indices()].std()/np.sqrt(len(ROI.get_indices()[0]))           
    else:
        serr=images[ROI.get_indices()].std()/np.sqrt(len(ROI.get_indices()[0]))       
    return serr    

def calc_ROI_std(ROI,images):
    """Calculate the mean of an ROI across a list of images"""  
    if isinstance(images,list):
        std=np.zeros(len(images))
        for jj in np.arange(len(images)):
            std[jj]=images[jj][ROI.get_indices()].std()          
    else:
        std=images[ROI.get_indices()].std()
    return std    
    
def load_ROIs(filename):
    """Load all of the ROIs within a file"""        
    f=open(filename,'r')
    roi_list=[]
    while True:
        try:
            roi_list.append(cPickle.load(f))
        except EOFError:
            break
    f.close()
    return roi_list

def read_dicoms(foldername, attributes=[]):
    """Read in all the dicom files in a folder and any necessary attributes"""   
    file_list= [fn for fn in os.listdir(foldername) if os.path.isfile(os.path.join(foldername, fn))]
    file_list.sort()
    image_list=[]
    dicom_list = []
    attribute_list=[]
    for file_name in file_list:    
        attribute_dict={}    
        try:
            dicom_obj=dicom.read_file(os.path.join(foldername, file_name))
            dicom_list.append(dicom_obj)
        except dicom.filereader.InvalidDicomError:
            print "Error, invalid dicom file: {}".format(os.path.join(foldername, file_name))
            continue            
            
        image_list.append(dicom_obj.pixel_array)
        attribute_dict['filename'] = file_name
        for item in attributes:
            attribute_dict[item]=getattr(dicom_obj,item,0)
        attribute_list.append(attribute_dict)
    return image_list, attribute_list, dicom_list     

def img_roi_signal(folders_to_process,attributes=[]):
    """Function to return an array with the mean and std of the signal
    within each ROI (defined in a file called 'rois' within each folder) 
    for each image within each folder. If attribute is not null function
    will also provide the value of the attribute (e.g. echo time) for 
    each image"""
    folder_image_list=[]
    roi_list=[]
    attribute_lists=[]
    max_images=0

    for ii,folder in enumerate(folders_to_process):
        image_list,attribute_list=read_dicoms(folder,attributes)
        folder_image_list.append(image_list)
        attribute_list=[entry[attributes[0]] for entry in attribute_list] 
        attribute_lists.append(attribute_list)
        roifile=folder+'/rois'
        try:
            rois=load_ROIs(roifile)   
        except IOError:
            error_string=str(folder) + ' does not have an ROI file.'                
            print error_string
            sys.exit(1)
        roi_list.append(rois)
        max_images=np.max([max_images,len(image_list)])
    
    mean_signal_mat=np.zeros([len(folders_to_process),len(rois),max_images])
    serr_signal_mat=np.zeros([len(folders_to_process),len(rois),max_images])
    
    for ii,images in enumerate(folder_image_list):      
        for jj,roi in enumerate(roi_list[ii]):
            mean_signal_mat[ii,jj,0:len(images)]=calc_ROI_mean(roi,images)
            serr_signal_mat[ii,jj,0:len(images)]=calc_ROI_serr(roi,images)
    return roi_list, folder_image_list, attribute_lists, mean_signal_mat, serr_signal_mat        

def SE_fit_new(te, signal, mean_noise=0, noise_floor='n', noise_factor=2):
    """Fit spin echo data.  Will not include points with intensity < noise_factor*mean_noise in the fit"""

    serr_signal=np.ones_like(signal)
        
    if noise_floor=='n':    
        spin_echo=fitting.model('M0*exp(-x/T2)',{'M0':signal.max(),'T2':100})
        spin_echo.fitpars['method']='leastsq'    
        #print mean_noise
        temp=np.where(signal<noise_factor*mean_noise)[0]
        if temp.size:
            stop_index=temp[0]
        else:
            stop_index=len(signal)
        print stop_index
        
        if stop_index<2:
            #print 'T2 is too short, need to fit for noise floor'
            print 'T2 too short, data point needs to be excluded'
            spin_echo=fitting.model('M0*exp(-x/T2)+a',{'M0':signal.max(),'T2':-1,'a':signal.min()}) 
            #spin_echo.fitpars['method']='leastsq'             
            #spin_echo.fit(te,signal)      
        else:
            spin_echo.fit(te[0:stop_index],signal[0:stop_index],serr_signal[0:stop_index])         
          
    else:   #noise floor=='y'
        print 'explicitly fitting for noise floor'
        spin_echo=fitting.model('M0*exp(-x/T2)+a',{'M0':signal.max(),'T2':100,'a':signal.min()}) 
#        spin_echo['a']=mean_noise
#        spin_echo['a'].freeze()        
        spin_echo.fit(te,signal) 
        #ci=spin_echo.conf(sigma=1) 
    return spin_echo

def T2_cpmg_process(folder_to_process,plot='y'):
    """Given a folder of images will process cpmg data and return
    fitted T2 values and associated uncertainties"""
    data=img_roi_signal([folder_to_process],['EchoTime'])
    rois=data[0][0]
    TEs=data[2][0]
    mean_signal_mat=data[3]
    serr_signal_mat=data[4]
        
    if plot=='y':
        plt.figure()
    spin_echo_fits=[]
    for jj in np.arange(len(rois)-2):
        mean_sig=mean_signal_mat[0,jj,:]
        #serr_sig=serr_signal_mat[0,jj,:]
        mean_noise=np.mean(mean_signal_mat[0,-2,:])
        try:
            spin_echo_fit = SE_fit_new( np.array(TEs[0:]), mean_sig[0:], mean_noise, 'n' )
            if plot=='y':
                TE_full=np.arange(0,400,1)
                plt.subplot(4,4,jj+1)
                plt.plot(np.array(TEs[0:]), mean_sig[0:],'o')
                plt.plot(TE_full,spin_echo_fit(TE_full))
            
            spin_echo_fits.append(spin_echo_fit)
        except RuntimeError: 
            print 'RuntimeError'
            spin_echo=fitting.model('M0*exp(-x/T2)+a',{'M0':0,'T2':0,'a':0}) 
            spin_echo_fits.append(spin_echo)
    return spin_echo_fits   
    
def read_bga(filename):
    """Function to read in blood gas analyzer data from a pickle"""     
    f=open(filename,'r')
    data=[]
    while True:
        try:
            data.append(cPickle.load(f))
        except EOFError:
            break
    f.close()
    sO2=data[0]
    Hct=data[1]
    metHb=data[2]
    return sO2, Hct, metHb
    
def calc_R_err(T,T_err):
    """Calculates error in R from T and T error (using 
    propagation of error)"""
    R_err=1000/T**2*T_err  
    return R_err   

def IR_fit(ti,signal,serr_signal=[]):
    """Fit inversion recovery data"""
    if serr_signal==[]:
        serr_signal=np.ones_like(signal) 
    inversion_recovery=fitting.model('abs(M0*(1-2*aa*exp(-x/T1)))',
                                     {'M0':signal.max(),'aa':1,'T1':1000})                                 
    #Initial T1 guess based on null time    
    T1_guess=-ti[np.where(signal==signal.min())]/np.log(0.5) 
    if(np.size(T1_guess)>1):
        inversion_recovery['T1']=T1_guess[0]
    else:
        inversion_recovery['T1']=T1_guess                                
    inversion_recovery.fit(np.array(ti),signal,serr_signal)  
    chi2=inversion_recovery.chi2()
    print chi2
    return inversion_recovery
    
def IR_process(folders_to_process):
    """Calculates T1s given a list of folders, each corresponding to an individual TI"""
    data=img_roi_signal(folders_to_process,['InversionTime'])
    TIs=np.concatenate(data[2])
    mean_signal_mat=data[3]
    serr_signal_mat=data[4]
    T1s=[]
    aas=[]
    plt.figure()
    rois=load_ROIs(folders_to_process[0]+'/rois')
    for kk,roi in enumerate(rois[0:-2]):
        print kk
        plt.subplot(4,4,kk+1)
        #chi2, chi_surf, M0, T1, T1_err, aa = IR_fit(TIs, mean_signal_mat[:,kk,0], serr_signal_mat[:,kk,0],'n')
        fit = IR_fit(TIs, mean_signal_mat[:,kk,0])
        plt.plot(TIs, mean_signal_mat[:,kk,0],'o')
        TI_full=np.arange(1,8000,5)           
        plt.plot(TI_full, fit(TI_full))            
        T1s.append(fit['T1'].value)
        aas.append(fit['aa'].value)
    return T1s, aas     
    
def T2_cpmg_bootstrap(folder,N=1000,plot='y'):
    """Given a folder of images will process cpmg data and return
    fitted T2 values and associated uncertainties.  Uncertainties 
    are obtained by bootstrapping pixels within each ROI"""  
    roifile=folder +'/rois'
    try:
        rois=load_ROIs(roifile)   
    except IOError:
        error_string=str(folder) + ' does not have an ROI file.'                
        print error_string
        sys.exit(1)
    image_list,TEs=read_dicoms(folder,['EchoTime'])
    tes=np.array([entry['EchoTime'] for entry in TEs])
    mean_noise=np.mean([img[rois[-2].get_indices()].mean() for img in image_list])
    T2s=[]
    T2_errs=[]
    if plot=='y':
        plt.figure()
    for kk,roi in enumerate(rois[0:-2]):
        print kk
        npix=len(roi.get_indices()[0])
        pixels=roi.get_indices()
        sig=[img[pixels].mean() for img in image_list]
        noise=[img[rois[-2].get_indices()].mean() for img in image_list]
        spin_echo_fit = SE_fit_new( tes, np.array(sig), mean_noise, 'n' )
        if plot=='y':
            TE_full=np.arange(0,2000,1)
            plt.subplot(4,4,kk+1)
            plt.plot(tes, sig, 'o')
            plt.plot(tes,noise,'ok')
            plt.semilogy(TE_full,spin_echo_fit(TE_full))        
        T2s.append(spin_echo_fit['T2'].value)        
        T2bs=[]
        for nn in np.arange(N):
            ind=np.random.randint(npix,size=npix)
            sig=[img[pixels[0][ind],pixels[1][ind]].mean() for img in image_list]
            try:
                spin_echo_fit = SE_fit_new( tes, np.array(sig), mean_noise, 'n' )
                T2bs.append(spin_echo_fit['T2'].value) 
            except RuntimeError: 
                print 'RuntimeError'
                T2bs.append(np.nan)
        T2_errs.append(np.std(T2bs))
    return T2s, T2_errs
    
def T1_ir_bootstrap(folders_to_process,N=1000,plot='n'):    
    """Given a folder of images will process IR data and return
    fitted T1 values and associated uncertainties.  Uncertainties 
    are obtained by bootstrapping pixels within each ROI"""  
    images=[read_dicoms(folder,['InversionTime']) for folder in folders_to_process]
    all_rois=[load_ROIs(folder+'/rois') for folder in folders_to_process]   
    TI_images=[img[0][0] for img in images]
    TIs=np.array([img[1][0]['InversionTime'] for img in images])
    T1s=[]
    T1_errs=[]
    for kk in np.arange(len(all_rois[0])-2):
        rois=[roi_list[kk] for roi_list in all_rois]
        sig=np.zeros(len(rois))
        T1bs=[]
        for mm, roi in enumerate(rois):
            sig[mm]=TI_images[mm][roi.get_indices()].mean()
        fit = IR_fit(TIs, sig) 
        T1s.append(fit['T1'].value)
        if plot=='y':
            plt.figure()
            plt.plot(TIs,sig,'o')
            TI_full=np.arange(0,8000,10)
            plt.plot(TI_full,fit(TI_full))               
        if N>0:  
            for nn in np.arange(N):
                for mm,roi in enumerate(rois):
                    npix=len(roi.get_indices()[0])
                    pixels=roi.get_indices()
                    ind=np.random.randint(npix,size=npix)
                    sig[mm]=TI_images[mm][pixels[0][ind],pixels[1][ind]].mean()                
                fit = IR_fit(TIs, sig) 
                T1bs.append(fit['T1'].value)
            #T1s.append(np.mean(T1bs))
            T1_errs.append(np.std(T1bs))
    return T1s, T1_errs
 
def check_T1_chi_surface(folders_to_process, sample):  
    """Check the chi-squared surface for a T1 fit"""
    images=[read_dicoms(folder,['InversionTime']) for folder in folders_to_process]
    all_rois=[load_ROIs(folder+'/rois') for folder in folders_to_process]   
    TI_images=[img[0][0] for img in images]
    TIs=np.array([img[1][0]['InversionTime'] for img in images])
    rois=[roi_list[sample] for roi_list in all_rois]
    sig=np.zeros(len(rois))
    for mm, roi in enumerate(rois):
        sig[mm]=TI_images[mm][roi.get_indices()].mean()
    fit = IR_fit(TIs, sig) 
    plt.figure()
    plt.plot(TIs,sig,'o')
    TI_full=np.arange(0,8000,10)
    plt.plot(TI_full,fit(TI_full))        
    #T1s=np.arange(fit['T1'].value-0.5*fit['T1'].value, fit['T1'].value+0.5*fit['T1'].value, fit['T1'].value/500)
    #aas=np.arange(fit['aa'].value-0.1*fit['aa'].value, fit['aa'].value+0.1*fit['aa'].value, fit['aa'].value/500)
    T1s=np.arange(800,2000,2)
    aas=np.arange(0.9,1.1,0.001)
    chi_surf=fit.contours(fit['T1'],fit['aa'],T1s,aas)
    return T1s, aas, chi_surf    
          
def make_case_report_csv(T2_dict,T1_dict,filename):
    book = xlwt.Workbook(encoding="utf-8")  
    sheet1 = book.add_sheet("sheet1")           
    key_list=T2_dict.keys()
    Hct_sO2s=np.array([[key[1],key[0],key[2]] for key in key_list])
    #Extract unique values
    Hct_sO2s_unique=np.unique(Hct_sO2s.view(np.dtype((np.void, Hct_sO2s.dtype.itemsize*Hct_sO2s.shape[1])))).view(Hct_sO2s.dtype).reshape(15,3)         
    for jj in np.arange(len(Hct_sO2s_unique)):
        T2_key_list=[]
        T2_value_list=[]
        for key in T2_dict:
            if (key[0]==Hct_sO2s_unique[jj][1] and key[1]==Hct_sO2s_unique[jj][0]):
                T2_key_list.append(key)
                T2_value_list.append(T2_dict[key])
        taus=np.array([key[3] for key in T2_key_list])
        indices=np.argsort(taus)
        
        T1_value_list=[]
        for key in T1_dict: 
            if (key[0]==Hct_sO2s_unique[jj][1] and key[1]==Hct_sO2s_unique[jj][0]):
                T1_value=(T1_dict[key]).nominal_value    
        if(jj==0):
            sheet1.write(jj,0,'Hct')
            sheet1.write(jj,1,'sO2')
            sheet1.write(jj,2,'met')
            sheet1.write(jj,3,'T1')
            for kk in np.arange(len(taus)):
                sheet1.write(jj,4+kk,taus[indices[kk]])
        
        sheet1.write(jj+1,0,Hct_sO2s_unique[jj][0])  
        sheet1.write(jj+1,1,100*Hct_sO2s_unique[jj][1])
        sheet1.write(jj+1,2,100*Hct_sO2s_unique[jj][2])
        sheet1.write(jj+1,3,T1_value)
        for kk in np.arange(len(taus)):
            sheet1.write(jj+1,4+kk,T2_value_list[indices[kk]].nominal_value)
    book.save(filename)     
    
def make_T2_report_csv(T2_dict_list,filename):
    book = xlwt.Workbook(encoding="utf-8")  
    sheets=['18apr2015','25apr2015','23may2015','6jun2015','27jun2015']
    for jj, T2_dict in enumerate(T2_dict_list):        
        sheet1 = book.add_sheet(sheets[jj])           
        key_list=T2_dict.keys()
        Hct_sO2s=np.array([[key[1],key[0],key[2]] for key in key_list])
        #Extract unique values
        Hct_sO2s_unique=np.unique(Hct_sO2s.view(np.dtype((np.void, Hct_sO2s.dtype.itemsize*Hct_sO2s.shape[1])))).view(Hct_sO2s.dtype).reshape(15,3)         
        for kk in np.arange(len(Hct_sO2s_unique)):
            T2_key_list=[]
            T2_value_list=[]
            for key in T2_dict:
                if (key[0]==Hct_sO2s_unique[kk][1] and key[1]==Hct_sO2s_unique[kk][0]):
                    T2_key_list.append(key)
                    T2_value_list.append(T2_dict[key])
            taus=np.array([key[3] for key in T2_key_list])
            indices=np.argsort(taus)
            
            if(kk==0):
                sheet1.write(kk,0,'Hct')
                sheet1.write(kk,1,'sO2')
                sheet1.write(kk,2,'met')

            sheet1.write(kk+1,0,Hct_sO2s_unique[kk][0])  
            sheet1.write(kk+1,1,100*Hct_sO2s_unique[kk][1])
            sheet1.write(kk+1,2,100*Hct_sO2s_unique[kk][2])
            
            for mm in np.arange(len(taus)):
                sheet1.write(kk+1,4+mm,T2_value_list[indices[mm]].nominal_value)
    book.save(filename)      
            
def calc_B0_map(folder,roi_mask):
    from human_SWI_data import phase_rescale
    import unwrap
    images,TE=read_dicoms(folder,'EchoTime')   
    phase_images = [phase_rescale(img) for img in images]
    roi_mask_inv=np.zeros_like(roi_mask)
    roi_mask_inv[roi_mask==1]=0
    roi_mask_inv[roi_mask==0]=1
    phase_images_unwrapped=[unwrap.unwrap(np.ma.masked_array(img,roi_mask_inv)) for img in phase_images]
    return phase_images_unwrapped
    
def rescale_and_unwrap_3D(img,roi_mask):
    #First rescale:
    img=(np.float32(img)-2048)*(np.pi/2048)       
    #Then unwrap:    
    import unwrap
    img_unwrapped=unwrap.unwrap(np.ma.masked_array(img,roi_mask))
    return img_unwrapped
    
def unwrap_3D(img, roi_mask):
    import unwrap
    img_unwrapped=unwrap.unwrap(np.ma.masked_array(img,roi_mask))
    return img_unwrapped
        
def calc_B1_map(folder1,folder2):
    """Given a pair of folders with 60, 120 degree excitation images, calculate a B1 map"""
    img1, FlipAngle1 = read_dicoms(folder1,['Flip Angle'])
    img2, FlipAngle2 = read_dicoms(folder2, ['FlipAngle'])
    img60=np.squeeze(np.array(img1).astype(float))
    img120=np.squeeze(np.array(img2).astype(float))
    
    FAmap=np.arccos(0.5*img120/img60)
    return img1,img2,FAmap
    
def T1_MOLLI_process(folder):    
    """Calculates T1s given a list of folders with MOLLI data"""
    data=img_roi_signal([folder],['InversionTime'])
    tis=np.array(data[2][0])
    mean_signal_mat=data[3]
    T1_list=[]
    rois=data[0][0]
    #plt.figure()
    for kk,roi in enumerate(rois[0:-2]):
        print kk
        #plt.figure()
        inversion_recovery = IR_fit(tis, mean_signal_mat[0,kk,:],[])
        plt.plot(tis, mean_signal_mat[0,kk,:],'o')
        TIs=np.arange(1,8000,5) 
        plt.plot(TIs,inversion_recovery(TIs))          
        T1_corr=inversion_recovery['T1'].value*(2*inversion_recovery['aa'].value-1)            
        T1_list.append(np.float64(T1_corr))
    return T1_list  
    
def T1_MOLLI_bootstrap(folder,N=1000):    
    """Given a folder of images will process IR data and return
    fitted T1 values and associated uncertainties.  Uncertainties 
    are obtained by bootstrapping pixels within each ROI"""  
    images,TI_dict=read_dicoms(folder,['InversionTime'])
    rois=load_ROIs(folder+'/rois')    
    TIs=np.array([entry['InversionTime'] for entry in TI_dict])
    sig=np.zeros(len(rois))
    T1s=[]
    T1_errs=[]
    for kk,roi in enumerate(rois[0:-2]):
        T1bs=[]
        npix=len(roi.get_indices()[0])
        print npix
        pixels=roi.get_indices()
        for nn in np.arange(N):
            ind=np.random.randint(npix,size=npix)
            sig=np.array([image[pixels[0][ind],pixels[1][ind]].mean() for image in images])
            fit = IR_fit(TIs, sig) 
            T1bs.append(fit['T1'].value*(2*fit['aa'].value-1))  
            #T1bs.append(fit['T1'].value)
        T1s.append(np.mean(T1bs))
        T1_errs.append(np.std(T1bs))
    return T1s, T1_errs   
    
def T1_MOLLI_bootstrap_everyTI(folder,N=1000):    
    """Given a folder of images will process IR data and return
    fitted T1 values and associated uncertainties.  Uncertainties 
    are obtained by bootstrapping pixels within each ROI"""  
    images,TI_dict=read_dicoms(folder,['InversionTime'])
    rois=load_ROIs(folder+'/rois')    
    TIs=np.array([entry['InversionTime'] for entry in TI_dict])
    sig=np.zeros(len(rois))
    T1s=[]
    T1_errs=[]
    for mm,roi in enumerate(rois[0:-2]):
        sig=[TI_image[roi.get_indices()].mean() for TI_image in images]               
        T1bs=[]
        npix=len(roi.get_indices()[0])
        print npix
        pixels=roi.get_indices()
        
        for nn in np.arange(N):
            sig_boot=np.zeros(len(TIs))
            for jj in np.arange(len(TIs)):        
                ind=np.random.randint(npix,size=npix)
                sig_boot[jj]=images[jj][pixels[0][ind],pixels[1][ind]].mean()
                fit = IR_fit(TIs, sig_boot) 
            T1bs.append(fit['T1'].value*(2*fit['aa'].value-1))  
            #T1bs.append(fit['T1'].value)
        T1s.append(np.mean(T1bs))
        T1_errs.append(np.std(T1bs))
    return T1s, T1_errs
    
def plot_T2_contour(esp, model):
    """Given a model (defined according to fitting.py) will plot constant sO2 contours as a function of Hct for 
    a given esp (in seconds)""" 
    Hct_full=np.arange(0.0,0.9,0.005)
    sO2_full=np.arange(0.0,1.1,0.005) 
    [Hct_mat,sO2_mat]=np.meshgrid(Hct_full,sO2_full)           
    ij=itertools.product(range(Hct_mat.shape[0]),range(Hct_mat.shape[1]))
    R2_ans=np.zeros_like(Hct_mat)
    xx=np.zeros(3)
    for element in ij:
        xx=np.array([Hct_mat[element],sO2_mat[element],esp,0.0])
        R2_ans[element[0],element[1]]=model(xx)    
    #plt.figure() 
    #C = plt.contour(Hct_mat,1000/R2_ans,100*sO2_mat,np.arange(0,110,10),vmin=0,vmax=110,linewidths=1.5,linestyles='solid',cmap=matplotlib.cm.jet)  
    C = plt.contour(Hct_mat,1000/R2_ans,100*sO2_mat,np.r_[0:40:20,30:110:10],vmin=0,vmax=110,linewidths=1.5,linestyles='solid',cmap=matplotlib.cm.cool)      
    plt.clabel(C, C.levels[::2], colors='k', inline=True, fmt="%0.1f", fontsize=14)                 

def plot_T1_contour(model,metHb=0.005):
    """Given a model (defined according to fitting.py) will plot constant sO2 contours as a function of Hct for 
    a given esp (in seconds)""" 
    Hct_full=np.arange(0.0,1.01,0.01)
    sO2_full=np.arange(-0.1,1.1,0.01) 
    [Hct_mat,sO2_mat]=np.meshgrid(Hct_full,sO2_full)           
    ij=itertools.product(range(Hct_mat.shape[0]),range(Hct_mat.shape[1]))
    R1_ans=np.zeros_like(Hct_mat)
    xx=np.zeros(3)
    for element in ij:
        xx=[Hct_mat[element],sO2_mat[element],metHb]
        R1_ans[element[0],element[1]]=model(xx)    
    #plt.figure() 
    #C = plt.contour(100*sO2_mat,1000/R1_ans,Hct_mat,np.arange(0.0,1.1,0.05),vmin=0,vmax=0.95,linewidths=1.5,cmap=matplotlib.cm.jet)  
    C = plt.contour(100*sO2_mat,1000/R1_ans,Hct_mat,np.r_[0:0.3:0.05,0.3:1.1:0.1],vmin=0,vmax=0.95,linewidths=1.5,cmap=matplotlib.cm.jet)  
    #C = plt.contour(Hct_mat,1000/R1_ans,sO2_mat,np.arange(-0.5,1.05,0.05),vmin=0,vmax=1.1,linewidths=1.5,cmap=matplotlib.cm.cool)     
    plt.clabel(C, C.levels[::2], colors='k', inline=True, fmt="%0.2f", fontsize=14)     

def retrieve_dict_data(T2_dict, esp=0, metHb_threshold=0.02):
    R2s=[]
    R2errs=[]
    sO2s=[]
    Hcts=[] 
    esps=[]
    metHbs=[]
    pt_nos=[]  #keep track of patient number
    T2s=[]
    T2errs=[]
    
    if (esp==0):  #retrieve all the data    
        dict_data=T2_dict.items()
        Hcts.append([item[0][1] for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        sO2s.append([item[0][0] for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        esps.append([0.001*item[0][3] for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        R2s.append([1000/item[1].nominal_value for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold]) 
        R2errs.append([(1000/item[1]).std_dev for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        metHbs.append([item[0][2] for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        pt_nos.append([item[0][4] for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        T2s.append([item[1].nominal_value for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
        T2errs.append([item[1].std_dev for item in dict_data if item[1].nominal_value!=-1.0 and item[0][2]<metHb_threshold])
    else: 
        for key in T2_dict:
            if(key[3]==esp and T2_dict[key].nominal_value!=-1 and key[2]<metHb_threshold):
                Hcts.append(key[1])
                esps.append(0.001*key[3])
                T2s.append(T2_dict[key].nominal_value)
                T2errs.append(T2_dict[key].std_dev)
                R2=1000/T2_dict[key]
                R2s.append(R2.nominal_value)
                R2errs.append(R2.std_dev)
                sO2s.append(key[0])
                metHbs.append(key[2])
                pt_nos.append(key[4])
                
    return sO2s, Hcts, esps, R2s, R2errs, metHbs, pt_nos, T2s, T2errs   

def retrieve_T1_dict_data(T1_dict, metHb_threshold=0.02):
    R1s=[]
    R1errs=[]
    sO2s=[]
    Hcts=[] 
    metHbs=[]
    pt_nos=[]
    T1s=[]
    T1errs=[]

    dict_data=T1_dict.items()
    Hcts.append([item[0][1] for item in dict_data if item[0][2]<metHb_threshold])
    sO2s.append([item[0][0] for item in dict_data if item[0][2]<metHb_threshold])
    metHbs.append([item[0][2] for item in dict_data if item[0][2]<metHb_threshold])
    pt_nos.append([item[0][3] for item in dict_data if item[0][2]<metHb_threshold])
    R1s.append([1000/item[1].nominal_value for item in dict_data if item[0][2]<metHb_threshold]) 
    R1errs.append([(1000/item[1]).std_dev for item in dict_data if item[0][2]<metHb_threshold])
    T1s.append([item[1].nominal_value for item in dict_data if item[0][2]<metHb_threshold])
    T1errs.append([item[1].std_dev for item in dict_data if item[0][2]<metHb_threshold])
  
    return sO2s, Hcts, metHbs, R1s, R1errs, pt_nos, T1s, T1errs

def retrieve_chi_dict_data(chi_dict, metHb_threshold=0.02):
    chis=[]
    sO2s=[]
    Hcts=[] 
    metHbs=[]
    pt_nos=[]
   
    dict_data=chi_dict.items()
    Hcts.append([item[0][1] for item in dict_data if item[0][2]<metHb_threshold])
    sO2s.append([item[0][0] for item in dict_data if item[0][2]<metHb_threshold])
    metHbs.append([item[0][2] for item in dict_data if item[0][2]<metHb_threshold])
    pt_nos.append([item[0][3] for item in dict_data if item[0][2]<metHb_threshold])
    chis.append([item[1] for item in dict_data if item[0][2]<metHb_threshold])
    return sO2s, Hcts, metHbs, pt_nos, chis

def calc_R2(x,a1,a2,a3,b1,b2,c1):
    R2=(a1+a2*x[0]+a3*x[0]**2) + (b1*x[0]+b2*x[0]**2)*(1-x[1]) + (c1*x[0]*(1-x[0]))*(1-x[1])**2      
    return R2
    
def get_T2_prep_times_VB17(dicom_list):
    prep_times_re = r'sWiPMemBlock.alFree\[2\][ ]+=[ ]([.\d]+)'
    prep_times_value_re = r'sWiPMemBlock.adFree\[\d+\][ ]+=[ ]([.\d]+)'
    
    for dicom_hdr in dicom_list:    
        out=(dicom_hdr[(0x0029,0x1020)])    
        prep_time_search = re.findall(prep_times_re,str(out.value))     
        
        if not prep_time_search:
            return []   
        else:
            num_preps=int(prep_time_search[0]) 
            value_list=re.findall(prep_times_value_re, str(out.value))[0:num_preps]
            prep_times=[float(prep_time)-0.001 for prep_time in value_list]
            return prep_times
    
def get_T2_prep_times_VE11(dicom_list):
    prep_times_re = r'sPrepPulses.adT2PrepDuration.__attribute__.size[ ]+=[ ]+([\d]+)'  
    prep_times_value_re = r'sPrepPulses.adT2PrepDuration\[\d+\][ ]+=[ ]+([.\d]+)'
    
    for dicom_hdr in dicom_list:
        out=repr(dicom_hdr[(0x0029,0x1020)].value)
        out=out.replace('\\t',' ')
        out=out.replace('\\n','\n')    
        prep_time_search = re.findall(prep_times_re,out)
        if not prep_time_search:
            return []
        else:
            num_preps=int(prep_time_search[0])
            value_list=re.findall(prep_times_value_re, out)[0:num_preps]
            prep_times=[float(prep_time)-0.001 for prep_time in value_list]
            return prep_times
    
def plateau_detect(sig, slope_threshold=-5, factor=2):
    """To eliminate the noise floor from a subsequent T2 fit, returns the index
    of the first element more than double the plateau value.  A plateau is defined
    as a region with slope>slope_threshold"""
    
    sig_diff=np.diff(sig)  
    plateau=np.where(sig_diff>-5)[0]
    
    if plateau.size:
       plateau_start=plateau[0]-1
       plateau_value=sig[plateau_start]
       threshold_value=plateau_value*2
    else:
       threshold_value=0
       
    temp=np.where(sig<threshold_value)[0]
    if temp.size:
        stop_index=temp[0]
    else:
        stop_index=len(sig) 
        
    return stop_index    
    
def bootstrap_fit(folder, roi,  model, x, iterations=100):
    image_list,TEs=read_dicoms(folder,['EchoTime'])    
    
    npix=len(roi.get_indices()[0])
    pixels=roi.get_indices()
    #sig=[img[pixels].mean() for img in image_list]
    #model.fit(x,sig)
    best_fit_pars=[p.value for p in model.pars]  #keep original best fit pars   
    print best_fit_pars
    
    
    par_records=[]
    for nn in np.arange(iterations):
        sig=[img[pixels].mean() for img in image_list]
        for qq,par in enumerate(model.pars):
            model[par.name].value=best_fit_pars[qq]  #reset to best fit pars
        ind=np.random.randint(npix,size=npix)
        sig=[img[pixels[0][ind],pixels[1][ind]].mean() for img in image_list]
        try:
            model.fit(x,sig)
        except RuntimeError:
            continue
        par_records.append([p.value for p in model.pars])
        
       # plt.plot(TE_full,model(TE_full),'r',alpha=0.5)
    
    par_uncertainties=np.zeros(len(best_fit_pars))
    for jj in np.arange(len(best_fit_pars)):
         par_hist=[record[jj] for record in par_records]  
         par_uncertainties[jj]=np.std(par_hist)
         
     
    return best_fit_pars, par_uncertainties 
    
def protocol_file_dump(foldername,filename):
    file_list=os.listdir(foldername)
    dicom_filename = foldername +'/' + file_list[0]
    hdr=dicom.read_file(dicom_filename)
    out=repr(hdr[(0x0029,0x1020)].value)
    out=out.replace('\\t','\t')
    out=out.replace('\\n','\n')
    f=open(filename,'w')
    f.write(out)
    f.close()
    
    return out
    
                

   

            