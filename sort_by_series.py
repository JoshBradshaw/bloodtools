# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:04:13 2015

@author: sportnoy
"""

import ROI_tools
import os
import shutil
import dicom
import sys

#parent_folder='/home/sportnoy/MRI_data/blood/25apr2015_trio/'

def folder_sort(parent_folder):
    
    scan_list=os.listdir(parent_folder)
    scan_list.sort()
    dicom_list=[scan for scan in scan_list if '.IMA' in scan]


    image_list,attributes=ROI_tools.read_dicoms(parent_folder,['SeriesDescription','EchoTime'])

    #series_descriptions=[attribute['SeriesDescription'] for attribute in attributes]
                
    for ii,scan in enumerate(dicom_list):
        series_desc=dicom.read_file(parent_folder+scan).SeriesDescription
        series_no=dicom.read_file(parent_folder+scan).SeriesNumber        
        foldername=parent_folder + series_desc + '_' + str(series_no)
        src=parent_folder + scan
        dst=foldername + '/' + scan
        try:
            os.mkdir(foldername)
            shutil.copy(src,dst)
        except OSError:
            shutil.copy(src,dst)

if __name__ == '__main__':
    foldername = (sys.argv[-1])
    folder_sort(foldername + '/')    