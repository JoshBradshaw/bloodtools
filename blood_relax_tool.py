# -*- coding: utf-8 -*-
from __future__ import division
"""
Created on Tue May 24 13:09:08 2016

@author: Josh
"""

# std lib imports
import os
import traceback
import types
from functools import wraps
import cPickle
import math
# anaconda module imports
from PyQt4 import QtGui, QtCore
import qtawesome
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
# other python files that should be in the same dir
import ROI
import blood_tools
import fitting

def QTSlotExceptionRationalizer(*args):
    """
    for strange, mostly undocumented reasons, the default behaviour of PyQt
    is to silence all exceptions, which makes the debugging process roughly 
    equivilant in difficulty to trying to speak portugeuse when you don't know 
    any portugeuse. Inserting this decorator on every function is a kludgy way
    to make PyQt catch exceptions as you would except.
    
    DO NOT WRAP ANY METHODS OTHER THAN SLOT HANDLERS WITH THIS DECORATOR, IT
    BREAKS FUNCTION RETURNS!!!
    """
    if len(args) == 0 or isinstance(args[0], types.FunctionType):
        args = []
    @QtCore.pyqtSlot(*args)
    def slotdecorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                func(*args)
            except:
                print "Uncaught Exception in slot"
                traceback.print_exc()
        return wrapper

    return slotdecorator

class ROISelectPlot(QtGui.QWidget):
    """
    The plot canvas for the image, in the left pane of the GUI. This is where
    the plot is displayed, and the ROI is selected.
    
    To implement:
     - single ROI selection for multiple slices
     - multiple ROI selection for multiple slices
     - next/prev slice advancement with a progress indicator
    """
    @QTSlotExceptionRationalizer("bool")    
    def __init__(self, parent=None):
        super(ROISelectPlot, self).__init__(parent)
        # initialize the plot area        
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=False)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)
        self.im = None
        
    @QTSlotExceptionRationalizer("bool")
    def make_image(self, im, vmin=5, vmax=95):
        # only turn autoscale on when setting the image so that ROI changes won't tweak the autoscale
        self.im = im
        self.axes.set_autoscale_on(True)
        self.mpl_im = self.axes.imshow(im, vmin=np.percentile(im, vmin),vmax=np.percentile(im, vmax), cmap='gray', origin='upper')
        self.axes.set_autoscale_on(False)
        self.figure.canvas.draw()
        
    def get_axes(self):
        return self.axes
        
    def get_mpl_im(self):
        return self.mpl_im
        
    def get_figure(self):
        return self.figure
        
class ColourROISelectPlot(ROISelectPlot):
    """
    Identical to ROI select plot, except colourized
    """
    @QTSlotExceptionRationalizer("bool")
    def make_image(self, im, vmin=5, vmax=95):      
        # only turn autoscale on when setting the image so that ROI changes won't tweak the autoscale
        self.im = im
        self.axes.set_autoscale_on(True)
        self.mpl_im = self.axes.imshow(im, vmin=np.percentile(im, vmin),vmax=np.percentile(im, vmax), cmap='jet', origin='upper')
        self.axes.set_autoscale_on(False)
        self.figure.canvas.draw()

class T2CurvePlot(QtGui.QWidget):
    """
    The plot that displays the datapoints, the fitted monoexponential T2 curve, 
    and gives the value for T2 based on this fit.
    """
    @QTSlotExceptionRationalizer("bool")    
    def __init__(self, parent=None):
        super(T2CurvePlot, self).__init__(parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)
        
class MainWindow(QtGui.QWidget):
    @QTSlotExceptionRationalizer("bool")
    def __init__(self):
        # draw the interface
        self.init_gui()
        # initialize these to empty lists so that the slice select display will work before loading images
        self.images, self.image_attributes = [], []    
        # used for keeping track of what slice is displayed
        self.image_index = 0
        self.image_filename = None
        self.image_ROIs = {}
        self.image_filename_list = []
        self.grey_activeROI = None
        self.color_activeROI = None
        self.roi_path = None
        self.directory = ""
        self.vmin = 5
        self.vmax = 95
        self.grey_roi_patch = None
        self.color_roi_patch = None
        self.controls_enabled(False)
        self.included_slices = []

    @QTSlotExceptionRationalizer("bool")
    def init_gui(self):
        QtGui.QWidget.__init__(self, parent=None)
        # button for opening dicom directory
        self.button_load = QtGui.QPushButton(qtawesome.icon('fa.folder-open-o'), '')
        self.button_run = QtGui.QPushButton('Fit Data') 
        self.button_draw_roi = QtGui.QPushButton('Draw ROI') 
        self.button_exclude_slice = QtGui.QPushButton('Exclude Slice') 
        # first/last prev/next buttons for scrolling through the dicoms
        self.button_image_fwd = QtGui.QPushButton(qtawesome.icon('fa.chevron-right'), '')
        self.button_image_fwd.setObjectName('slice_fwd')
        self.button_image_bwd = QtGui.QPushButton(qtawesome.icon('fa.chevron-left'), '')
        self.button_image_bwd.setObjectName('slice_bwd')
        
        self.button_image_first = QtGui.QPushButton(qtawesome.icon('fa.step-backward'), '')
        self.button_image_first.setObjectName('slice_first')
        self.button_image_last = QtGui.QPushButton(qtawesome.icon('fa.step-forward'), '')
        self.button_image_last.setObjectName('slice_last')
        self.slice_label = QtGui.QLabel('(00/00)')
        # set whether the ROI should apply to a single slice or all of the slices
        self.combo_roi_scope = QtGui.QComboBox()   
        self.combo_roi_scope.addItems(['This Slice','All Slices'])
        self.combo_roi_scope.setCurrentIndex(0)
        # choose ROI shape
        self.combo_roi_style = QtGui.QComboBox()
        self.combo_roi_style.addItems(['Polygon','Circle', 'Ellipse'])
        self.combo_roi_style.setCurrentIndex(0)
        # fit either a basic T1, or basic T2 fit
        self.combo_relax_label = QtGui.QLabel('Fit Type')
        self.combo_relax = QtGui.QComboBox()
        self.combo_relax.addItems(['T1', 'T2'])
        self.combo_relax.setCurrentIndex(0)
        
        self.plot_im = ROISelectPlot(self)
        self.color_plot_im = ColourROISelectPlot(self)
        self.plot_graph = T2CurvePlot(self)
        
        layout_top = QtGui.QHBoxLayout()
        layout_top.addSpacing(10)
        layout_top.addWidget(self.button_load)
        
        layout_top.addStretch()
        layout_top.addWidget(QtGui.QLabel('Change Slice:'))

        layout_top.addWidget(self.button_image_first)
        layout_top.addWidget(self.button_image_bwd)
        layout_top.addWidget(self.button_image_fwd) 
        layout_top.addWidget(self.button_image_last)
        layout_top.addWidget(self.slice_label)
        
        layout_top.addStretch()
        layout_top.addWidget(QtGui.QLabel('ROI Style:'))
        layout_top.addWidget(self.combo_roi_style)
        layout_top.addWidget(QtGui.QLabel('Apply ROI to:'))
        layout_top.addWidget(self.combo_roi_scope)
        layout_top.addWidget(self.button_draw_roi)
        layout_top.addWidget(self.button_exclude_slice)
        layout_top.addStretch()
        layout_top.addWidget(self.combo_relax_label)
        layout_top.addWidget(self.combo_relax)
        layout_top.addWidget(self.button_run)
        layout_top.addSpacing(10)
        
        layout_mid = QtGui.QHBoxLayout()
        layout_mid.addWidget(self.plot_im)
        layout_mid.addWidget(self.color_plot_im)
        layout_mid.addWidget(self.plot_graph)
        
        self.vmin_window_slider = QtGui.QSlider(orientation=QtCore.Qt.Horizontal, minimum=0, maximum=100)
        self.vmin_window_slider.setValue(5)
        self.vmax_window_slider = QtGui.QSlider(orientation=QtCore.Qt.Horizontal, minimum=0, maximum=100)
        self.vmax_window_slider.setValue(95)
        layout_slider1 = QtGui.QHBoxLayout()
        layout_slider1.addWidget(QtGui.QLabel('Window Min:'))
        layout_slider1.addSpacing(3)
        layout_slider1.addWidget(self.vmin_window_slider)
        
        layout_slider2 = QtGui.QHBoxLayout()
        layout_slider2.addWidget(QtGui.QLabel('Window Max:'))
        layout_slider2.addWidget(self.vmax_window_slider)        
        
        layout_main = QtGui.QVBoxLayout()
        layout_main.addLayout(layout_top)
        layout_main.addLayout(layout_mid)
        layout_main.addLayout(layout_slider1)
        layout_main.addLayout(layout_slider2)
        self.setLayout(layout_main)
        
        self.button_load.pressed.connect(self.choose_dir)
        self.button_run.pressed.connect(self.process_data)
        self.button_draw_roi.pressed.connect(self.start_roi)
        self.button_exclude_slice.pressed.connect(self.exclude_current_slice)
        self.button_image_first.pressed.connect(self.change_image)
        self.button_image_last.pressed.connect(self.change_image)
        self.button_image_fwd.pressed.connect(self.change_image)
        self.button_image_bwd.pressed.connect(self.change_image)
        self.vmin_window_slider.valueChanged.connect(self.set_image_window)
        self.vmax_window_slider.valueChanged.connect(self.set_image_window)

    @QTSlotExceptionRationalizer("bool")
    def exclude_current_slice(self, *e):
        if self.included_slices[self.image_index]:
            self.included_slices[self.image_index] = False
            self.clear_roi()
            self.roi_controls_enable(False)
        else:
            self.included_slices[self.image_index] = True
            self.roi_controls_enable(True)
            self.load_roi()
            
    @QTSlotExceptionRationalizer("bool")
    def roi_controls_enable(self, enable=True):
        if enable:
            self.button_exclude_slice.setText('Exclude Slice')
        else:
            self.button_exclude_slice.setText('Include Slice')
            
        self.combo_roi_style.setEnabled(enable)
        self.combo_relax.setEnabled(enable)
        self.combo_roi_scope.setEnabled(enable)
        self.button_draw_roi.setEnabled(enable)

    @QTSlotExceptionRationalizer("bool")
    def controls_enabled(self, enable=True):
        self.button_run.setEnabled(enable)
        self.button_draw_roi.setEnabled(enable)
        self.button_image_first.setEnabled(enable)
        self.button_image_last.setEnabled(enable)
        self.button_image_fwd.setEnabled(enable)
        self.button_image_bwd.setEnabled(enable)
        self.combo_roi_style.setEnabled(enable)
        self.combo_relax.setEnabled(enable)
        self.combo_roi_scope.setEnabled(enable)
    
    @QTSlotExceptionRationalizer("bool")
    def set_image_window(self, *e):
        self.vmin = self.vmin_window_slider.value()
        self.vmax = self.vmax_window_slider.value()
        
        if self.vmin <= self.vmax and self.plot_im.im is not None:
            im_vmin = np.percentile(self.plot_im.im, self.vmin)
            im_vmax = np.percentile(self.plot_im.im, self.vmax)

            self.plot_im.mpl_im.set_clim(im_vmin, im_vmax)
            self.color_plot_im.mpl_im.set_clim(im_vmin, im_vmax)
            
            self.plot_im.figure.canvas.draw()
            self.color_plot_im.figure.canvas.draw()
        else: # matplotlib will throw an error if the window is negative
            pass
    
    @QTSlotExceptionRationalizer("bool")
    def clear_roi(self):        
        # remove the ROI from the screen, but do not delete it until it is
        # overwritten by another ROI
        # remove the ROI patches created when loading if necessary
        if self.grey_roi_patch is not None:        
            self.grey_roi_patch.remove()
            self.grey_roi_patch = None
        if self.color_roi_patch is not None:
            self.color_roi_patch.remove()
            self.color_roi_patch = None
            
        if self.color_activeROI is not None:
            self.color_activeROI.remove()
            axes = self.color_plot_im.get_axes()
            axes = []
            self.color_activeROI = None
        if self.grey_activeROI is not None: # check if there is an ROI
            self.grey_activeROI.remove()
            axes = self.plot_im.get_axes()
            axes = []
            self.grey_activeROI = None

        grey_figure = self.plot_im.get_figure()
        grey_figure.canvas.draw()
        color_figure = self.color_plot_im.get_figure()
        color_figure.canvas.draw()
        
    @QTSlotExceptionRationalizer("bool")
    def change_image(self):
        # serialize existing ROIs to file, this is quick+dirty b/c I haven't
        # figured out how to detect the ROI complete event in this class
        # so I just save them when the user changes images
        sender_btn = self.sender().objectName()
        num_images = len(self.images)
        
        # prevent the buttons from raising div by zero exceptions when no images loaded
        if num_images > 0:
            if sender_btn == 'slice_first':
                self.image_index = 0
            elif sender_btn == 'slice_bwd':
                ind = self.image_index - 1
                self.image_index = ind % num_images
            elif sender_btn == 'slice_fwd':
                ind = self.image_index + 1
                self.image_index = ind % num_images
            else:
                self.image_index = num_images - 1
        
            self.image_filename = self.image_filename_list[self.image_index]
            # display the slice selection label, with zero padding to keep the toolbar from shifting around
            num, demon = str(self.image_index+1).rjust(2, '0'), str(num_images).rjust(2, '0')
            # display previous ROI if it exists
            self.clear_roi()
            self.plot_im.mpl_im.set_data(self.images[self.image_index])
            self.color_plot_im.mpl_im.set_data(self.images[self.image_index])
            self.set_image_window()

            self.plot_im.figure.canvas.draw()
            self.color_plot_im.figure.canvas.draw()
            self.load_roi()
            self.slice_label.setText("{}/{}".format(num, demon))
            
    def save_analysis(self):
        """save ROIs and any other essential settings to a .ROIs file"""
        to_save = {}
        to_save['ROIs'] = self.image_ROIs
        to_save['included_slices'] = self.included_slices        
        
        with open(self.roi_path, 'w') as f:
            cPickle.dump(to_save, f)
            
    def load_prev_analysis(self):
        quit_msg = "Would you like to reload your previous ROIs for this series?"
        reply = QtGui.QMessageBox.question(self, 'Message', 
                         quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
    
        if reply == QtGui.QMessageBox.Yes:
            # load image ROIs if possible
            if os.path.isfile(self.roi_path):
                with open(self.roi_path, 'r') as f:
                    to_load = cPickle.load(f)
                    self.image_ROIs = to_load['ROIs']
                    self.included_slices = to_load['included_slices']
    
    @QTSlotExceptionRationalizer("bool")
    def grey_roi_complete_callback(self):
        roi_scope = self.get_roi_scope()
        
        if roi_scope == "All Slices":
            for img_fn in self.image_filename_list:
                self.image_ROIs[img_fn] = self.grey_activeROI
        else:
            self.image_ROIs[self.image_filename] = self.grey_activeROI
            
        self.save_analysis()
        
        self.clear_roi()
        self.load_roi()

    @QTSlotExceptionRationalizer("bool")
    def color_roi_complete_callback(self): 
        roi_scope = self.get_roi_scope()       
        
        if roi_scope == "All Slices":
            for img_fn in self.image_filename_list:
                self.image_ROIs[img_fn] = self.color_activeROI
        else:
            self.image_ROIs[self.image_filename] = self.color_activeROI
        # save ROI files
        self.save_analysis()
        
        self.clear_roi()
        self.load_roi()
        
    def get_roi_scope(self):
        return self.combo_roi_scope.currentText()

    def get_roi_style(self):
        current_text = self.combo_roi_style.currentText()
        return current_text
    
    def get_relax_type(self):
        return self.combo_relax.currentText()
    
    @QTSlotExceptionRationalizer("bool")
    def choose_dir(self, *event):        
        """opens a directory choose dialog box, allows the user to select their
        dicom series of interest and loads that series."""  
        self.plot_graph.axes.clear()
        self.plot_graph.figure.canvas.draw()
        
        if self.directory:
            out = QtGui.QFileDialog.getExistingDirectory(directory=os.path.split(self.directory)[0], caption='MRI Data Directory')
        else:
            out = QtGui.QFileDialog.getExistingDirectory(caption='MRI Data Directory')
        
        if out:
            self.image_index = 0
            self.image_filename = None
            self.image_ROIs = {}
            self.image_filename_list = []
            self.grey_activeROI = None
            self.color_activeROI = None
            self.roi_path = None
            self.images, self.image_attributes = [], []
            self.plot_graph.axes.clear()
            self.clear_roi()
            self.directory = out
            self.roi_path = os.path.join(self.directory, '.ROIs')
            self.images, self.image_attributes, self.dicom_list = blood_tools.read_dicoms(out, ['InversionTime'])
            # initiate included slices to be all True
            self.included_slices = [True for _ in range(len(self.images))]            
            
            if not self.images:
                error = QtGui.QErrorMessage()
                error.showMessage('The selected directory does not contain a DICOM series which this widget is capable of loading')
                error.exec_()
                return
            
            for attributes in self.image_attributes:
                self.image_filename_list.append(attributes['filename'])

            self.image_filename = self.image_filename_list[self.image_index]
            self.plot_im.make_image(self.images[self.image_index], self.vmin, self.vmax)
            self.controls_enabled(True)
            self.color_plot_im.make_image(self.images[self.image_index], self.vmin, self.vmax)
            num, demon = '01', str(len(self.images)).rjust(2, '0')
            self.slice_label.setText("{}/{}".format(num, demon))
            
            self.load_prev_analysis()
            self.load_roi()
        else: # user hit the cancel or x button to leave the dialog
            pass
    
    @QTSlotExceptionRationalizer("bool")
    def load_roi(self):
        """if the active image has a previously drawn ROI, this method reloads its"""
        if not self.included_slices[self.image_index]:
            self.roi_controls_enable(False)
            return
        else:
            self.roi_controls_enable(True)
        
        if self.image_filename in self.image_ROIs:
            self.grey_activeROI = self.image_ROIs[self.image_filename]
            self.color_activeROI = self.image_ROIs[self.image_filename]
            self.grey_roi_patch = self.grey_activeROI.draw(self.plot_im.axes, self.plot_im.figure, 'red')
            self.color_roi_patch = self.color_activeROI.draw(self.color_plot_im.axes, self.color_plot_im.figure, 'black')
    
    @QTSlotExceptionRationalizer("bool")
    def start_roi(self):
        """create a new ROI for the image"""
        if not len(self.images) > 0:
            error = QtGui.QErrorMessage()
            error.showMessage('You must a load a series of images before drawing the ROI')
            error.exec_()
            return
           
        roi_style = self.get_roi_style().lower() # style names are lowercase in ROI.py
        
        self.clear_roi()
        # create an ROI object for both images, keep the one that calls the complete callback first
        self.grey_activeROI = ROI.new_ROI(self.plot_im.get_mpl_im(), 
            self.plot_im.get_axes(), self.plot_im.get_figure(), 
            roi_style, 'red', self.grey_roi_complete_callback)
            
        self.color_activeROI = ROI.new_ROI(self.color_plot_im.get_mpl_im(), 
            self.color_plot_im.get_axes(), self.color_plot_im.get_figure(), 
            roi_style, 'black', self.color_roi_complete_callback)
        
    @QTSlotExceptionRationalizer("bool")
    def process_data(self, *event):
        # check that user has drawn all of the required ROIs
        for ii, fn in enumerate(self.image_filename_list):
            if self.included_slices[ii] and not fn in self.image_ROIs: 
                error = QtGui.QErrorMessage()
                error.showMessage("You must draw an ROI on every included slice before you can fit the data. Use the 'All Slices' option if the ROIs are conincident accross the slices")
                error.exec_()
                return
                
        image_attributes = []
        images = []
        roi_list = []
        dicom_list = []
        
        for image_index, include_slice in enumerate(self.included_slices):
            image_filename = self.image_filename_list[image_index]
            if include_slice:
                image_attributes.append(self.image_attributes[image_index])
                images.append(self.images[image_index])
                roi_list.append(self.image_ROIs[image_filename])
                dicom_list.append(self.dicom_list[image_index])
                
                
        if not len(self.images) > 0:
            error = QtGui.QErrorMessage()
            error.showMessage('You must load a series of dicom images before fitting the data')
            error.exec_()
            return

        # todo add error message if images or ROI not loaded
        relaxation_type = self.get_relax_type()      
        axes = self.plot_graph.axes
        axes.clear()

        if relaxation_type == 'T1':
            ti, signal = get_T1_decay_signal(self.image_attributes, self.images, roi_list, self.included_slices, log_scale=False)
            if not len(ti) or ti[0] == 0:
                error = QtGui.QErrorMessage()
                error.showMessage('Failed to find T1 recovery times for this dataset, ensure that this is a T1 series')
                error.exec_()
                return
            
            axes.plot(ti, signal, 'ro')
            axes.set_xlabel('inversion time (ms)')
            axes.set_ylabel('signal intensity')
            
            inversion_recovery = fitting.model('abs(M0*(1-2*aa*exp(-x/T1)))', {'M0':signal.max(),'aa':1,'T1':1000}) 
            #Initial T1 guess based on null time    
            T1_guess=-ti[np.where(signal==signal.min())]/np.log(0.5) 
            if(np.size(T1_guess)>1):
                inversion_recovery['T1']=T1_guess[0]
            else:
                inversion_recovery['T1']=T1_guess  
                
            inversion_recovery.fit(ti, signal)
            fix_x_points=np.arange(0,20000,1)
            axes.plot(fix_x_points,inversion_recovery(fix_x_points), 'r-')          
            T1_corr=inversion_recovery['T1'].value*(2*inversion_recovery['aa'].value-1)
            axes.text(0.3, 0.9, "T1 Value: {}ms".format(round(T1_corr)), transform=axes.transAxes)
        elif relaxation_type == 'T2':    
            x, y = get_T2_decay_signal(self.dicom_list, self.images, roi_list, self.included_slices)

            if not len(x):
                error = QtGui.QErrorMessage()
                error.showMessage('Failed to find T2 recovery times for this dataset, ensure that this is a T2 series')
                error.exec_()
                return            
            
            axes.plot(x, y, 'ro')
            axes.set_xlabel('TE time (ms)')
            axes.set_ylabel('log(signal intensity)')
                        
            first_coeff, zeroth_coeff = np.polyfit(x, y, 1)
            start, stop = x[0], x[-1]
            fit_x_points = np.linspace(start, stop, 1000)
            fit_y_points = [first_coeff * x_ + zeroth_coeff for x_ in fit_x_points]
            axes.plot(fit_x_points, fit_y_points, 'b-')
            t2_value = -1/first_coeff
            axes.text(0.3, 0.9, "T2 Value: {}ms".format(round(t2_value)), transform=axes.transAxes)
        else:
            raise NotImplementedError("This mapping type's fitting algorithm has not been implemented yet") 
        axes.set_xlim(xmin=0)
        self.plot_graph.figure.canvas.draw()
            
def get_T2_decay_signal(dicom_list, image_list, roi_list, included_slices, log_scale=True):
    """Simplest case: makes a T2 monoexponential fit given a directory of T2 
    DICOMs and a single ROI mask.
    """
    # get T2 prep times
    # VE11 is the new Siemens software, VE17 is the old version
    VE17_prep_times = blood_tools.get_T2_prep_times_VB17(dicom_list)
    VE11_prep_times = blood_tools.get_T2_prep_times_VE11(dicom_list)  
    
    if VE11_prep_times:
        prep_times = VE11_prep_times
    else:
        prep_times = VE17_prep_times
        
    for ii, slice_included in enumerate(included_slices):
        if not slice_included:
            prep_times.pop(ii)        
    
    signal = []
    log_signal = []
    for image, roi in zip(image_list, roi_list):
        mean_signal_magnitude_over_ROI = blood_tools.calc_ROI_mean(roi, image)
        signal.append(mean_signal_magnitude_over_ROI)
        log_signal.append(math.log(mean_signal_magnitude_over_ROI))
    if log_scale:
        return prep_times, log_signal
    else:
        return prep_times, signal
        
def get_T1_decay_signal(image_attributes, image_list, roi_list, included_slices, log_scale=True):
    """Simplest case: makes a T2 monoexponential fit given a directory of T2 
    DICOMs and a single ROI mask.
    """
    inversion_times = np.array([att['InversionTime'] for att in image_attributes])  
    
    for ii, slice_included in enumerate(included_slices):
        if not slice_included:
            inversion_times.pop(ii)    
    
    signal = []
    log_signal = []
    for image, roi in zip(image_list, roi_list):
        mean_signal_magnitude_over_ROI = blood_tools.calc_ROI_mean(roi, image)
        signal.append(mean_signal_magnitude_over_ROI)
        log_signal.append(math.log(mean_signal_magnitude_over_ROI))
    
    if log_scale:
        return inversion_times, np.array(log_signal)
    else:
        return inversion_times, np.array(signal)

def main():
    app = QtGui.QApplication.instance() or QtGui.QApplication([])
    win = MainWindow()
    win.show()
    app.exec_()
    return win

if __name__ == '__main__':
    win = main()