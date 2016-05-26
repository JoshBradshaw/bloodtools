# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:09:08 2016

@author: Josh
"""

from PyQt4 import QtGui, QtCore
import qtawesome
import traceback
import types
from functools import wraps
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import ROI
import blood_tools
import numpy as np
from pprint import pprint

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
        self.image_ROI = None
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)
        
    @QTSlotExceptionRationalizer("bool")
    def make_image(self, im):
        self.clear_roi()
        self.axes.clear()
        
        # only turn autoscale on when setting the image so that ROI changes won't tweak the autoscale
        self.axes.set_autoscale_on(True)
        self.mpl_im = self.axes.imshow(im, vmin=np.percentile(im, 5),vmax=np.percentile(im, 95), cmap='gray')
        self.figure.canvas.draw()
        self.axes.set_autoscale_on(False)
        
    @QTSlotExceptionRationalizer("bool")
    def clear_roi(self):        
        if self.image_ROI is not None: # check if there is an ROI
            self.image_ROI.remove()
            self.image_ROI.disconnect()
            self.axes.lines = []
            self.figure.canvas.draw()
            self.image_ROI = None
            
    def get_roi(self):
        return self.roi
        
    def get_axes(self):
        return self.axes
        
    def get_mpl_im(self):
        return self.mpl_im
        
    def get_figure(self):
        return self.figure
        
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

    @QTSlotExceptionRationalizer("bool")
    def init_gui(self):
        QtGui.QWidget.__init__(self, parent=None)
        
        button_load = QtGui.QPushButton(qtawesome.icon('fa.folder-open-o'), '')
        button_run = QtGui.QPushButton('Fit Data') 
        button_draw_roi = QtGui.QPushButton('Draw ROI') 
        
        button_slice_fwd = QtGui.QPushButton(qtawesome.icon('fa.chevron-right'), '')
        button_slice_fwd.setObjectName('slice_fwd')
        button_slice_bwd = QtGui.QPushButton(qtawesome.icon('fa.chevron-left'), '')
        button_slice_bwd.setObjectName('slice_bwd')
        
        button_slice_first = QtGui.QPushButton(qtawesome.icon('fa.step-backward'), '')
        button_slice_first.setObjectName('slice_first')
        button_slice_last = QtGui.QPushButton(qtawesome.icon('fa.step-forward'), '')
        button_slice_last.setObjectName('slice_last')
        self.slice_label = QtGui.QLabel('(00/00)')
        
        self.combo_roi_scope = QtGui.QComboBox()   
        self.combo_roi_scope.addItems(['This Slice','All Slices'])
        self.combo_roi_scope.setCurrentIndex(0)

        self.combo_roi_style = QtGui.QComboBox()
        self.combo_roi_style.addItems(['Polygon','Circle', 'Ellipse'])
        self.combo_roi_style.setCurrentIndex(0)
        
        combo_relax_label = QtGui.QLabel('Fit Type')
        self.combo_relax = QtGui.QComboBox()
        self.combo_relax.addItems(['T1', 'T2'])
        self.combo_relax.setCurrentIndex(0)
        
        self.plot_im = ROISelectPlot(self)
        self.plot_graph = T2CurvePlot(self)
        
        layout_top = QtGui.QHBoxLayout()
        layout_top.addSpacing(10)
        layout_top.addWidget(button_load)
        
        layout_top.addStretch()
        layout_top.addWidget(QtGui.QLabel('Change Slice:'))

        layout_top.addWidget(button_slice_first)
        layout_top.addWidget(button_slice_bwd)
        layout_top.addWidget(button_slice_fwd) 
        layout_top.addWidget(button_slice_last)
        layout_top.addWidget(self.slice_label)
        
        layout_top.addStretch()
        layout_top.addWidget(QtGui.QLabel('ROI Style:'))
        layout_top.addWidget(self.combo_roi_style)
        layout_top.addWidget(QtGui.QLabel('Apply ROI to:'))
        layout_top.addWidget(self.combo_roi_scope)
        layout_top.addWidget(button_draw_roi)
        layout_top.addStretch()
        layout_top.addWidget(combo_relax_label)
        layout_top.addWidget(self.combo_relax)
        layout_top.addWidget(button_run)
        
        layout_mid = QtGui.QHBoxLayout()
        layout_mid.addWidget(self.plot_im)
        layout_mid.addWidget(self.plot_graph)
        layout_top.addSpacing(10)
        
        layout_main = QtGui.QVBoxLayout()
        layout_main.addLayout(layout_top)
        layout_main.addLayout(layout_mid)
        self.setLayout(layout_main)
        
        button_load.clicked.connect(self.choose_dir)
        button_run.clicked.connect(self.process_data)
        button_draw_roi.pressed.connect(self.start_roi)
        button_slice_first.pressed.connect(self.change_image)
        button_slice_last.pressed.connect(self.change_image)
        button_slice_fwd.pressed.connect(self.change_image)
        button_slice_bwd.pressed.connect(self.change_image)

    @QTSlotExceptionRationalizer("bool")
    def change_image(self):
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
        
        # display the slice selection label, with zero padding to keep the toolbar from shifting around
        num, demon = str(self.image_index+1).rjust(2, '0'), str(num_images).rjust(2, '0')
        self.slice_label.setText("{}/{}".format(num, demon))
        self.plot_im.make_image(self.images[self.image_index])
        
    def get_roi_scope(self):
        return self.combo_roi_scope.currentText()

    def get_roi_style(self):
        current_text = self.combo_roi_style.currentText()
        return current_text
    
    def get_relax_type(self):
        return self.combo_relax.currentText()
    
    @QTSlotExceptionRationalizer("bool")
    def choose_dir(self, event):
        out = QtGui.QFileDialog.getExistingDirectory(caption='MRI Data Directory')
        if out:
            self.directory = out
            self.images, self.image_attributes = blood_tools.read_dicoms(out,[])
            self.plot_im.make_image(self.images[0])
            
            num, demon = '01', str(len(self.images)).rjust(2, '0')
            self.slice_label.setText("{}/{}".format(num, demon))
        else: # user hit the cancel or x button to leave the dialog
            pass
    
    @QTSlotExceptionRationalizer("bool")
    def start_roi(self):
        if not len(self.images) > 0:
            error = QtGui.QErrorMessage()
            error.showMessage('You must a load a series of images before drawing the ROI')
            error.exec_()        
        
        roi_style = self.get_roi_style().lower() # style names are lowercase in ROI.py
        roi_scope = self.get_roi_scope()
        
        self.plot_im.clear_roi()
        self.plot_im.image_ROI = ROI.new_ROI(self.plot_im.get_mpl_im(), 
        self.plot_im.get_axes(), self.plot_im.get_figure(), shape=roi_style)
        
    @QTSlotExceptionRationalizer("bool")
    def process_data(self, event):
        # todo add error message if images or ROI not loaded
        relaxation_type = self.get_relax_type()         
        roi = self.plot_im.get_roi()
        relax = self.combo_relax.currentText()
        # todo - implement this method
        
        if relaxation_type == 'T1':
            raise NotImplementedError("T1 mapping has not been implemented yet")            
        elif relaxation_type == 'T2':    
            x, y = make_T2_fit_from_directory_and_roi(self.directory, self.images, roi)
            self.plot_graph.axes.clear()
            self.plot_graph.axes.plot(x, y, 'ro')
            self.plot_graph.figure.canvas.draw()
            #x, y, parameters = make_fit_from_data(x, y, relax)
            #self.plot_graph.axes.plot(x, y, 'r-')
            #self.plot_graph.axes.figure.text(figx, figy, parameters)
        else:
            raise NotImplementedError("This mapping type's fitting algorithm has not been implemented yet") 
            
def make_T2_fit_from_directory_and_roi(dicom_directory, images, roi):
    """Simplest case: makes a T2 monoexponential fit given a directory of T2 
    DICOMs and a single ROI mask.
    """
    # get T2 prep times
    # VE11 is the new Siemens software, VE17 is the old version
    VE17_prep_times = blood_tools.get_T2_prep_times_VB17(dicom_directory)     
    VE11_prep_times = blood_tools.get_T2_prep_times_VE11(dicom_directory)  
    
    if VE11_prep_times:
        prep_times = VE11_prep_times
    else:
        prep_times = VE17_prep_times
    
    mean_signal_magnitude = blood_tools.calc_ROI_mean(roi, images)
    
    print "Prep times: {}".format(prep_times)
    print "T2 signal intensities: {}".format(mean_signal_magnitude)    
    return prep_times, mean_signal_magnitude


def main():
    app = QtGui.QApplication.instance() or QtGui.QApplication([])
    win = MainWindow()
    win.show()
    app.exec_()
    return win

if __name__ == '__main__':
    win = main()