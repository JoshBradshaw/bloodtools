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
import os,sys
import pylab as plt

def MyPyQtSlot(*args):
    """
    for strange, mostly undocumented reasons, the default behaviour of PyQt
    is to silence all exceptions, which makes the debugging process roughly 
    equivilant in difficulty to trying to speak portugeuse when you don't know 
    any portugeuse. Inserting this decorator on every function is a kludgy way
    to make PyQt catch exceptions as you would except.
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
    @MyPyQtSlot("bool")    
    def __init__(self, parent=None):
        super(ROISelectPlot, self).__init__(parent)
        self.r = None
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        
        self.canvas = FigureCanvas(self.figure)
        #self.toolbar = NavigationToolbar2QT(self.canvas, self)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        #layout.addWidget(self.toolbar)
        self.setLayout(layout)
        
    @MyPyQtSlot("bool")
    def make_image(self, im):
        self.clear_roi()
        self.axes.clear()

        self.im = self.axes.imshow(im, vmin=np.percentile(im, 5),vmax=np.percentile(im, 95), cmap='gray')
        self.figure.canvas.draw()
        
    @MyPyQtSlot("bool")
    def new_roi(self):
        self.r = ROI.new_ROI(self.im, self.axes, self.figure)
        
    @MyPyQtSlot("bool")
    def clear_roi(self):
        print "clear ROI called!!"        
        
        if self.r is not None: # check if there is an ROI
            self.r.patch.remove()
            self.r.disconnect()
            self.axes.lines = []
            self.r = None
            self.figure.canvas.draw()
            
    @MyPyQtSlot("bool")
    def get_roi(self):
        return self.r
        
class T2CurvePlot(QtGui.QWidget):
    """
    The plot that displays the datapoints, the fitted monoexponential T2 curve, 
    and gives the value for T2 based on this fit.
    """
    @MyPyQtSlot("bool")    
    def __init__(self, parent=None):
        super(T2CurvePlot, self).__init__(parent)
        self.r = None
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        
        self.canvas = FigureCanvas(self.figure)
        #self.toolbar = NavigationToolbar2QT(self.canvas, self)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        #layout.addWidget(self.toolbar)
        self.setLayout(layout)
        
class MainWindow(QtGui.QWidget):
    @MyPyQtSlot("bool")
    def __init__(self):
        QtGui.QWidget.__init__(self, parent=None)
        
        button_load = QtGui.QPushButton(qtawesome.icon('fa.folder-open-o'), '')
        self.button_roi = QtGui.QCheckBox("ROI")
        button_run = QtGui.QPushButton(qtawesome.icon('fa.play'), '')
        self.button_roi.setChecked(False)
        self.combo_relax = QtGui.QComboBox()
        self.combo_relax.addItems(['T1', 'T2'])
        self.combo_relax.setCurrentIndex(0)
        
        self.plot_im = ROISelectPlot(self)
        self.plot_graph = T2CurvePlot(self)
        
        layout_top = QtGui.QHBoxLayout()
        layout_top.addWidget(button_load)
        layout_top.addWidget(self.button_roi)
        layout_top.addWidget(self.combo_relax)
        layout_top.addWidget(button_run)
        layout_top.addStretch()
        
        layout_mid = QtGui.QHBoxLayout()
        layout_mid.addWidget(self.plot_im)
        layout_mid.addWidget(self.plot_graph)
        
        layout_main = QtGui.QVBoxLayout()
        layout_main.addLayout(layout_top)
        layout_main.addLayout(layout_mid)
        self.setLayout(layout_main)
        
        button_load.clicked.connect(self.choose_dir)
        self.button_roi.stateChanged.connect(self.start_roi)
        button_run.clicked.connect(self.process_data)
    
    @MyPyQtSlot("bool")
    def choose_dir(self, event):
        out = QtGui.QFileDialog.getExistingDirectory(caption='MRI Data Directory')
        
        if out:
            self.directory = out
            self.images, self.image_attributes = blood_tools.read_dicoms(out,[])
            self.plot_im.make_image(self.images[0])
        else: # user hit the cancel or x button to leave the dialog
            pass
    
    @MyPyQtSlot("bool")
    def start_roi(self, state):
        if state:
            self.plot_im.new_roi()
        else:
            self.plot_im.clear_roi()
    
    @MyPyQtSlot("bool")
    def process_data(self, event):
        roi = self.plot_im.r 
        
        relax = self.combo_relax.currentText()
        # todo - implement this method
        x, y = make_T2_fit_from_directory_and_mask(self.directory, self.images, roi)
        self.plot_graph.axes.clear()
        self.plot_graph.axes.plot(x, y, 'ro')
        self.plot_graph.figure.canvas.draw()
        #x, y, parameters = make_fit_from_data(x, y, relax)
        #self.plot_graph.axes.plot(x, y, 'r-')
        #self.plot_graph.axes.figure.text(figx, figy, parameters)
        
                    
def make_T2_fit_from_directory_and_mask(dicom_directory, images, roi):
    """
    Simplest case: makes a T2 monoexponential fit given a directory of T2 
    DICOMs and a single ROI mask.
    """
    # get T2 prep times
    # VE11 is the new Siemens software, VE17 is the old version
    VE17_prep_times = blood_tools.get_T2_prep_times_VB17(dicom_directory)     
    VE11_prep_times = blood_tools.get_T2_prep_times_VE11(dicom_directory)

    print "VE11 preps: {}".format(VE11_prep_times)
    print "VE17 preps: {}".format(VE17_prep_times)    
    
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