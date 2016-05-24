
"""
Created on Mon May 16 10:47:18 2016

@author: mdurant
"""
import sys
from PyQt4 import QtGui, QtCore
import qtawesome
from pprint import pprint

import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import ROI
import blood_tools

class PlotWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        super(PlotWidget, self).__init__(parent)
        self.r = None
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)
        #self.make_image(np.zeros((10, 10)))
        
        
    def make_image(self, im):
        self.clear_roi()
        self.axes.clear()
        print im
        self.im = self.axes.imshow(im, vmin=np.percentile(im, 5),
                                    vmax=np.percentile(im, 95), cmap='gray')
        #self.im = self.axes.imshow(im)

    def new_roi(self):
        self.r = ROI.new_ROI(self.im)

    def clear_roi(self):
        if self.r is not None: # check if there is an ROI
            self.r.patch.remove()
            self.r.disconnect()
            self.axes.lines.clear()

    
    def get_roi(self):
        return self.r.get_mask()


class MainWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self, parent=None)
        
        button_load = QtGui.QPushButton(qtawesome.icon('fa.folder-open-o'), '')
        self.button_roi = QtGui.QCheckBox("ROI")
        button_run = QtGui.QPushButton(qtawesome.icon('fa.play'), '')
        self.button_roi.setChecked(False)
        self.combo_relax = QtGui.QComboBox()
        self.combo_relax.addItems(['T1', 'T2'])
        self.combo_relax.setCurrentIndex(0)
        
        self.plot_im = PlotWidget(self)
        self.plot_graph = PlotWidget(self)
        
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
    
    def choose_dir(self):
        out = QtGui.QFileDialog.getExistingDirectory(caption='MRI Data Directory')
        
        if out:
            self.directory = out
            #image = make_image_from_directory(out)
            images, attributes = blood_tools.read_dicoms(out,[])
            self.plot_im.make_image(images[0])
        else:
            pass
    
    def start_roi(self, state):
        if state:
            self.plot_im.new_roi()
        else:
            self.plot_im.clear_roi()

    def process_data(self):
        roi_mask = self.plot_im.get_roi()
        relax = self.combo_relax.currentText()
        x, y = make_graph_from_dir_and_mask(self.directory, roi_mask)
        self.plot_graph.axes.clear()
        self.plot_graph.axes.plot(x, y, 'k+')
        x, y, parameters = make_fit_from_data(x, y, relax)
        self.plot_graph.axes.plot(x, y, 'r-')
        self.plot_graph.axes.figure.text(figx, figy, parameters)

def main():
    app = QtGui.QApplication.instance() or QtGui.QApplication([])
    win = MainWindow()
    win.show()
    app.exec_()
    return win

if __name__ == '__main__':
    win = main()