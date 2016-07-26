"""
Created on Thu Jul 21 14:15:11 2016

@author: Josh

This is an implementation of the math in Sharon Portney's manuscript.

Notes:
1. The physical range of Hct with measurement error is: -0.05 < Hct < 1.05
2. The physical range of sO2 with measurement error is: -0.1 < sO2 < 1.1

These uncertainties were determined by using measuremenents with up to 5% error
in a Monte Carlo simulation.

There is an overlap in the solution space for certain {R1, R2} vals. When 
multiple solutions arise, simply choosing the one that has a real part within 
the acceptable range is sufficient.

A normal fetal sO2 is 0.85, so in some cases it may be necessary to rule out 
solutions with sO2 > 0.95
"""

# std lib imports
from __future__ import division
import os
import traceback
import types
import math
# anaconda module imports
from PyQt4 import QtGui, QtCore
import numpy as np
import csv
from inspect import getsourcefile
from pprint import pprint
from functools import wraps

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

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

class MainWindow(QtGui.QWidget):
    @QTSlotExceptionRationalizer("bool")
    def __init__(self):
        # draw the interface
        self.init_gui()
        self.input_type_changed()

    @QTSlotExceptionRationalizer("bool")
    def init_gui(self):
        QtGui.QWidget.__init__(self, parent=None)
        # button for opening dicom directory
        self.t1_input = QtGui.QLineEdit()
        self.t1_input.setValidator(QtGui.QDoubleValidator())
        self.t1_input.setMaxLength(6)
        
        self.t2_input = QtGui.QLineEdit()
        self.t2_input.setValidator(QtGui.QDoubleValidator())
        self.t2_input.setMaxLength(6)
        
        self.tau_180_input = QtGui.QLineEdit()
        self.tau_180_input.setValidator(QtGui.QDoubleValidator())
        self.tau_180_input.setMaxLength(6)
        
        self.sO2_input = QtGui.QLineEdit()
        self.sO2_input.setValidator(QtGui.QDoubleValidator())
        self.sO2_input.setMaxLength(6)
        
        self.Hct_input = QtGui.QLineEdit()
        self.Hct_input.setValidator(QtGui.QDoubleValidator())
        self.Hct_input.setMaxLength(6)
        
        self.sO2_output = QtGui.QLineEdit()
        self.sO2_output.setValidator(QtGui.QDoubleValidator())
        
        self.Hct_output = QtGui.QLineEdit()
        self.Hct_output.setValidator(QtGui.QDoubleValidator())
        
        self.base_field_strength_selector = QtGui.QComboBox()   
        self.base_field_strength_selector.addItems(['1.5T','3T'])
        self.base_field_strength_selector.setCurrentIndex(0)
        
        self.case_type_selector = QtGui.QComboBox()   
        self.case_type_selector.addItems(['Fetal','Adult'])
        self.case_type_selector.setCurrentIndex(0)
        
        self.input_type_selector = QtGui.QComboBox()   
        self.input_type_selector.addItems(['T1_T2', 'T1_sO2', 'T2_Hct'])
        self.input_type_selector.setCurrentIndex(0)
        
        self.button_solve = QtGui.QPushButton('SOLVE') 
        
        layout_main = QtGui.QGridLayout()
        layout_main.addWidget(QtGui.QLabel('INPUTS'), 1, 1, 1, 2)
        layout_main.addWidget(QtGui.QLabel('Input Type:'), 2, 1)
        layout_main.addWidget(self.input_type_selector, 2, 2)
        
        layout_main.addWidget(QtGui.QLabel('T1 (ms)'), 3, 1)
        layout_main.addWidget(self.t1_input, 3, 2)
        layout_main.addWidget(QtGui.QLabel('T2 (ms)'), 4, 1)
        layout_main.addWidget(self.t2_input, 4, 2)
        layout_main.addWidget(QtGui.QLabel('tau_180 (ms)'), 5, 1)
        layout_main.addWidget(self.tau_180_input, 5, 2)
        layout_main.addWidget(QtGui.QLabel('sO2 (decimal)'), 6, 1)
        layout_main.addWidget(self.sO2_input, 6, 2)
        layout_main.addWidget(QtGui.QLabel('Hct (decimal)'), 7, 1)
        layout_main.addWidget(self.Hct_input, 7, 2)
        layout_main.addWidget(QtGui.QLabel('Field Strength:'), 8, 1)
        layout_main.addWidget(self.base_field_strength_selector, 8, 2)
        layout_main.addWidget(QtGui.QLabel('Vessel Type:'), 9, 1)
        layout_main.addWidget(self.case_type_selector, 9, 2)
        
        layout_main.addWidget(self.button_solve, 10,1,1,4)
        layout_main.addWidget(QtGui.QLabel("Created by Josh Bradshaw, using Sharon Portney's manuscript.\nCite: \nMIT Licensed"), 11, 1, 1, 4)
        
        layout_main.addWidget(QtGui.QLabel('OUTPUTS'), 1, 3, 1, 2)
        layout_main.addWidget(QtGui.QLabel('sO2 (decimal)'), 2, 3)
        layout_main.addWidget(self.sO2_output, 2, 4)
        layout_main.addWidget(QtGui.QLabel('Hct (decimal)'), 3, 3)
        layout_main.addWidget(self.Hct_output, 3, 4)
        
        self.setLayout(layout_main)
        self.button_solve.pressed.connect(self.solve)
        self.input_type_selector.currentIndexChanged.connect(self.input_type_changed)
    
    @QTSlotExceptionRationalizer("bool")
    def solve(self):
        input_type = self.input_type_selector.currentText()
        field_strength = self.base_field_strength_selector.currentText()
        case_type = self.case_type_selector.currentText()
        model_params = get_model_parameters(field_strength, case_type)
        
        t1_text = self.t1_input.text().strip()
        t2_text = self.t2_input.text().strip()
        sO2_text = self.sO2_input.text().strip()
        Hct_text = self.Hct_input.text().strip()
        tau_180_text = self.tau_180_input.text().strip() 
        
        if input_type == 'T1_T2':
            t1, t2, tau_180 = text_to_num(t1_text, t2_text, tau_180_text)
            # ms to s conversion            
            tau_180 = tau_180 / 1000
            Hct_vals, SO2_vals = Hct_sO2_from_T1_T2(model_params, t1, t2, tau_180)
            hct_output_str = " | ".join('{:.3f}'.format(val) for val in Hct_vals)         
            sO2_output_str = " | ".join('{:.3f}'.format(val) for val in SO2_vals)
            self.Hct_output.setText(hct_output_str)
            self.sO2_output.setText(sO2_output_str)
        elif input_type == 'T1_sO2':
            t1, sO2 = text_to_num(t1_text, sO2_text)
            # ms to s conversion       
            Hct = Hct_from_T1_sO2(model_params, t1, sO2)
            hct_output_str = '{}'.format(Hct)
            self.Hct_output.setText(hct_output_str)
        elif input_type == 'T2_Hct':
            t2, Hct, tau_180 = text_to_num(t2_text, Hct_text, tau_180_text)
            tau_180 = tau_180 / 1000            
            # ms to s conversion
            sO2_vals = sO2_from_T2_Hct(model_params, t2, Hct, tau_180)
            sO2_output_str = " | ".join('{:.3f}'.format(val) for val in sO2_vals)
            self.sO2_output.setText(sO2_output_str)
        else:
            raise ValueError('Invalid input type')
    
    @QTSlotExceptionRationalizer("bool")
    def input_type_changed(self, *e):
        input_type = self.input_type_selector.currentText()
        
        if input_type == 'T1_T2':
            self.t1_input.setEnabled(True)
            self.t2_input.setEnabled(True)
            self.sO2_input.setEnabled(False)
            self.sO2_input.setText('')
            self.Hct_input.setEnabled(False)
            self.Hct_input.setText('')
            self.tau_180_input.setEnabled(True)
            self.Hct_output.setEnabled(True)
            self.sO2_output.setEnabled(True)
        elif input_type == 'T1_sO2':
            self.t1_input.setEnabled(True)
            self.t2_input.setEnabled(False)
            self.t2_input.setText('')
            self.sO2_input.setEnabled(True)
            self.Hct_input.setEnabled(False)
            self.Hct_input.setText('')
            self.tau_180_input.setEnabled(False)
            self.tau_180_input.setText('')
            self.Hct_output.setEnabled(True)
            self.sO2_output.setEnabled(False)
            self.sO2_output.setText('')
        elif input_type == 'T2_Hct':
            self.t1_input.setEnabled(False)
            self.t1_input.setText('')
            self.t2_input.setEnabled(True)
            self.sO2_input.setEnabled(False)
            self.sO2_input.setText('')
            self.Hct_input.setEnabled(True)
            self.tau_180_input.setEnabled(True)
            self.Hct_output.setEnabled(False)
            self.Hct_output.setText('')
            self.sO2_output.setEnabled(True)
        else:
            raise ValueError('Invalid input type')
            
def text_to_num(*values):
    # todo: catch exceptions and display warning message if value out of range
    converted_vals = []
    for val in values:
        converted_vals.append(float(val))
        
    return converted_vals
            
def get_model_parameters(field_strength, vessel_type):
    GAMMA = 42.576 # Hz/Tesla    
    f_s = float(field_strength.replace('T', ''))
    expected_fieldnames = ['R1_plas','R1_ery_0','r1_prime_dHB','R2_plas',
    'R2_dia_plus_R2_oxy','R2_deoxy_minus_R2_oxy','w_dia_plus_w_oxy',
    'w_deo_minus_w_oxy','tau']        

    dir_name = os.path.dirname(getsourcefile(lambda:0))
    parameter_fp = os.path.join(dir_name, 'fit_parameters', field_strength, vessel_type.lower() + '.csv')
    
    parameters = {}        
    
    with open(parameter_fp, 'r') as f:
        reader = csv.reader(f)

        for row in reader:
            p_name = row[0].strip()
            try:
                p_val = float(row[1].strip())
            except:
                err_msg = 'Parameter: {} is not properly defined in file: {}'.format(p_name, parameter_fp)
                raise ValueError(err_msg)
            parameters[p_name] = p_val
            
    # sanity check: TODO catch this exception in the GUI
    for exp_par in expected_fieldnames:
        if not exp_par in parameters:
            err_msg = 'Parameter: {} is not properly defined in file: {}'.format(exp_par, parameter_fp)
            raise ValueError(err_msg)
    
    # some not entirely obvious unit conversions
    # convert wdia + woxy and wdeo - woxy from ppm to rads/sec
    parameters['w_dia_plus_w_oxy'] = 2 * math.pi * GAMMA * parameters['w_dia_plus_w_oxy'] * f_s
    parameters['w_deo_minus_w_oxy'] = 2 * math.pi * GAMMA * parameters['w_deo_minus_w_oxy'] * f_s
    parameters['tau'] = parameters['tau'] / 1000
    
    pprint(parameters)
    return parameters
        
def intermediate_constants(parameters, tau_180):
    p = Bunch(parameters)
    
    mu = p.tau * (1 - 2*p.tau*math.tanh(tau_180/(2*p.tau))/tau_180)
    K0 = p.R2_plas
    K1 = p.R2_dia_plus_R2_oxy + p.R2_deoxy_minus_R2_oxy + mu*math.pow(p.w_dia_plus_w_oxy + p.w_deo_minus_w_oxy, 2)
    K2 = -mu*math.pow((p.w_dia_plus_w_oxy + p.w_deo_minus_w_oxy), 2)
    K3 = -p.R2_deoxy_minus_R2_oxy - 2*mu*p.w_deo_minus_w_oxy*(p.w_dia_plus_w_oxy + p.w_deo_minus_w_oxy)
    K4 = 2*mu*p.w_deo_minus_w_oxy * (p.w_dia_plus_w_oxy + p.w_deo_minus_w_oxy)
    K5 = mu*math.pow(p.w_deo_minus_w_oxy, 2)
    M0 = p.R1_plas
    M1 = p.R1_ery_0 - p.R1_plas + p.r1_prime_dHB
    M2 = -p.r1_prime_dHB
    return (K0, K1, K2, K3, K4, K5, M0, M1, M2)

def SO2_from_T1_Hct(parameters, T1, Hct):
    # equation 9
    p = Bunch(parameters)
    R1 = 1000/T1

    sO2 = (R1 - p.R1_plas - Hct * (p.R1_ery_0 - p.R1_plas + p.r1_prime_dHB) ) / (-p.r1_prime_dHB*Hct)
    return sO2
    
def sO2_from_T2_Hct(parameters, T2, Hct, tau_180):
    # equation 11
    R2 = 1000/T2
    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(parameters, tau_180)    
    
    polynomial_coeff = [K5 *Hct - K5*math.pow(Hct, 2), K3*Hct + K4 * math.pow(Hct, 2), K0 + K1 * Hct + K2 * math.pow(Hct, 2) - R2]
    sO2_roots = np.real(np.roots(polynomial_coeff))
    return sO2_roots
    
def Hct_from_T1_sO2(parameters, T1, sO2):
    # equation 10
    R1 = 1000/T1
    p = Bunch(parameters)
    
    Hct = (R1 - p.R1_plas) / (p.R1_ery_0 - p.R1_plas + p.r1_prime_dHB*(1-sO2))
    return Hct

def Hct_sO2_from_T1_T2(parameters, T1, T2, tau_180):
    # using equations 7 and 9
    R1 = 1000/T1 # 1/s (to be set)
    R2 = 1000/T2 # 1/s (to be set)

    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(parameters, tau_180)
    
    A1 = K2*M2**2 - K4*M1*M2 - K5*M1**2
    B1 = K1*M2**2 - K3*M1*M2 + K4*(R1*M2 - M0*M2) + K5*(M1**2 + 2*R1*M1 - 2*M0*M1)
    C1 = K0*M2**2 + K3*(R1*M2 - M0*M2) + K5*(2*M0*M1 - 2*M1*R1 - (R1-M0)**2) - R2*M2**2
    D1 = K5*(R1 - M0)**2
    
    # equation 7
    polynomial_coefficients = [A1, B1, C1, D1]    
    roots = np.roots(polynomial_coefficients)
    Hct_vals = np.real(roots)
    
    SO2_vals = []
    for hct_root in Hct_vals:
        SO2_vals.append(SO2_from_T1_Hct(parameters, T1, hct_root))
    
    return Hct_vals, SO2_vals

def alt_Hct_sO2_from_T1_T2(parameters, T1, T2, tau_180):
    # using equations 8 and 10
    R1 = 1000/T1 # 1/s (to be set)
    R2 = 1000/T2 # 1/s (to be set)
    
    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(parameters, tau_180)
    
    A2 = K5*(R1*M2 - M0*M2)
    B2 = K5*(R1*M1 - M0*M1 - (M0-R1)**2) + K3*(R1*M2 - M0*M2) + K0*M2**2 - R2*M2**2
    C2 = K4*(M0 - R1)**2 + K3*(R1*M1 - M0*M1) + K1*(R1*M2 - M0*M2) + 2*K0*M1*M2 - 2*R2*M1*M2
    D2 = K2*(M0 - R1)**2 + K1*(R1*M1 - M0*M1) + K0*M1**2 - R2*M1**2
    
    # equation 7
    polynomial_coefficients = [A2, B2, C2, D2]  
    roots = np.roots(polynomial_coefficients)
    SO2_vals = np.real(roots)
    
    Hct_vals = []
    for so2_root in SO2_vals:
        Hct_vals.append(Hct_from_T1_sO2(parameters, T1, so2_root))
    return Hct_vals, SO2_vals
        
def main():
    app = QtGui.QApplication.instance() or QtGui.QApplication([])
    win = MainWindow()
    win.setWindowTitle('sO2 and Hematocrit Calculator')
    win.show()
    app.exec_()
    return win

if __name__ == '__main__':
    win = main()