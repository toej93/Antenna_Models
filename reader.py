import csv
import matplotlib.pyplot as pyp
import pandas as pd
import numpy as np
import scipy.constants
import numpy.fft as fft
import cmath

def read_file(finame):
    fi = open(finame, 'r')
    rdr = csv.reader(fi, delimiter=',', skipinitialspace=True)
    table = []
    for row in rdr:
        #print(row)
        freq = float(row[0])*1000
        impedance = (1+float(row[1]))/(1-float(row[1]))
        row = {'freq':freq, 'impedance':impedance}
        table.append(row)
    return pd.DataFrame(table)

def impedance(trans_coeff):
    Z=(trans_coeff+1)/(trans_coeff-1)
    return Z

def effective_height(freq, gain, phase):
    phase = np.unwrap(phase)
    cgain = gain * np.cos(phase*np.pi/180.) + 1j * gain * np.sin(phase*np.pi/180.)
    heff2 = 50./377. * cgain * scipy.constants.c**2 / (np.pi *freq**2 )
    heff = complex_square_root(heff2)
    return heff

wipld=read_file("S11_bicone.csv")
pyp.plot(wipld['freq'], wipld['impedance'])
pyp.show()

