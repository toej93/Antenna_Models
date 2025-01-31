#!/usr/bin/python

import csv
import matplotlib.pyplot as pyp
import pandas as pd
import numpy as np
import scipy.constants
import numpy.fft as fft
import cmath

hPol = "HPol_XFDTD_CurrentMod.txt"
bicone = "ARA_bicone6in_output.txt"
dipole = "ARA_dipoletest1_output.txt"
LPDA = "ARIANNA_LPDA.txt"
LPDA_ice="ARIANNA_LPDA_firn.txt"
vpol_full = "VPol_XFDTD_FullPlate.txt"
vpol_cross = "VPol_XFDTD_Crossbar.txt"
bicone_w="ARA_bicone.txt"


#get_ipython().run_line_magic('matplotlib', 'inline')
pyp.rcParams['font.size']=18
pyp.rcParams['legend.fontsize']=10

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def read_ARA_antenna_model(finame):
    fi = open(finame, 'r')
    rdr = csv.reader(fi, delimiter=' ', skipinitialspace=True)
    freq = 0.
    SWR = 999.
    impedance = 0.
    table = []
    for row in rdr:
    #    print(row)
        if( row[0] == 'freq'):
            freq = float(row[2])
        elif( row[0] == 'SWR'):
            SWR = float(row[2])
        elif( row[0] == 'Impedance'):
            impedance = float(row[2])
        elif( row[0] == 'Theta' ):
            continue
        else:
            theta = float(row[0])
            phi = float(row[2])
            gaindb= float(row[4])
            gain = float(row[6])
            phase = float(row[8])
            row = {'freq':freq, 'SWR':SWR, 'theta':theta, 'phi':phi, 'gaindb':gaindb, 'gain':gain, 'phase':phase, 'impedance':impedance}
            #print(row)
            table.append(row)
    return pd.DataFrame(table)

def read_file(finame):
    fi = open(finame, 'r')
    rdr = csv.reader(fi, delimiter=' ', skipinitialspace=True)
    table = []
    for row in rdr:
    #    print(row)
        freq = float(row[0])
        heff = float(row[1])
        row = {'freq':freq, 'heff':heff}
        table.append(row)
    return pd.DataFrame(table)

        
def read_ARA_antenna_model_updated(finame):
    fi = open(finame, 'r')
    rdr = csv.reader(fi, delimiter=' ', skipinitialspace=True)
    freq = 0.
    SWR = 999.
    table = []
    for row in rdr:
#        print(row)
        if( row[0] == 'freq'):
            freq = float(row[2])
        elif( row[0] == 'SWR'):
            
            SWR = float(row[2])
        elif( row[0] == 'Theta' ):
            continue
        else:
            theta = float(row[0])
            phi = float(row[1])
            gaindb= float(row[2])
            gain = float(row[3])
            phase = float(row[4])
            row = {'freq':freq, 'SWR':SWR, 'theta':theta, 'phi':phi, 'gaindb':gaindb, 'gain':gain, 'phase':phase}
            table.append(row)
    return pd.DataFrame(table)

def complex_square_root(c):
    r = abs(c)
    theta = np.angle(c)
    #unwrap the phases
    theta = np.unwrap(theta)
    return np.sqrt(r)*( np.cos( theta/2 ) + 1j*np.sin(theta/2) )

def effective_height(freq, gain, phase):
    phase = np.unwrap(phase)
    cgain = gain * np.cos(phase*np.pi/180.) + 1j * gain * np.sin(phase*np.pi/180.)
    #heff2 = impedance/377. * cgain * scipy.constants.c**2 / (np.pi *freq**2 )#1.78=n_ice #non-constant Z
    heff2 = 50/377. * cgain * scipy.constants.c**2 / (1.78*np.pi *freq**2)  #constant Z
    heff = 0.5*complex_square_root(heff2)
    return heff

#A_e = np.sqrt(np.array(gains)/(np.pi)*(c/(frequencies*1e6))**2*np.array(self.Re_Z)/(Z_0))


vpol_cross_ = read_ARA_antenna_model_updated(vpol_cross)
hPol_ = read_ARA_antenna_model_updated(hPol)
bicone_ = read_ARA_antenna_model(bicone)
dipole_ = read_ARA_antenna_model(dipole)
LPDA_ = read_ARA_antenna_model(LPDA)
LPDA_ice_ = read_ARA_antenna_model(LPDA_ice)
vpol_full_ = read_ARA_antenna_model_updated(vpol_full)
bicone_w_ = read_ARA_antenna_model(bicone_w)
#wipld=read_file("bicone_wipld_Z.csv")

from matplotlib import cm
from numpy import linspace
"""
start = 0.0
stop = 1.0
number_of_lines= 10
cm_subsection = linspace(start, stop, number_of_lines)
colors = [ cm.jet(x) for x in cm_subsection ]

for i, f in enumerate(drange(300.01, 1000, 100)):
    boresight = tab[tab.phi==90]
    boresight = boresight[boresight.freq==f]
    color = colors[i]
    heff = effective_height(boresight.freq*1e6, boresight.gain, boresight.phase)
    pyp.plot(boresight.theta, abs(heff), color=color, label=str(f))
    
pyp.title("ARA Bicone Effective Height")
pyp.xlabel("Frequency (MHz)")
pyp.ylabel("Effective Height (m)")
pyp.ylim(0,0.3)
#pyp.figure(2)
#pyp.title("Boresight Effective Height")
#pyp.xlabel("Frequency (MHz)")
#pyp.ylabel("Effective Height (m)")
pyp.legend(title='Freq (MHz)', loc=[1.,0.001])
pyp.show()

"""
"""
start = 0.0
stop = 1.0
number_of_lines= 10
cm_subsection = linspace(start, stop, number_of_lines)
colors = [ cm.jet(x) for x in cm_subsection ]
for i,t in enumerate(range(0, 100, 10)):
    ant = tab[tab.phi==90.]
    ant = ant[ant.theta==t]
    color = colors[i]
    heff = effective_height(ant.freq*1e6, ant.gain, ant.phase)
    pyp.plot(ant.freq, abs(heff), color=color, label=str(t))
pyp.title("ARA Bicone (WIPL-D) Effective Height ($\\phi=90$)")
pyp.xlabel("Frequency (MHz)")
pyp.ylabel("Effective Height (m)")
pyp.ylim(0,1)
pyp.legend(title='$\\theta (deg)$', loc=[1.05,0.001])
pyp.show()    
"""

#Plot effective height
boresight_vc = vpol_cross_[vpol_cross_.phi==90.]
boresight_vc = boresight_vc[boresight_vc.theta==90.]
boresight_vf = vpol_full_[vpol_full_.phi==90.]
boresight_vf = boresight_vf[boresight_vf.theta==90.]
boresight_v = bicone_[bicone_.phi==90.]
boresight_v = boresight_v[boresight_v.theta==90.]
boresight_vw = bicone_w_[bicone_w_.phi==90.]
boresight_vw = boresight_vw[boresight_vw.theta==90.]

"""

heff_vc = effective_height(boresight_vc.freq*1e6, boresight_vc.gain, boresight_vc.phase)
heff_vf = effective_height(boresight_vf.freq*1e6, boresight_vf.gain, boresight_vf.phase)
heff_v = effective_height(boresight_v.freq*1e6, boresight_v.gain, boresight_v.phase)
heff_vw = effective_height(boresight_vw.freq*1e6, boresight_vw.gain, boresight_vw.phase)

#p = np.poly1d(np.polyfit(boresight_vc.freq, abs(heff_vc), 8))


pyp.plot(boresight_v.freq, abs(heff_v),'g--', linewidth=2,label='Old model')
#pyp.plot(boresight_vf.freq, abs(heff_vf),'b:',linewidth=4, label='XFDTD_FullBar')
pyp.plot(boresight_vc.freq, abs(heff_vc),'r-.',linewidth=4, label='New model (Chiba U)')
#pyp.plot(boresight_vw.freq, abs(heff_vw),'k',linewidth=2, label='WIPLD (constant $Z$)')
#pyp.plot(wipld['freq']*1000, wipld['heff'],'m',linewidth=2, label='WIPLD (non-constant $Z$)')
#pyp.plot(p)
pyp.legend(loc='upper center')
pyp.title("ARA Bicone Effective Height")
pyp.xlabel("Frequency (MHz)")
pyp.ylabel("Heff (m)")
#pyp.ylim(0,1)
pyp.grid()
#pyp.xlim(50,1100)
pyp.show()    
"""
#Plot gains
tab = LPDA_[LPDA_.phi==90.]
tab = tab[tab.theta==90]
tab2 = LPDA_ice_[LPDA_ice_.phi==90.]
tab2 = tab2[tab2.theta==90]
heff = effective_height(tab2.freq*1e6, tab2.gain, tab2.phase)
#heff_ice = effective_height(tab2.freq*1e6, tab2.gain, tab2.phase, tab2.impedance)
#pyp.plot(.freq, abs(heff),'r',linewidth=2,label='Yue interpolation')
pyp.plot(tab2.freq, abs(heff), 'b--',linewidth=2,label='Heff (abs(S11)')
pyp.title("Bicone gain")
pyp.xlabel("Frequency (MHz)")
pyp.ylabel("Gain [arb]")
#pyp.ylim(0,1)
pyp.xlim(0,1000)
wipld=read_file("Heff_LPDA_boresight.txt")
pyp.plot(wipld['freq'], wipld['heff'],'g', linewidth=2,label='From Christian')
pyp.grid()
pyp.legend(loc='upper center')
pyp.show()    
