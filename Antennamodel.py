import ROOT
from array import array
import sys
import numpy as np
from scipy.interpolate import interp1d
import geometryHelper


class AntennaModel(object):

    def __init__(self,f = "WIPLD_antennamodel_firn.root",nfreq = 110):
    
        """ Use the antenna model as provided by WIPL-D,
        Antenna (in simulations) is assumed to be aligned with x-axis, boresight towards positive x, tines are
        in the x-y plane. 
        Phi is the angle in antenna model in xy-plane, phi = 0 is boresight and mathematically
        positive. Theta is angle to z-axis. Range [-90,90], antenna at theta = 0 
        The angles azimuth and zenith are ARIANNA convention (normal Cartesian, x East) and 
        determine the arrival direction. 
        All angles as inputs and outputs are expected in degrees.
        The input is a ROOT tree that has been filled from the WIPL-D simulations using makeROOTTree.py
        All frequencies are provided in MHz.
        This class provides custom functions to rotate antenna model into place for given antenna
        orientation. 
        """
        
        """First tree in file contains general values such as gain and currents. """

        self.nt = ROOT.TChain("AntTree")
        self.nt.Add(f)

        self.N = array('i',[0])
        self.thetas = array('f',[0])
        self.phis = array('f',[0])

        self.frequencies = array('f',nfreq*[0.])
       # self.frequencies = array('f',[0])        
        self.gains       = array('d',nfreq*[0.])
        self.Re_phi      = array('d',nfreq*[0.])
        self.Im_phi      = array('d',nfreq*[0.])
        self.Re_theta    = array('d',nfreq*[0.])
        self.Im_theta    = array('d',nfreq*[0.])

        self.nt.SetBranchAddress("N", self.N)
        self.nt.SetBranchAddress("thetas", self.thetas)
        self.nt.SetBranchAddress("phis", self.phis)
        self.nt.SetBranchAddress("gains", self.gains)
        self.nt.SetBranchAddress("frequencies",self.frequencies)
        self.nt.SetBranchAddress("Re_phi", self.Re_phi)
        self.nt.SetBranchAddress("Im_phi", self.Im_phi)
        self.nt.SetBranchAddress("Re_theta",self.Re_theta)
        self.nt.SetBranchAddress("Im_theta",self.Im_theta)

        self.nt.BuildIndex("thetas", "phis")
        
        """ Second tree in file contains only Z as function of frequency, 
        it is the same for all phi and theta."""
        
        self.ntt = ROOT.TChain("ZTree")
        self.ntt.Add(f)
        
        self.NN = array('i',[0])
        self.Re_Z  = array('d',nfreq*[0.])
        self.Im_Z  = array('d',nfreq*[0.])
        
        self.ntt.SetBranchAddress("N", self.NN)
        self.ntt.SetBranchAddress("Re_Z", self.Re_Z)
        self.ntt.SetBranchAddress("Im_Z", self.Im_Z)
        
        self.ntt.GetEntry(0)


    def LoadData(self, theta, phi):
        """Loading data into variables from tree. Perform checking of chosen directions. """
        theta = int(round(theta))
        
        phi = int(round(phi))
        n = self.nt.GetEntryNumberWithIndex(theta,phi)
        if (n>-1):
            self.nt.GetEntry(n)    
            if not (ROOT.TSnMath.AlmostEqualRelativeDbl(theta, self.thetas[0]) | 
                    ROOT.TSnMath.AlmostEqualRelativeDbl(phi, self.phis[0])):
                sys.exit()    
        else:
            sys.exit() 
        
        
        
    def CheckBoundaries(self,theta,phi):
        """ Wrapping of angles according to angel convention of WIPL-D. """
        if theta > 90:
            theta -= 180
        if theta < -90:
            theta += 180
        if phi < 0:
            phi += 360
        if phi > 360:
            phi -= 360
        
        return theta, phi

    def GetVector(self,zenith,azimuth):
        """Initialize vector for rotation."""
        v = ROOT.TVector3(1,0,0)
        v.SetPhi(np.deg2rad(azimuth))
        v.SetTheta(np.deg2rad(zenith))
        return v  
        
    def Cartesian2WPLD(self,v):
        """ WIPL-D theta range different than Cartesian"""         
        theta = 90. - np.rad2deg(v.Theta())
        phi = np.rad2deg(v.Phi())  
        theta, phi = self.CheckBoundaries(theta,phi)
        return theta, phi  
        
    def WPLD2Cartesian(self,v):
        """ WIPL-D theta range different than Cartesian"""         
        theta = 90. - np.rad2deg(v.Theta())
        phi = np.rad2deg(v.Phi())  
        theta, phi = self.CheckBoundaries(theta,phi)
        return theta, phi     
    
    def RotateEfieldXYZ(self,vec,antenna):
        """ Use antenna tag as switch to rotate antenna response.
        Currently implemented: 
        DEW (downward East West aligned) 
        DNS (downward North South aligned)
        UPN (upward pointing North, tines East West aligned, station 32 )
        UPW (upward pointing West, tines North South aligned, station 32)
        HCRH (HCR Tower, horizontal pol)
        HCRV (upward pointing West, tines North South aligned)
        UEW (upward East West aligned)
        UNS (upward North South aligned) """
        
        if antenna == "DEW":
            m1 = geometryHelper.RotY(np.pi/2)
            vec_n = np.dot(m1, vec)
            m2 = geometryHelper.RotZ(-1*np.pi/2)
            vec_n = np.dot(m2, vec_n)
        elif antenna == "DNS":    
            m = geometryHelper.RotY(np.pi/2)
            vec_n = np.dot(m, vec)
        elif antenna == "UPN":
            m1 = geometryHelper.RotZ(np.pi/2)
            vec_n = np.dot(m1, vec)
            m2 = geometryHelper.RotX(np.pi/4.)
            vec_n = np.dot(m2, vec_n)
        elif antenna == "UPW":
            m = geometryHelper.RotY(np.pi+np.pi/4.)
            vec_n = np.dot(m, vec)
        elif antenna == "HCRH":
            m = geometryHelper.RotZ(np.pi/2)
            vec_n = np.dot(m, vec)
        elif antenna == "HCRV":
            m1 = geometryHelper.RotY(np.pi/2)
            vec_n = np.dot(m1, vec)
            m2 = geometryHelper.RotZ(np.pi/2.)
            vec_n = np.dot(m2, vec_n)    
        elif antenna == "UEW":
            m1 = geometryHelper.RotY(np.pi/2)
            vec_n = np.dot(m1, vec)
            m2 = geometryHelper.RotZ(-1*np.pi/2)
            vec_n = np.dot(m2, vec_n)
        elif antenna == "UNS":
            m = geometryHelper.RotY(-np.pi/2)
            vec_n = np.dot(m, vec)        
            
        else:
            vec_n = 0
        return vec_n
        
    def GetThetaAndPhi(self,zenith,azimuth,antenna):
        """ Use antenna tag as switch to get theta and phi to call antenna model"""
        if antenna == "DEW":
            theta, phi = self.GetThetaAndPhiForDownwardEWaligned(zenith,azimuth)
        elif antenna == "DNS":    
            theta, phi = self.GetThetaAndPhiForDownwardNSaligned(zenith,azimuth)
        elif antenna == "UPN":
            theta, phi = self.GetThetaAndPhiForUpwardNorth(zenith,azimuth)
        elif antenna == "UPW":
            theta, phi = self.GetThetaAndPhiForUpwardWest(zenith,azimuth)
        elif antenna == "UPWF":
            theta, phi = self.GetThetaAndPhiForUpwardWestField(zenith,azimuth)
        elif antenna == "UPNF":
            theta, phi = self.GetThetaAndPhiForUpwardNorthField(zenith,azimuth)
        elif antenna == "HCRH":
            theta, phi = self.GetThetaAndPhiForHCRHpol(zenith,azimuth)
        elif antenna == "HCRV":
            theta, phi = self.GetThetaAndPhiForHCRVpol(zenith,azimuth)
        elif antenna == "UEW":
            theta, phi = self.GetThetaAndPhiForUpwardEWaligned(zenith,azimuth)
        elif antenna == "UNS":    
            theta, phi = self.GetThetaAndPhiForUpwardNSaligned(zenith,azimuth)

        else:
            return theta, p1hi
    
    """ Add here possible new antenna configurations.
        Define function that turns the WIPL-D standard antenna (nose pointing in positive x, tines in xy-plane)
        to the desired location by using rotations."""
        
    def GetThetaAndPhiFlexible(self,zenith,azimuth,azi_orientation=90,zen_orientation=100):
        """ LPDA for various rotations part from straight down or straight up"""
        if ((zen_orientation == 0) or (zen_orientation == 180)):
            return -1, -1
        v = self.GetVector(zenith,azimuth)  
        v.RotateZ(np.deg2rad(-1*azi_orientation))
        v.RotateY(np.deg2rad((90-zen_orientation)))
        theta, phi = self.Cartesian2WPLD(v)
        return theta, phi
        
    def GetThetaAndPhiForDownwardEWaligned(self,zenith,azimuth):
        """ LPDA  pointing downward with tines EW aligned """
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateZ(np.pi/2)
        v.RotateY(-1*np.pi/2)
        theta, phi = self.Cartesian2WPLD(v)
        return theta, phi
        
    def GetThetaAndPhiForDownwardNSaligned(self,zenith,azimuth):
        """ LPDA  pointing downward with tines NS aligned """
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateY(-1*np.pi/2)
        theta, phi = self.Cartesian2WPLD(v)
        return theta, phi
            
    def GetThetaAndPhiForUpwardNorth(self,zenith,azimuth):
        """ LPDA  pointing upwards North at 45 degrees with tines EW aligned """
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateZ(-1*np.pi/2.) #=/- here?
        v.RotateY(np.pi/4.)
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
        
    def GetThetaAndPhiForUpwardWest(self,zenith,azimuth):
        """ LPDA  pointing upward west at 45 degrees with tines NS aligned """
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateY(-1*(np.pi+np.pi/4.))
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi            
        
    def GetThetaAndPhiForUpwardNorthField(self,zenith,azimuth):
        """ LPDA  pointing upwards North at 40 degrees with tines EW aligned (slightly differently aligned)"""
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateZ(-1*np.pi/2.) 
        v.RotateY(np.deg2rad(40))
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
        
    def GetThetaAndPhiForUpwardWestField(self,zenith,azimuth):
        """ LPDA  pointing upward west at 40 degrees with tines NS aligned (slightly differently aligned)"""
        v = self.GetVector(zenith,azimuth)
        #Rotate to antenna orientation
        v.RotateY(-1*(np.pi+np.deg2rad(40)))
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi   
        
    def GetThetaAndPhiForHCRVpol(self,zenith,azimuth):
        """ LPDA mounted on HCR tower towards north, vertical polarization """
        v = self.GetVector(zenith,azimuth)
        v.RotateX(-1*np.pi/2.)
        v.RotateZ(-1*np.pi/2.)
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
        
    def GetThetaAndPhiForHCRHpol(self,zenith,azimuth):
        """ LPDA mounted on HCR tower towards north, horizontal polarization """
        v = self.GetVector(zenith,azimuth)
        v.RotateZ(-1*np.pi/2.)
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
    
    def GetThetaAndPhiForUpwardEWaligned(self,zenith,azimuth):
        v = self.GetVector(zenith,azimuth)
        v.RotateZ(np.pi/2)
        v.RotateY(np.pi/2)
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
        
    def GetThetaAndPhiForUpwardNSaligned(self,zenith,azimuth):
        v = self.GetVector(zenith,azimuth)
        v.RotateY(np.pi/2)
        theta, phi = self.Cartesian2WPLD(v)  
        return theta, phi
        
    def GetFrequencies(self,theta,phi):
        """ Providing all frequencies available """
        self.LoadData(theta,phi)
        #Convert frequencies to MHz as default unit
        f = np.array(self.frequencies)*1e3
        return f

    def GetGains(self,theta,phi):
        """ Providing absolute gain (Not dB) """
        self.LoadData(theta,phi)
        return np.array(self.gains)
        
    def InterpolateGains(self,theta,phi,frequencies):
        """ Provide interpolated gains for given frequencies """
        f =  self.GetFrequencies(theta,phi)
        g = self.GetGains(theta, phi)

        if np.min(frequencies) < np.min(f):
            round(np.min(f)),round(np.max(f))
            f[0] = np.min(frequencies)
            
        if np.max(frequencies) > np.max(f):
            round(np.min(f)),round(np.max(f))) 
            f[-1] = np.max(frequencies)
            
        F_Gains = interp1d(f, g)
        Gains = F_Gains(frequencies)
        
        return Gains
    

    def GetEffectiveHeight(self,gains,frequencies,c=2.99792e8,Z_0=119.9169*np.pi):
        """ Antenna gain converted to Antenna effective height, see also thesis Kamlesh """
        A_e = np.sqrt(np.array(gains)/(np.pi)*(c/(frequencies*1e6))**2*np.array(self.Re_Z)/(Z_0))
        return A_e
        
    def GetVoltages(self,theta,phi):
        """ WIPL-D output of complex voltages for E_theta and E_phi """
        self.LoadData(theta,phi)
        return   np.array(self.Re_theta), np.array(self.Im_theta), np.array(self.Re_phi), np.array(self.Im_phi)
    
    
    def InterpolateVoltages(self,theta,phi,frequencies):
        """ Provide interpolated voltages for given frequencies """
        f =  self.GetFrequencies(theta,phi)
        re_t, im_t, re_p, im_p  = self.GetVoltages(theta,phi)
        
        if np.min(frequencies) < np.min(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f)))
            f[0] = np.min(frequencies)
            
        if np.max(frequencies) > np.max(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f))) 
            f[-1] = np.max(frequencies)         
    
        complex_phi     = re_p + 1j*im_p
        complex_theta   = re_t + 1j*im_t
    
        Phase_phi   = interp1d(f, np.unwrap(np.angle(complex_phi)))
        Phase_theta = interp1d(f, np.unwrap(np.angle(complex_theta)))
    
        Abs_phi = interp1d(f, np.abs(complex_phi))
        Abs_theta = interp1d(f, np.abs(complex_theta)) 
        
        V_theta = Abs_theta(frequencies)*np.exp(1j*Phase_theta(frequencies))  
        V_phi   = Abs_phi(frequencies)*np.exp(1j*Phase_phi(frequencies))
        
        return  V_theta, V_phi

    def GetElectricFieldFactor(self,theta,phi, frequencies,c=2.997e8,Z_0=119.9169*np.pi):
        """ Provide interpolated factors for e-field multiplication. 
        This is the quantity needed for neutrino/CR analysis."""
        wavelength = c / (frequencies*1e6) #Frequencies as default in MHz
        f =  self.GetFrequencies(theta,phi)
        
        if np.min(frequencies) < np.min(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f)))
            f[0] = np.min(frequencies)
            
        if np.max(frequencies) > np.max(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f))) 
            f[-1] = np.max(frequencies)

        V_phi, V_theta = self.InterpolateVoltages(theta,phi,frequencies)
        
        complex_Z = np.array(self.Re_Z) + 1j*np.array(self.Im_Z)
        
        Phase_Z = interp1d(f, np.unwrap(np.angle(complex_Z)))
        Abs_Z = interp1d(f, np.abs(complex_Z))
        
        Z = Abs_Z(frequencies)*np.exp(1j*Phase_Z(frequencies))
                
        E_phi   = ( 2*wavelength * Z* V_phi) / (Z_0) / 1j
        E_theta = ( 2*wavelength * Z* V_theta) / (Z_0) / 1j
    
    def GetElectricFieldFactorAnt(self,zenith,azimuth, frequencies,antenna, c=2.997e8,Z_0=119.9169*np.pi):
        """ Provide interpolated factors for e-field multiplication.
        This is the quantity needed for neutrino/CR analysis."""
        wavelength = c / (frequencies*1e6) #Frequencies as default in MHz
        
        theta, phi = self.GetThetaAndPhi(zenith,azimuth,antenna)
        
        f =  self.GetFrequencies(theta,phi)
        
        if np.min(frequencies) < np.min(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f)))
            f[0] = np.min(frequencies)
            
        if np.max(frequencies) > np.max(f):
            if DEBUG:
                round(np.min(f)),round(np.max(f))) 
            f[-1] = np.max(frequencies)

        V_phi, V_theta = self.InterpolateVoltages(theta,phi,frequencies)
        
        complex_Z = np.array(self.Re_Z) + 1j*np.array(self.Im_Z)
        
        Phase_Z = interp1d(f, np.unwrap(np.angle(complex_Z)))
        Abs_Z = interp1d(f, np.abs(complex_Z))
        
        Z = Abs_Z(frequencies)*np.exp(1j*Phase_Z(frequencies))
                
        E_phi   = ( 2*wavelength * Z* V_phi) / (Z_0) / 1j
        E_theta = ( 2*wavelength * Z* V_theta) / (Z_0) / 1j
        
        E_xyz =  geometryHelper.onsky2ari(vec=np.array([E_theta,E_phi,np.zeros(E_theta.shape[0])]),az=np.deg2rad(phi),zen=np.deg2rad(90-theta))
        
        E_xyz_new = self.RotateEfieldXYZ(E_xyz,antenna)
 
        E_onsky = geometryHelper.ari2onsky(vec=E_xyz_new,az = np.deg2rad(azimuth),zen=np.deg2rad(zenith))
        
        E_theta = E_onsky[0]
        E_phi = E_onsky[1]
        return  E_theta, E_phi
    
##########################################################
#### Everything below here is debug output ###############
#### Might be useful to clean up in later revision #######
##########################################################

DEBUG = True
plot = 3



if DEBUG:        
    import matplotlib.pyplot as plt    
    from mpl_toolkits.mplot3d import Axes3D
    
    
  #  if plot == 9:
        
   #     a = AntennaModel("WIPLD_antennamodel_air.root")
        
        #zenith = 
    
    if plot == 8:
        a = AntennaModel("WIPLD_antennamodel_new_firn.root",nfreq=219)
        
        
        
    if plot == 7:
        zenith = 74.9
        azimuth = 47.29
        antenna='UPW'        
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        
        azimuths = np.arange(0,360,2)
        plt.figure()
        
#         for inx in xrange(len(frequencies)):
        
        t,p = a.GetThetaAndPhi(zenith, azimuth, antenna)
        frequencies =  a.GetFrequencies(t,p)
        E_theta, E_phi = a.GetElectricFieldFactorAnt(zenith,azimuth, frequencies,antenna, c=2.997e8,Z_0=119.9169*np.pi)        
        plt.scatter(frequencies,np.angle(E_theta)-np.angle(E_phi))
        plt.axhline(np.pi)
        plt.axhline(np.pi/2.)
        plt.axhline(np.pi/4.)
        plt.axhline(0)
#         plt.figure()
#         plt.plot(frequencies,np.abs(E_theta),label='E theta ')
#         plt.plot(frequencies,np.abs(E_phi),label='E phi')
#         plt.legend()


        plt.figure()
        plt.plot(frequencies,np.unwrap(np.angle(E_theta)),label='E theta ')
        plt.plot(frequencies,np.unwrap(np.angle(E_phi)),label='E phi')
        plt.legend()
        plt.show()
        
    if plot == 6:
        zenith = 180
        azimuth = 90
        antenna='DEW'        
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        t,p = a.GetThetaAndPhi(zenith, azimuth, antenna)
        frequencies =  a.GetFrequencies(t,p)
        pr = a.GetElectricFieldFactorAnt(zenith,azimuth, frequencies,antenna, c=2.997e8,Z_0=119.9169*np.pi)
        
        g = a.InterpolateGains(t,p,frequencies)
        
        A_eff = a.GetEffectiveHeight(g,frequencies)
        
        plt.figure()
        plt.plot(frequencies,np.abs(pr[0]),label='E_theta')
        plt.plot(frequencies,np.abs(pr[1]),label='E_phi')
        plt.plot(frequencies,A_eff,label='A_eff')
        plt.xlabel("Frequencies [MHz]")
        plt.ylabel("Antenna factor [m]")
        plt.legend()
        plt.title("Azimuth {} degres".format(azimuth))
        #plt.savefig("Debug_antenna_factor_{}.png".format(azimuth))
        plt.show()
        E = 1 + np.exp(1j*0.1)
        
        
    if plot ==5:
        zenith = 45
        azimuth = 0
        antenna='UPW'
        
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        
        t,p = a.GetThetaAndPhi(zenith, azimuth, antenna)
        frequencies =  a.GetFrequencies(t,p)
        a.GetElectricFieldFactorAnt(zenith,azimuth, frequencies,antenna, c=2.997e8,Z_0=119.9169*np.pi)
    
    if plot == 4:
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        
        zenith = 55
        azimuth = 100
        
        t, p = a.GetThetaAndPhiFlexible(zenith,azimuth,azi_orientation=90,zen_orientation=45)
        ta, pa = a.GetThetaAndPhiForUpwardNorth(zenith,azimuth)
#         ta, pa = a.GetThetaAndPhiForUpwardWest(zenith,azimuth)
        
    
    if plot == 3:
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        
        zenith = 30
        azimuth = 90
        f =  a.GetFrequencies(zenith,azimuth)
       # rt, it, rp, ip = a.GetVoltages(zenith,azimuth)
        
       # plt.figure()
       # plt.plot(f,np.abs(rt),label='rt')
       # plt.plot(f,np.abs(rp),label='rp')
       # plt.legend()
       # plt.show()
        
    
    if plot == 2:
        a = AntennaModel("WIPLD_antennamodel_firn.root")
        
        f =  a.GetFrequencies(20,30)
        zenith = 90
        azimuth = 90
        
        E_phi_1, E_theta_1 = a.GetElectricFieldFactorAnt(zenith,azimuth, f,"DEW")
        zenith = 90
        azimuth = 180
        
        E_phi_2, E_theta_2 = a.GetElectricFieldFactorAnt(zenith,azimuth, f,"DNS")
        
        
        plt.figure()
        plt.plot(f, np.abs(E_theta_1))
        plt.plot(f, np.abs(E_theta_2))
        
        plt.figure()
        plt.plot(f, np.abs(E_phi_1))
        plt.plot(f, np.abs(E_phi_2))
        plt.title("phi")
        
        plt.figure()
        plt.plot(f, np.unwrap(np.angle(E_theta_1)))
        plt.plot(f, np.unwrap(np.angle(E_theta_2)))
        
        plt.figure()
        plt.plot(f, np.unwrap(np.angle(E_phi_1)))
        plt.plot(f, np.unwrap(np.angle(E_phi_2)))
        
        
        plt.show()
        
        
    if plot == 1:
    

        
        def DrawWIPLD(ax,col='r'):
            ax.plot([0,0.5],[0,0],[0,0],c=col)
            ax.plot([0.2,0.2],[-0.3,0.3],[0,0],c=col)
            ax.plot([0.25,0.25],[-0.25,0.25],[0,0],c=col)
            ax.plot([0.3,0.3],[-0.2,0.2],[0,0],c=col)
            ax.plot([0.35,0.35],[-0.15,0.15],[0,0],c=col)
            ax.plot([0.4,0.4],[-0.1,0.1],[0,0],c=col)
        
        def DrawUpwardNorth(ax,col='r'):
            ax.plot([0,0],[0,np.cos(np.deg2rad(45))*0.5],[0,np.sin(np.deg2rad(45))*0.5],c=col)
            ax.plot([-0.3,0.3],[np.cos(np.deg2rad(45))*0.2,np.cos(np.deg2rad(45))*0.2], [np.sin(np.deg2rad(45))*0.2,np.sin(np.deg2rad(45))*0.2],c=col)
            ax.plot([-0.25,0.25],[np.cos(np.deg2rad(45))*0.25,np.cos(np.deg2rad(45))*0.25], 
                    [np.sin(np.deg2rad(45))*0.25,np.sin(np.deg2rad(45))*0.25],c=col)
            ax.plot([-0.2,0.2],[np.cos(np.deg2rad(45))*0.3,np.cos(np.deg2rad(45))*0.3],[np.sin(np.deg2rad(45))*0.3,np.sin(np.deg2rad(45))*0.3],c=col)
            ax.plot([-0.15,0.15],[np.cos(np.deg2rad(45))*0.35,np.cos(np.deg2rad(45))*0.35], 
                    [np.sin(np.deg2rad(45))*0.35,np.sin(np.deg2rad(45))*0.35],c=col)
            ax.plot([-0.1,0.1],[np.cos(np.deg2rad(45))*0.4,np.cos(np.deg2rad(45))*0.4],[np.sin(np.deg2rad(45))*0.4,np.sin(np.deg2rad(45))*0.4],c=col)
        
        def DrawUpwardWest(ax,col='r'):
            ax.plot([0,-1*np.cos(np.deg2rad(45))*0.5],[0,0],[0,np.sin(np.deg2rad(45))*0.5],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.2,-1.*np.cos(np.deg2rad(45))*0.2],[-0.3,0.3],
                    [np.sin(np.deg2rad(45))*0.2,np.sin(np.deg2rad(45))*0.2],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.2,-1.*np.cos(np.deg2rad(45))*0.2],[-0.3,0.3],  
                    [np.sin(np.deg2rad(45))*0.2,np.sin(np.deg2rad(45))*0.2],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.25,-1.*np.cos(np.deg2rad(45))*0.25], [-0.25,0.25],         
                    [np.sin(np.deg2rad(45))*0.25,np.sin(np.deg2rad(45))*0.25],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.3,-1.*np.cos(np.deg2rad(45))*0.3],[-0.2,0.2], 
                    [np.sin(np.deg2rad(45))*0.3,np.sin(np.deg2rad(45))*0.3],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.35,-1.*np.cos(np.deg2rad(45))*0.35],[-0.15,0.15],
                    [np.sin(np.deg2rad(45))*0.35,np.sin(np.deg2rad(45))*0.35],c=col)
            ax.plot([-1.*np.cos(np.deg2rad(45))*0.4,-1.*np.cos(np.deg2rad(45))*0.4],[-0.1,0.1], 
                    [np.sin(np.deg2rad(45))*0.4,np.sin(np.deg2rad(45))*0.4],c=col)
        
        def DrawDownNS(ax,col='c'):
            ax.plot([0,0],[0,0],[0,-0.5],c=col)
            ax.plot([0,0],[-0.3,0.3],[-0.2,-0.2],c=col)
            ax.plot([0,0],[-0.25,0.25],[-0.25,-0.25],c=col)
            ax.plot([0,0],[-0.2,0.2],[-0.3,-0.3],c=col)
            ax.plot([0,0],[-0.15,0.15],[-0.35,-0.35],c=col)
            ax.plot([0,0],[-0.1,0.1],[-0.4,-0.4],c=col)

        def DrawDownEW(ax,col='m'):
            ax.plot([0,0],[0,0],[0,-0.5],c=col)
            ax.plot([-0.3,0.3],[0,0],[-0.2,-0.2],c=col)
            ax.plot([-0.25,0.25],[0,0],[-0.25,-0.25],c=col)
            ax.plot([-0.2,0.2],[0,0],[-0.3,-0.3],c=col)
            ax.plot([-0.15,0.15],[0,0],[-0.35,-0.35],c=col)
            ax.plot([-0.1,0.1],[0,0],[-0.4,-0.4],c=col)

        def DrawVector(ax,zenith,azimuth,col='y'):
            zenith = np.deg2rad(zenith)
            azimuth = np.deg2rad(azimuth)
            ax.plot([0,np.sin(zenith)*np.cos(azimuth)],[0,np.sin(zenith)*np.sin(azimuth)],[0,np.cos(zenith)],c=col)
            
            
        a = AntennaModel("WIPLD_antennamodel_firn.root")
            
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot([-1,1],[0,0],[0,0],c='k')
        ax.plot([0,0],[-1,1],[0,0],c='k')
        ax.plot([0,0],[0,0],[-1,1],c='k')
        
        zenith = 45
        azimuth = 90
        DrawVector(ax,zenith,azimuth,col='g')
        t, p = a.GetThetaAndPhiForUpwardNorth(zenith,azimuth)
        
        
        DrawVector(ax,90-t,p,col='r')
        DrawWIPLD(ax,col='r')
        
#         DrawUpwardNorth(ax,col='g')
#         DrawUpwardWest(ax,col='b')
        
        DrawUpwardNorth(ax,col='g')
#         DrawDownNS(ax)
        plt.show()
    
    
    if plot == 0:
    
        names = {"WIPLD_antennamodel_air.root":"air", "WIPLD_antennamodel_firn.root":"firn"}
        f1, ax1 = plt.subplots()
        f2,ax2 = plt.subplots()
        f3,ax3 = plt.subplots()
        f4,ax4 =  plt.subplots()
    
        for filename in ["WIPLD_antennamodel_air.root"]:#,"WIPLD_antennamodel_firn.root"]:
    
            a = AntennaModel(f=filename)


            zenith  = 180
            azimuth = 90
            ta, pa = a.GetThetaAndPhiForDownwardEWaligned(zenith,azimuth)

            f =  a.GetFrequencies(ta,pa)

    
            re_t, im_t, re_p, im_p = a.GetVoltages(ta,pa)
    
            tv = re_t + 1j*im_t
            pv = re_p + 1j*im_p
    
            mag_t = np.absolute(tv)
            mag_p = np.absolute(pv)

            V_theta, V_phi = a.InterpolateVoltages(ta,pa,f)

            ax1.set_xlabel("Frequency [MHz]")
    

            ax2.plot(f,a.Re_Z,label="RE {}".format(names[filename]))
            ax2.plot(f,a.Im_Z,label="IM {}".format(names[filename]))
        
    
            E_theta, E_phi = a.GetElectricFieldFactor(ta,pa, f,c=2.997e8,Z_0=119.9169*np.pi)

            unwrapped_phase = np.unwrap(np.angle(E_phi))
            grp = []
            for i in xrange(V_theta.shape[0]-1):
                grp.append((unwrapped_phase[i+1]-unwrapped_phase[i])/((f[i+1]-f[i])*2*np.pi*1e6*-1))

            ax1.plot(f[:-1],grp)
        
        
            ax3.plot(f,np.abs(E_phi),label='E_phi, {}'.format(names[filename]))

            ax3.plot(f,np.abs(E_theta),label='E_theta, {}'.format(names[filename]))
            ax3.set_yscale('log')
            ax3.set_ylabel("Multiplicative factor [m]")
            ax3.set_xlabel("Frequency [MHz]")   
            plt.axvline(100,linestyle='--',c='k')
        
            ax4.plot(f,np.angle(E_phi),label='E_phi, {}'.format(names[filename]))
            ax4.plot(f,np.angle(E_theta),label='E_theta, {}'.format(names[filename]))
        
        ax2.legend()
        ax3.legend()
        
        plt.show()
