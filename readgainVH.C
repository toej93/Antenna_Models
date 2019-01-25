#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "TMath.h"

const Double_t spedci=2.99792458e8;// in m/s

Double_t* GaintoHeffV(Double_t freq, Double_t nz, Double_t ant_theta, Double_t ant_phi);
Double_t* GetGainV(Double_t freq, Double_t theta, Double_t phi);

Double_t* GaintoHeffH(Double_t freq, Double_t nz, Double_t ant_theta, Double_t ant_phi);
Double_t* GetGainH(Double_t freq, Double_t theta, Double_t phi);

const Int_t freq_step = 60;
const Int_t ang_step = 2664;
Double_t Vgain[freq_step][ang_step];
Double_t Vphase[freq_step][ang_step];
Double_t VFreq[freq_step];

Double_t Hgain[freq_step][ang_step];
Double_t Hphase[freq_step][ang_step];
Double_t HFreq[freq_step];

void readgainVH(){
  
  //string filenameV="/home/uzair/Documents/rasta/VPol_XFDTD_FeedPlate.txt";
  string filenameV="/home/uzair/Documents/raytrace/heff/VPol_XFDTD_Crossbar.txt";
  ifstream NecOutV( filenameV.c_str() );
    
  string line;
  if ( NecOutV.is_open() ) {
    while (NecOutV.good() ) {
      
      for (Int_t i=0; i<freq_step; i++) {
	getline (NecOutV, line);
	if ( line.substr(0, line.find_first_of(":")) == "freq ") {
	  VFreq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
	  //cout<<"freq["<<i<<"] = "<<VFreq[i]<<" MHz"<<endl;
	  getline (NecOutV, line); //read SWR
	  
	  getline (NecOutV, line); //read names
	  
	  for (Int_t j=0; j<ang_step; j++) {
	    getline (NecOutV, line); //read data line
	    //Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
	    ////for_feedplate
	    //Vgain[i][j] = atof( line.substr( 25, 34 ).c_str() ); 
	    //Vphase[i][j] = atof( line.substr( 35 ).c_str() ); 

	    ////for crossbar file
	     Vgain[i][j] = atof( line.substr( 35, 45 ).c_str() );
	     Vphase[i][j] = atof( line.substr( 45 ).c_str() );  
	     //cout<<Vgain[i][j]<<" "<<Vphase[i][j]<<" "<<Vgain[i][j]<<" "<<Vgain[i][j]<<endl;
	    
	  }// end ang_step
          
	}// end check freq label
        
      }// end freq_step
      
    }// end while NecOutV.good
    NecOutV.close();
  }// end if file open

  //for (Int_t i=200; i<801; i++) {
  //Double_t heff=GaintoHeffV(400,1.78, 75, 45);			
    //cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz, "<<" VGain : "<<Vgain[i][j]<<", VPhase : "<<Vphase[i][j]<<" Heff : "<<heff<<endl;
    //cout<<"freq = "<<400<<" MHz, "<<" Heff : "<<heff<<endl;
    //}
  
  string filenameH="/home/uzair/Documents/rasta/HPol_XFDTD_CurrentMod.txt";
  ifstream NecOutH( filenameH.c_str() );
    
  //string line;
  if ( NecOutH.is_open() ) {
    while (NecOutH.good() ) {
      
      for (Int_t i=0; i<freq_step; i++) {
  	getline (NecOutH, line);
  	if ( line.substr(0, line.find_first_of(":")) == "freq ") {
  	  HFreq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
  	  //cout<<"freq["<<i<<"] = "<<HFreq[i]<<" MHz"<<endl;
  	  getline (NecOutH, line); //read SWR
	  
  	  getline (NecOutH, line); //read names
	  
  	  for (Int_t j=0; j<ang_step; j++) {
  	    getline (NecOutH, line); //read data line
  	    //Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
  	    Hgain[i][j] = atof( line.substr( 25, 33 ).c_str() );  // read gain (not dB)
  	    Hphase[i][j] = atof( line.substr( 34 ).c_str() );  // read gain (not dB)
	    //cout<<Hgain[i][j]<<" "<<Hphase[i][j]<<endl;
	    
	    // Hgain[i][j] = atof( line.substr( 35, 45 ).c_str() );  // read gain (not dB)
  	    // Hphase[i][j] = atof( line.substr( 45 ).c_str() );  // read gain (not dB)
	  }// end ang_step
          
  	}// end check freq label
        
      }// end freq_step
      
    }// end while NecOutH.good
    NecOutH.close();
  }// end if file open
  
  
}
    
// }// end ReadVgain

Double_t* GaintoHeffV(Double_t freq, Double_t nz, Double_t ant_theta, Double_t ant_phi){
  
  Double_t *output=new Double_t[2];
  Double_t Zl=50;//Ohms-> Load resistance or radiation resistance
  Double_t Z0=120*TMath::Pi();//Ohms-> The impedance of free space
  Double_t *getres=GetGainV(freq,ant_theta,ant_phi);
    
  output[0]=2*sqrt((getres[0]/4*TMath::Pi())*pow(spedci/(nz*freq*pow(10,6)),2)*(Zl/(Z0)));
  output[1]=getres[1];
  delete getres;
  
  return output; 

}

Double_t* GetGainV(Double_t freq, Double_t theta, Double_t phi){

  Double_t *output=new Double_t[2];
  Double_t finalG=0;
  Double_t finalPH=0;
  if(fmod(theta,5.0)==0 && fmod(phi,5.0)==0){
    //cout<<"both are divisible 5"<<endl;
    
    Int_t phibin=phi/5.0;
    Int_t thetabin=theta/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
    
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(VFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=VFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=VFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=VFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    
    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;
    Int_t iibin=thetabin+phibin*37;
    
    Double_t Gup_freq=Vgain[abvfreqbin][iibin];    
    Double_t Glow_freq=Vgain[belfreqbin][iibin];

    Double_t PHup_freq=Vphase[abvfreqbin][iibin];    
    Double_t PHlow_freq=Vphase[belfreqbin][iibin];
    
    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }
  
  if(fmod(theta,5.0)==0 && fmod(phi,5.0)!=0){
    //cout<<"Only theta is divisible 5"<<endl;
    
    Double_t phidiff1=0;
    Double_t phidiff2=10000;
    Double_t abvphi=0;
    Double_t belphi=0;
    Int_t abvphibin=0;
    Int_t belphibin=0;
    
    //find closest phi
    for (Int_t i=0; i<72; i++) {
      phidiff1=fabs((Double_t)(i*5.0)-phi);
      if (phidiff1<phidiff2){
	phidiff2=phidiff1;
	abvphi=i*5.0;
	abvphibin=i;
      }
    }
    
    if(phi<abvphi){
      abvphi=abvphi;
      belphi=(abvphibin-1)*5.0;
      abvphibin=abvphibin;
      belphibin=abvphibin-1;
    }
    
    if(phi>abvphi){
      belphi=abvphi;
      abvphi=(abvphibin+1)*5.0;
      belphibin=abvphibin;
      abvphibin=abvphibin+1;
    }
    
    //cout<<"the phi bins are "<<phidiff2<<" "<<abvphi<<" "<<belphi<<" "<<abvphibin<<" "<<belphibin<<endl;
    
    Int_t thetabin=theta/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
  
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(VFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=VFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=VFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=VFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }

    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;

    Int_t iibin=thetabin+abvphibin*37;
    Int_t ijbin=thetabin+belphibin*37;
    
    Double_t Gup_freq=Vgain[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t Glow_freq=Vgain[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));

    Double_t PHup_freq=Vphase[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHlow_freq=Vphase[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
  
    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }

  if(fmod(theta,5.0)!=0 && fmod(phi,5.0)==0){
    //cout<<"Only phi is divisible 5"<<endl;
    
    Double_t thetadiff1=0;
    Double_t thetadiff2=10000;
    Double_t abvtheta=0;
    Double_t beltheta=0;
    Int_t abvthetabin=0;
    Int_t belthetabin=0;
    
    //find closest theta
    for (Int_t i=0; i<37; i++) {
      thetadiff1=fabs((Double_t)(i*5.0)-theta);
      if (thetadiff1<thetadiff2){
	thetadiff2=thetadiff1;
	abvtheta=i*5.0;
	abvthetabin=i;
      }
    }
    
    if(theta<abvtheta){
      abvtheta=abvtheta;
      beltheta=(abvthetabin-1)*5.0;
      abvthetabin=abvthetabin;
      belthetabin=abvthetabin-1;
    }
    
    if(theta>abvtheta){
      beltheta=abvtheta;
      abvtheta=(abvthetabin+1)*5.0;
      belthetabin=abvthetabin;
      abvthetabin=abvthetabin+1;
    }
    
    Int_t phibin=phi/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
  
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(VFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=VFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=VFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=VFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }

    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;
    
    Int_t iibin=abvthetabin+phibin*37;
    Int_t ijbin=belthetabin+phibin*37;
    
    Double_t Gup_freq=Vgain[abvfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Vgain[abvfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t Glow_freq=Vgain[belfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Vgain[belfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));

    Double_t PHup_freq=Vphase[abvfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Vphase[abvfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHlow_freq=Vphase[belfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Vphase[belfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
  
    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }

  if(fmod(theta,5.0)!=0 && fmod(phi,5.0)!=0){
    //cout<<"Both are not divisible 5"<<endl;
    
    Double_t phidiff1=0;
    Double_t phidiff2=10000;
    Double_t abvphi=0;
    Double_t belphi=0;
    Int_t abvphibin=0;
    Int_t belphibin=0;  
    
    //find closest phi
    for (Int_t i=0; i<72; i++) {
      phidiff1=fabs((Double_t)(i*5.0)-phi);
      if (phidiff1<phidiff2){
	phidiff2=phidiff1;
	abvphi=i*5.0;
	abvphibin=i;
      }
    }
    
    if(phi<abvphi){
      abvphi=abvphi;
      belphi=(abvphibin-1)*5.0;
      abvphibin=abvphibin;
      belphibin=abvphibin-1;
    }
    
    if(phi>abvphi){
      belphi=abvphi;
      abvphi=(abvphibin+1)*5.0;
      belphibin=abvphibin;
      abvphibin=abvphibin+1;
    }
    
    //cout<<"the phi bins are "<<phidiff2<<" "<<abvphi<<" "<<belphi<<" "<<abvphibin<<" "<<belphibin<<endl;
    
    Double_t thetadiff1=0;
    Double_t thetadiff2=10000;
    Double_t abvtheta=0;
    Double_t beltheta=0;
    Int_t abvthetabin=0;
    Int_t belthetabin=0;
  
    //find closest theta
    for (Int_t i=0; i<37; i++) {
      thetadiff1=fabs((Double_t)(i*5.0)-theta);
      if (thetadiff1<thetadiff2){
	thetadiff2=thetadiff1;
	abvtheta=i*5.0;
	abvthetabin=i;
      }
    }
    
    if(theta<abvtheta){
      abvtheta=abvtheta;
      beltheta=(abvthetabin-1)*5.0;
      abvthetabin=abvthetabin;
      belthetabin=abvthetabin-1;
    }
    
    if(theta>abvtheta){
      beltheta=abvtheta;
      abvtheta=(abvthetabin+1)*5.0;
      belthetabin=abvthetabin;
      abvthetabin=abvthetabin+1;
    }
    
    //cout<<thetadiff2<<" "<<abvtheta<<" "<<beltheta<<" "<<abvthetabin<<" "<<belthetabin<<endl;
  
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
    
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(VFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=VFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=VFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=VFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    
    //upfreq
    //cout<<" For up freq "<<endl;
    
    //uptheta
    Int_t iibin=abvthetabin+abvphibin*37;
    Int_t ijbin=abvthetabin+belphibin*37;
    
    Double_t Gup_theta=Vgain[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHup_theta=Vphase[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    //low theta
    Int_t jibin=belthetabin+abvphibin*37;
    Int_t jjbin=belthetabin+belphibin*37;
    
    Double_t Glow_theta=Vgain[abvfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[abvfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHlow_theta=Vphase[abvfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[abvfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    Double_t Gup_freq=Gup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Glow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHup_freq=PHup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+PHlow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    
    //downfreq
    //cout<<" For low freq "<<endl;

    //uptheta
    Gup_theta=Vgain[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    PHup_theta=Vphase[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
  
    //low theta
    Glow_theta=Vgain[belfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vgain[belfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    PHlow_theta=Vphase[belfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Vphase[belfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    //cout<<"Gup and Glow thetas "<<Gup_theta<<" "<<Glow_theta<<endl;
    
    Double_t Glow_freq=Gup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Glow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHlow_freq=PHup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+PHlow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));

    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }
  //cout<<"final G "<<finalG<<endl;

  output[0]=finalG;
  output[1]=finalPH;
  
  return output;
  
}

Double_t* GaintoHeffH(Double_t freq, Double_t nz, Double_t ant_theta, Double_t ant_phi){

  Double_t *output=new Double_t[2];
  Double_t Zl=50;//Ohms-> Load resistance or radiation resistance
  Double_t Z0=120*TMath::Pi();//Ohms-> The impedance of free space
  Double_t *getres=GetGainH(freq,ant_theta,ant_phi);
    
  output[0]=2*sqrt((getres[0]/4*TMath::Pi())*pow(spedci/(nz*freq*pow(10,6)),2)*(Zl/(Z0)));
  output[1]=getres[1];
  delete getres;
  
  return output; 
}

Double_t* GetGainH(Double_t freq, Double_t theta, Double_t phi){
  Double_t *output=new Double_t[2];
  Double_t finalG=0;
  Double_t finalPH=0;
  
  if(fmod(theta,5.0)==0 && fmod(phi,5.0)==0){
    //cout<<"both are divisible 5"<<endl;
    
    Int_t phibin=phi/5.0;
    Int_t thetabin=theta/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
    
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(HFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=HFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=HFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=HFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    
    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;
    Int_t iibin=thetabin+phibin*37;
    
    Double_t Gup_freq=Hgain[abvfreqbin][iibin];
    Double_t Glow_freq=Hgain[belfreqbin][iibin];

    Double_t PHup_freq=Hphase[abvfreqbin][iibin];
    Double_t PHlow_freq=Hphase[belfreqbin][iibin];
   
    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));

  }
  
  if(fmod(theta,5.0)==0 && fmod(phi,5.0)!=0){
    //cout<<"Only theta is divisible 5"<<endl;
    
    Double_t phidiff1=0;
    Double_t phidiff2=10000;
    Double_t abvphi=0;
    Double_t belphi=0;
    Int_t abvphibin=0;
    Int_t belphibin=0;
    
    //find closest phi
    for (Int_t i=0; i<72; i++) {
      phidiff1=fabs((Double_t)(i*5.0)-phi);
      if (phidiff1<phidiff2){
	phidiff2=phidiff1;
	abvphi=i*5.0;
	abvphibin=i;
      }
    }
    
    if(phi<abvphi){
      abvphi=abvphi;
      belphi=(abvphibin-1)*5.0;
      abvphibin=abvphibin;
      belphibin=abvphibin-1;
    }
    
    if(phi>abvphi){
      belphi=abvphi;
      abvphi=(abvphibin+1)*5.0;
      belphibin=abvphibin;
      abvphibin=abvphibin+1;
    }
    
    //cout<<"the phi bins are "<<phidiff2<<" "<<abvphi<<" "<<belphi<<" "<<abvphibin<<" "<<belphibin<<endl;
    
    Int_t thetabin=theta/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
  
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(HFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=HFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=HFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=HFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }

    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;

    Int_t iibin=thetabin+abvphibin*37;
    Int_t ijbin=thetabin+belphibin*37;
    
    Double_t Gup_freq=Hgain[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t Glow_freq=Hgain[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));  

    Double_t PHup_freq=Hphase[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHlow_freq=Hphase[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));  

    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }

  if(fmod(theta,5.0)!=0 && fmod(phi,5.0)==0){
    //cout<<"Only phi is divisible 5"<<endl;
    
    Double_t thetadiff1=0;
    Double_t thetadiff2=10000;
    Double_t abvtheta=0;
    Double_t beltheta=0;
    Int_t abvthetabin=0;
    Int_t belthetabin=0;
    
    //find closest theta
    for (Int_t i=0; i<37; i++) {
      thetadiff1=fabs((Double_t)(i*5.0)-theta);
      if (thetadiff1<thetadiff2){
	thetadiff2=thetadiff1;
	abvtheta=i*5.0;
	abvthetabin=i;
      }
    }
    
    if(theta<abvtheta){
      abvtheta=abvtheta;
      beltheta=(abvthetabin-1)*5.0;
      abvthetabin=abvthetabin;
      belthetabin=abvthetabin-1;
    }
    
    if(theta>abvtheta){
      beltheta=abvtheta;
      abvtheta=(abvthetabin+1)*5.0;
      belthetabin=abvthetabin;
      abvthetabin=abvthetabin+1;
    }
    
    Int_t phibin=phi/5.0;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
  
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(HFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=HFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=HFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=HFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }

    //cout<<"abvfreq "<<abvfreq<<" belfreq "<<belfreq<<endl;
    
    Int_t iibin=abvthetabin+phibin*37;
    Int_t ijbin=belthetabin+phibin*37;
    
    Double_t Gup_freq=Hgain[abvfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Hgain[abvfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t Glow_freq=Hgain[belfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Hgain[belfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));

    Double_t PHup_freq=Hphase[abvfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Hphase[abvfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHlow_freq=Hphase[belfreqbin][iibin]*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Hphase[belfreqbin][ijbin]*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
  
    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }

  if(fmod(theta,5.0)!=0 && fmod(phi,5.0)!=0){
    //cout<<"Both are not divisible 5"<<endl;
    
    Double_t phidiff1=0;
    Double_t phidiff2=10000;
    Double_t abvphi=0;
    Double_t belphi=0;
    Int_t abvphibin=0;
    Int_t belphibin=0;  
    
    //find closest phi
    for (Int_t i=0; i<72; i++) {
      phidiff1=fabs((Double_t)(i*5.0)-phi);
      if (phidiff1<phidiff2){
	phidiff2=phidiff1;
	abvphi=i*5.0;
	abvphibin=i;
      }
    }
    
    if(phi<abvphi){
      abvphi=abvphi;
      belphi=(abvphibin-1)*5.0;
      abvphibin=abvphibin;
      belphibin=abvphibin-1;
    }
    
    if(phi>abvphi){
      belphi=abvphi;
      abvphi=(abvphibin+1)*5.0;
      belphibin=abvphibin;
      abvphibin=abvphibin+1;
    }
    
    //cout<<"the phi bins are "<<phidiff2<<" "<<abvphi<<" "<<belphi<<" "<<abvphibin<<" "<<belphibin<<endl;
    
    Double_t thetadiff1=0;
    Double_t thetadiff2=10000;
    Double_t abvtheta=0;
    Double_t beltheta=0;
    Int_t abvthetabin=0;
    Int_t belthetabin=0;
  
    //find closest theta
    for (Int_t i=0; i<37; i++) {
      thetadiff1=fabs((Double_t)(i*5.0)-theta);
      if (thetadiff1<thetadiff2){
	thetadiff2=thetadiff1;
	abvtheta=i*5.0;
	abvthetabin=i;
      }
    }
    
    if(theta<abvtheta){
      abvtheta=abvtheta;
      beltheta=(abvthetabin-1)*5.0;
      abvthetabin=abvthetabin;
      belthetabin=abvthetabin-1;
    }
    
    if(theta>abvtheta){
      beltheta=abvtheta;
      abvtheta=(abvthetabin+1)*5.0;
      belthetabin=abvthetabin;
      abvthetabin=abvthetabin+1;
    }
    
    //cout<<thetadiff2<<" "<<abvtheta<<" "<<beltheta<<" "<<abvthetabin<<" "<<belthetabin<<endl;
  
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
    
    //find closest freq
    for (Int_t i=0; i<freq_step; i++) {
      freqdiff1=fabs(HFreq[i]-freq);
      if (freqdiff1<freqdiff2){
	freqdiff2=freqdiff1;
	abvfreq=HFreq[i];
	abvfreqbin=i;
      }
    }
    
    if(freq<abvfreq){
      abvfreq=abvfreq;
      belfreq=HFreq[abvfreqbin-1];
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(freq>abvfreq){
      belfreq=abvfreq;
      abvfreq=HFreq[abvfreqbin+1];
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    
    //upfreq
    //cout<<" For up freq "<<endl;
    
    //uptheta
    Int_t iibin=abvthetabin+abvphibin*37;
    Int_t ijbin=abvthetabin+belphibin*37;
    
    Double_t Gup_theta=Hgain[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHup_theta=Hphase[abvfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[abvfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    //low theta
    Int_t jibin=belthetabin+abvphibin*37;
    Int_t jjbin=belthetabin+belphibin*37;
    
    Double_t Glow_theta=Hgain[abvfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[abvfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    Double_t PHlow_theta=Hphase[abvfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[abvfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    //cout<<"Gup and Glow thetas "<<Gup_theta<<" "<<Glow_theta<<endl;
    
    Double_t Gup_freq=Gup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Glow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHup_freq=PHup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+PHlow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    
    //downfreq
    //cout<<" For low freq "<<endl;

    //uptheta
    Gup_theta=Hgain[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    PHup_theta=Hphase[belfreqbin][iibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[belfreqbin][ijbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
  
    //low theta
    Glow_theta=Hgain[belfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hgain[belfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    PHlow_theta=Hphase[belfreqbin][jibin]*(1-fabs(phi-abvphi)/(abvphi-belphi))+Hphase[belfreqbin][jjbin]*(1-fabs(phi-belphi)/(abvphi-belphi));
    
    //cout<<"Gup and Glow thetas "<<Gup_theta<<" "<<Glow_theta<<endl;
    
    Double_t Glow_freq=Gup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+Glow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));
    Double_t PHlow_freq=PHup_theta*(1-fabs(theta-abvtheta)/(abvtheta-beltheta))+PHlow_theta*(1-fabs(theta-beltheta)/(abvtheta-beltheta));

    //cout<<"Gup and Glow freqs "<<Gup_freq<<" "<<Glow_freq<<endl;
    
    finalG=Gup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+Glow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
    finalPH=PHup_freq*(1-fabs(freq-abvfreq)/(abvfreq-belfreq))+PHlow_freq*(1-fabs(freq-belfreq)/(abvfreq-belfreq));
  }
  //cout<<"final G "<<finalG<<endl;

  output[0]=finalG;
  output[1]=finalPH;
  
  return output;
  
}

