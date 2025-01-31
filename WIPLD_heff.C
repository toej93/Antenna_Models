#include "TComplex.h"
#include "TFile.h"
#include "TTree.h"
Double_t *WIPLD_heff(Float_t mytheta,Float_t myphi,Float_t myfreq,TFile *ariannaAntennaFile, TTree *AntTree, TTree *ZTree){

  Double_t * output=new Double_t[10];
  for(Int_t io=0;io<10;io++){
    output[io]=0;
  }

  // Float_t mytheta =54.0;
  // Float_t myphi = 64.0;
  // Float_t myfreq = 257.6;

  Int_t N = 200;
  const Int_t N2 = 240;

  // Double_t thetastp=2.0;
  // Double_t phistp=8.0;  
  // TFile *ariannaAntennaFile = new TFile("WIPLD_ARA_bicone_v1.root");

  Double_t thetastp=1.0;
  Double_t phistp=1.0;  
  //TFile *ariannaAntennaFile = new TFile("WIPLD_antennamodel_firn.root");

  // TTree *AntTree = (TTree*)ariannaAntennaFile->Get("AntTree");
  // TTree *ZTree = (TTree*)ariannaAntennaFile->Get("ZTree");
  
  Float_t freq[N2];
  Double_t gains[N2];
  Double_t v_re_phi[N2];
  Double_t v_im_phi[N2];
  Double_t v_re_theta[N2];
  Double_t v_im_theta[N2];
  Float_t thetas;
  Float_t phis;

  Double_t z_re[N2];
  Double_t z_im[N2];
  
  TComplex Z[N2];
  
  Float_t Z_mag[N2];
  Float_t V_theta_mag[N2];
  Float_t V_theta_phase[N2];
  Float_t V_phi_mag[N2];
  Float_t V_phi_phase[N2];
   
  Float_t pi = TMath::Pi();
  
  AntTree->SetBranchAddress("frequencies", &freq);
  AntTree->SetBranchAddress("gains", &gains);
  AntTree->SetBranchAddress("Re_phi", &v_re_phi);
  AntTree->SetBranchAddress("Im_phi", &v_im_phi);
  AntTree->SetBranchAddress("Re_theta", &v_re_theta);
  AntTree->SetBranchAddress("Im_theta", &v_im_theta);
  AntTree->SetBranchAddress("thetas", &thetas);
  AntTree->SetBranchAddress("phis", &phis);
  AntTree->BuildIndex("thetas","phis");
  
  // fill the antenna impedance. 
  // This tree only has one entry.
  ZTree->SetBranchAddress("Im_Z", &z_im);
  ZTree->SetBranchAddress("Re_Z", &z_re);
  ZTree->GetEntry(0);
  
  Float_t h_phi[2][4];
  Float_t h_theta[2][4];
  TComplex h_thetaC[2];
  TComplex h_phiC[2];	
  TComplex V_theta[2];
  TComplex V_phi[2];

  Double_t G_theta[2];
  Double_t G_phi[2];

  Double_t G_thetaF[2];
  Double_t G_phiF[2];

  Double_t PHG_thetaF[2];
  Double_t PHG_phiF[2];

  Double_t THG_thetaF[2];
  Double_t THG_phiF[2];

  Double_t TPG_thetaF;
  Double_t TPG_phiF;

  Double_t z0 = 377.;
  Double_t v0 = 1.0;

  Double_t h_phiF[2][4];
  Double_t h_thetaF[2][4];

  Double_t PHh_phiF[2][4];
  Double_t PHh_thetaF[2][4];

  Double_t THh_phiF[4];
  Double_t THh_thetaF[4];

  Double_t TPh_phiF[4];
  Double_t TPh_thetaF[4];

  Double_t wavelength;
  TComplex norm;
  Int_t entry;
  Double_t A_eff_theta;
  Double_t A_eff_phi;
  ///////////////////////////////////////////////////////////////////////
  if(fmod(mytheta,thetastp)==0 && fmod(myphi,phistp)==0){
    //cout<<"both are divisible"<<endl;

    entry = AntTree->GetEntryNumberWithIndex(mytheta,myphi);
    AntTree->GetEntry(entry);
    //cout<<"entry is "<<entry<<endl;

    Int_t phibin=myphi/phistp;
    Int_t thetabin=mytheta/thetastp;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;
    
    //find closest freq
    for (Int_t i=0; i<N; i++) {
      freqdiff1=fabs(freq[i]*pow(10,3)-myfreq);
      if (freqdiff1<freqdiff2){
  	freqdiff2=freqdiff1;
  	abvfreq=freq[i]*pow(10,3);
  	abvfreqbin=i;
      }
    }
    
    if(myfreq<abvfreq){
      abvfreq=abvfreq;
      belfreq=freq[abvfreqbin-1]*pow(10,3);
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(myfreq>abvfreq){
      belfreq=abvfreq;
      abvfreq=freq[abvfreqbin+1]*pow(10,3);
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    //cout<<abvfreq<<" "<<belfreq<<" "<<abvfreqbin<<" "<<belfreqbin<<endl; 

    //calculate heff for abv freq.
    V_theta[0] = TComplex(v_re_theta[abvfreqbin], v_im_theta[abvfreqbin]);
    V_phi[0] = TComplex(v_re_phi[abvfreqbin], v_im_phi[abvfreqbin]);
  
    wavelength = 0.299792/((Double_t) freq[abvfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
    norm = TComplex(0., -2 * wavelength/z0/v0);

    Z[abvfreqbin] = TComplex(z_re[abvfreqbin], z_im[abvfreqbin]);
    h_thetaC[0] = norm * Z[abvfreqbin] * V_theta[0];
    h_phiC[0]   = norm * Z[abvfreqbin] * V_phi[0];

    A_eff_theta = pow(TComplex::Abs(h_thetaC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
    A_eff_phi   = pow(TComplex::Abs(h_phiC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
    G_theta[0] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
    G_phi[0] = 4. * pi * A_eff_phi / (wavelength * wavelength);
 
    //calculate heff for low freq.
    V_theta[1] = TComplex(v_re_theta[belfreqbin], v_im_theta[belfreqbin]);
    V_phi[1] = TComplex(v_re_phi[belfreqbin], v_im_phi[belfreqbin]);
  
    wavelength = 0.299792/((Double_t) freq[belfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
    norm = TComplex(0., -2 * wavelength/z0/v0);

    Z[belfreqbin] = TComplex(z_re[belfreqbin], z_im[belfreqbin]);
    h_thetaC[1] = norm * Z[belfreqbin] * V_theta[1];
    h_phiC[1]   = norm * Z[belfreqbin] * V_phi[1];

    A_eff_theta = pow(TComplex::Abs(h_thetaC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
    A_eff_phi   = pow(TComplex::Abs(h_phiC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
    G_theta[1] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
    G_phi[1] = 4. * pi * A_eff_phi / (wavelength * wavelength);
    
    //calculate heff individual components for up and low
    for(Int_t iud=0;iud<2;iud++){
      h_theta[iud][0] = (Double_t) h_thetaC[iud].Re();
      h_theta[iud][1] = (Double_t) h_thetaC[iud].Im();
      h_theta[iud][2] = (Double_t) TComplex::Abs(h_thetaC[iud]);
      h_theta[iud][3] = (Double_t) h_thetaC[iud].Theta();//180./TMath::Pi();
     
      h_phi[iud][0] = (Double_t) h_phiC[iud].Re();
      h_phi[iud][1] = (Double_t) h_phiC[iud].Im();
      h_phi[iud][2] = (Double_t) TComplex::Abs(h_phiC[iud]);
      h_phi[iud][3] = (Double_t) h_phiC[iud].Theta();//180./TMath::Pi();

      // cout<<" "<<h_theta[iud][0]<<" "<<h_theta[iud][1]<<" "<<h_theta[iud][2]<<" "<<h_theta[iud][3]<<endl;
      // cout<<" "<<h_phi[iud][0]<<" "<<h_phi[iud][1]<<" "<<h_phi[iud][2]<<" "<<h_phi[iud][3]<<endl;
    }

    //calculate the weighted average
    for(Int_t ig=0;ig<4;ig++){
      h_thetaF[0][ig]=h_theta[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_theta[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
      h_phiF[0][ig]=h_phi[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_phi[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
      if(ig<4){
	output[ig]=h_thetaF[0][ig];
      }
      if(ig>3){
	output[ig]=h_phiF[0][ig-4];
      }
    }

    G_thetaF[0]=G_theta[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_theta[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
    G_phiF[0]=G_phi[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phi[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

    output[8]=G_thetaF[0];
    output[9]=G_phiF[0];
    //cout<<" "<<h_thetaF[0][0]<<" "<<h_thetaF[0][1]<<" "<<h_thetaF[0][2]<<" "<<h_thetaF[0][3]<<endl;
    //cout<<" "<<h_phiF[0][0]<<" "<<h_phiF[0][1]<<" "<<h_phiF[0][2]<<" "<<h_phiF[0][3]<<endl;

    
   }

  ///////////////////////////////////////////////////////////////////////
  if(fmod(mytheta,thetastp)==0 && fmod(myphi,phistp)!=0){
    //cout<<"only theta is divisible"<<endl;

    //cout<<(myphi-fmod(myphi,phistp))/phistp<<" "<<myfreq<<endl;

    Float_t dou_phi[2];
    dou_phi[0]=(myphi-fmod(myphi,phistp))+phistp;
    dou_phi[1]=(myphi-fmod(myphi,phistp));

    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;

    AntTree->GetEntry(0);
    //find closest freq
    for (Int_t i=0; i<N; i++) {
      freqdiff1=fabs(freq[i]*pow(10,3)-myfreq);
      if (freqdiff1<freqdiff2){
  	freqdiff2=freqdiff1;
  	abvfreq=freq[i]*pow(10,3);
  	abvfreqbin=i;
      }
    }
    
    if(myfreq<abvfreq){
      abvfreq=abvfreq;
      belfreq=freq[abvfreqbin-1]*pow(10,3);
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(myfreq>abvfreq){
      belfreq=abvfreq;
      abvfreq=freq[abvfreqbin+1]*pow(10,3);
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    //cout<<abvfreq<<" "<<belfreq<<" "<<abvfreqbin<<" "<<belfreqbin<<endl; 
    
    for(Int_t ipt=0;ipt<2;ipt++){
      entry = AntTree->GetEntryNumberWithIndex(mytheta,dou_phi[ipt]);
      AntTree->GetEntry(entry);
      //cout<<"entry is "<<entry<<endl;

      //calculate heff for abv freq.
      V_theta[0] = TComplex(v_re_theta[abvfreqbin], v_im_theta[abvfreqbin]);
      V_phi[0] = TComplex(v_re_phi[abvfreqbin], v_im_phi[abvfreqbin]);
  
      wavelength = 0.299792/((Double_t) freq[abvfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
      norm = TComplex(0., -2 * wavelength/z0/v0);

      Z[abvfreqbin] = TComplex(z_re[abvfreqbin], z_im[abvfreqbin]);
      h_thetaC[0] = norm * Z[abvfreqbin] * V_theta[0];
      h_phiC[0]   = norm * Z[abvfreqbin] * V_phi[0];

      A_eff_theta = pow(TComplex::Abs(h_thetaC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
      A_eff_phi   = pow(TComplex::Abs(h_phiC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
      G_theta[0] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
      G_phi[0] = 4. * pi * A_eff_phi / (wavelength * wavelength);
      
      //calculate heff for low freq.
      V_theta[1] = TComplex(v_re_theta[belfreqbin], v_im_theta[belfreqbin]);
      V_phi[1] = TComplex(v_re_phi[belfreqbin], v_im_phi[belfreqbin]);
  
      wavelength = 0.299792/((Double_t) freq[belfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
      norm = TComplex(0., -2 * wavelength/z0/v0);

      Z[belfreqbin] = TComplex(z_re[belfreqbin], z_im[belfreqbin]);
      h_thetaC[1] = norm * Z[belfreqbin] * V_theta[1];
      h_phiC[1]   = norm * Z[belfreqbin] * V_phi[1];

      A_eff_theta = pow(TComplex::Abs(h_thetaC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
      A_eff_phi   = pow(TComplex::Abs(h_phiC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
      G_theta[1] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
      G_phi[1] = 4. * pi * A_eff_phi / (wavelength * wavelength);
      
      //calculate heff individual components for up and low
      for(Int_t iud=0;iud<2;iud++){
	h_theta[iud][0] = (Double_t) h_thetaC[iud].Re();
	h_theta[iud][1] = (Double_t) h_thetaC[iud].Im();
	h_theta[iud][2] = (Double_t) TComplex::Abs(h_thetaC[iud]);
	h_theta[iud][3] = (Double_t) h_thetaC[iud].Theta();//180./TMath::Pi();
	
	h_phi[iud][0] = (Double_t) h_phiC[iud].Re();
	h_phi[iud][1] = (Double_t) h_phiC[iud].Im();
	h_phi[iud][2] = (Double_t) TComplex::Abs(h_phiC[iud]);
	h_phi[iud][3] = (Double_t) h_phiC[iud].Theta();//180./TMath::Pi();

	// cout<<" "<<h_theta[iud][0]<<" "<<h_theta[iud][1]<<" "<<h_theta[iud][2]<<" "<<h_theta[iud][3]<<endl;
	// cout<<" "<<h_phi[iud][0]<<" "<<h_phi[iud][1]<<" "<<h_phi[iud][2]<<" "<<h_phi[iud][3]<<endl;
      }

      //calculate the weighted average
      for(Int_t ig=0;ig<4;ig++){
	h_thetaF[ipt][ig]=h_theta[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_theta[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
	h_phiF[ipt][ig]=h_phi[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_phi[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
      }

      G_thetaF[ipt]=G_theta[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_theta[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
      G_phiF[ipt]=G_phi[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phi[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

      //cout<<" "<<h_thetaF[ipt][0]<<" "<<h_thetaF[ipt][1]<<" "<<h_thetaF[ipt][2]<<" "<<h_thetaF[ipt][3]<<endl;
      //cout<<" "<<h_phiF[ipt][0]<<" "<<h_phiF[ipt][1]<<" "<<h_phiF[ipt][2]<<" "<<h_phiF[ipt][3]<<endl;
      //cout<<" "<<endl;
    }

    for(Int_t ig=0;ig<4;ig++){
      PHh_thetaF[0][ig]=h_thetaF[0][ig]*(1-fabs(myphi-dou_phi[0])/(dou_phi[0]-dou_phi[1]))+h_thetaF[1][ig]*(1-fabs(myphi-dou_phi[1])/(dou_phi[0]-dou_phi[1]));
      PHh_phiF[0][ig]=h_phiF[0][ig]*(1-fabs(myphi-dou_phi[0])/(dou_phi[0]-dou_phi[1]))+h_phiF[1][ig]*(1-fabs(myphi-dou_phi[1])/(dou_phi[0]-dou_phi[1]));
      if(ig<4){
	output[ig]=PHh_thetaF[0][ig];
      }
      if(ig>3){
	output[ig]=PHh_phiF[0][ig-4];
      }
    }

     PHG_thetaF[0]=G_thetaF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_thetaF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
     PHG_phiF[0]=G_phiF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phiF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

     output[8]=PHG_thetaF[0];
     output[9]=PHG_phiF[0];

    //cout<<" "<<PHh_thetaF[0][0]<<" "<<PHh_thetaF[0][1]<<" "<<PHh_thetaF[0][2]<<" "<<PHh_thetaF[0][3]<<endl;
    //cout<<" "<<PHh_phiF[0][0]<<" "<<PHh_phiF[0][1]<<" "<<PHh_phiF[0][2]<<" "<<PHh_phiF[0][3]<<endl;
    
  }

  ////////////////////////////////////////////////////////////////////////////////////
  if(fmod(mytheta,thetastp)!=0 && fmod(myphi,phistp)==0){
    //cout<<"only phi is divisible"<<endl;
    
    //cout<<(mytheta-fmod(mytheta,thetastp))/thetastp<<endl;
    
    Float_t dou_theta[2];
    dou_theta[0]=(mytheta-fmod(mytheta,thetastp))+thetastp;
    dou_theta[1]=(mytheta-fmod(mytheta,thetastp));

    if(mytheta<0){
      dou_theta[1]=(mytheta-fmod(mytheta,thetastp))-thetastp;
      dou_theta[0]=(mytheta-fmod(mytheta,thetastp));
    }
    //cout<<dou_theta[0]<<" "<<dou_theta[1]<<endl;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;

    AntTree->GetEntry(0);
    //find closest freq
    for (Int_t i=0; i<N; i++) {
      freqdiff1=fabs(freq[i]*pow(10,3)-myfreq);
      if (freqdiff1<freqdiff2){
  	freqdiff2=freqdiff1;
  	abvfreq=freq[i]*pow(10,3);
  	abvfreqbin=i;
      }
    }
    
    if(myfreq<abvfreq){
      abvfreq=abvfreq;
      belfreq=freq[abvfreqbin-1]*pow(10,3);
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(myfreq>abvfreq){
      belfreq=abvfreq;
      abvfreq=freq[abvfreqbin+1]*pow(10,3);
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    //cout<<abvfreq<<" "<<belfreq<<" "<<abvfreqbin<<" "<<belfreqbin<<endl; 

    for(Int_t ipt=0;ipt<2;ipt++){
      entry = AntTree->GetEntryNumberWithIndex(dou_theta[ipt],myphi);
      AntTree->GetEntry(entry);
      //cout<<"entry is "<<entry<<endl;

      //calculate heff for abv freq.
      V_theta[0] = TComplex(v_re_theta[abvfreqbin], v_im_theta[abvfreqbin]);
      V_phi[0] = TComplex(v_re_phi[abvfreqbin], v_im_phi[abvfreqbin]);
  
      wavelength = 0.299792/((Double_t) freq[abvfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
      norm = TComplex(0., -2 * wavelength/z0/v0);

      Z[abvfreqbin] = TComplex(z_re[abvfreqbin], z_im[abvfreqbin]);
      h_thetaC[0] = norm * Z[abvfreqbin] * V_theta[0];
      h_phiC[0]   = norm * Z[abvfreqbin] * V_phi[0];

      A_eff_theta = pow(TComplex::Abs(h_thetaC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
      A_eff_phi   = pow(TComplex::Abs(h_phiC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
      G_theta[0] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
      G_phi[0] = 4. * pi * A_eff_phi / (wavelength * wavelength);
      
      //calculate heff for low freq.
      V_theta[1] = TComplex(v_re_theta[belfreqbin], v_im_theta[belfreqbin]);
      V_phi[1] = TComplex(v_re_phi[belfreqbin], v_im_phi[belfreqbin]);
  
      wavelength = 0.299792/((Double_t) freq[belfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
      norm = TComplex(0., -2 * wavelength/z0/v0);

      Z[belfreqbin] = TComplex(z_re[belfreqbin], z_im[belfreqbin]);
      h_thetaC[1] = norm * Z[belfreqbin] * V_theta[1];
      h_phiC[1]   = norm * Z[belfreqbin] * V_phi[1];

      A_eff_theta = pow(TComplex::Abs(h_thetaC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
      A_eff_phi   = pow(TComplex::Abs(h_phiC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
      G_theta[1] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
      G_phi[1] = 4. * pi * A_eff_phi / (wavelength * wavelength);
      
      //calculate heff individual components for up and low
      for(Int_t iud=0;iud<2;iud++){
	h_theta[iud][0] = (Double_t) h_thetaC[iud].Re();
	h_theta[iud][1] = (Double_t) h_thetaC[iud].Im();
	h_theta[iud][2] = (Double_t) TComplex::Abs(h_thetaC[iud]);
	h_theta[iud][3] = (Double_t) h_thetaC[iud].Theta();//180./TMath::Pi();
	
	h_phi[iud][0] = (Double_t) h_phiC[iud].Re();
	h_phi[iud][1] = (Double_t) h_phiC[iud].Im();
	h_phi[iud][2] = (Double_t) TComplex::Abs(h_phiC[iud]);
	h_phi[iud][3] = (Double_t) h_phiC[iud].Theta();//180./TMath::Pi();

	// cout<<" "<<h_theta[iud][0]<<" "<<h_theta[iud][1]<<" "<<h_theta[iud][2]<<" "<<h_theta[iud][3]<<endl;
	// cout<<" "<<h_phi[iud][0]<<" "<<h_phi[iud][1]<<" "<<h_phi[iud][2]<<" "<<h_phi[iud][3]<<endl;
      }

      //calculate the weighted average
      for(Int_t ig=0;ig<4;ig++){
	h_thetaF[ipt][ig]=h_theta[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_theta[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
	h_phiF[ipt][ig]=h_phi[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_phi[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
      }

      G_thetaF[ipt]=G_theta[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_theta[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
      G_phiF[ipt]=G_phi[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phi[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
      
      //cout<<" "<<h_thetaF[ipt][0]<<" "<<h_thetaF[ipt][1]<<" "<<h_thetaF[ipt][2]<<" "<<h_thetaF[ipt][3]<<endl;
      //cout<<" "<<h_phiF[ipt][0]<<" "<<h_phiF[ipt][1]<<" "<<h_phiF[ipt][2]<<" "<<h_phiF[ipt][3]<<endl;
      //cout<<" "<<endl;
    }  

    
    for(Int_t ig=0;ig<4;ig++){
      THh_thetaF[ig]=h_thetaF[0][ig]*(1-fabs(mytheta-dou_theta[0])/(dou_theta[0]-dou_theta[1]))+h_thetaF[1][ig]*(1-fabs(mytheta-dou_theta[1])/(dou_theta[0]-dou_theta[1]));
      THh_phiF[ig]=h_phiF[0][ig]*(1-fabs(mytheta-dou_theta[0])/(dou_theta[0]-dou_theta[1]))+h_phiF[1][ig]*(1-fabs(mytheta-dou_theta[1])/(dou_theta[0]-dou_theta[1]));
      if(ig<4){
	output[ig]=THh_thetaF[ig];
      }
      if(ig>3){
	output[ig]=THh_phiF[ig-4];
      }
    }

    THG_thetaF[0]=G_thetaF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_thetaF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
    THG_phiF[0]=G_phiF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phiF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

     output[8]=THG_thetaF[0];
     output[9]=THG_phiF[0];
    
    //cout<<" "<<THh_thetaF[0]<<" "<<THh_thetaF[1]<<" "<<THh_thetaF[2]<<" "<<THh_thetaF[3]<<endl;
    //cout<<" "<<THh_phiF[0]<<" "<<THh_phiF[1]<<" "<<THh_phiF[2]<<" "<<THh_phiF[3]<<endl;

  }

  /////////////////////////////////////////////////////////////////////////////
  if(fmod(mytheta,thetastp)!=0 && fmod(myphi,phistp)!=0){
    //cout<<"none are divisible"<<endl;

    //cout<<(myphi-fmod(myphi,phistp))/phistp<<" "<<myfreq<<endl;

    Float_t dou_phi[2];
    dou_phi[0]=(myphi-fmod(myphi,phistp))+phistp;
    dou_phi[1]=(myphi-fmod(myphi,phistp));
    
    //cout<<(mytheta-fmod(mytheta,thetastp))/thetastp<<endl;
    
    Float_t dou_theta[2];
    dou_theta[0]=(mytheta-fmod(mytheta,thetastp))+thetastp;
    dou_theta[1]=(mytheta-fmod(mytheta,thetastp));

    if(mytheta<0){
      dou_theta[1]=(mytheta-fmod(mytheta,thetastp))-thetastp;
      dou_theta[0]=(mytheta-fmod(mytheta,thetastp));
    }
    //cout<<dou_theta[0]<<" "<<dou_theta[1]<<endl;
    
    Double_t freqdiff1=0;
    Double_t freqdiff2=10000;
    Double_t abvfreq=0;
    Double_t belfreq=0;
    Int_t abvfreqbin=0;
    Int_t belfreqbin=0;

    AntTree->GetEntry(0);
    //find closest freq
    for (Int_t i=0; i<N; i++) {
      freqdiff1=fabs(freq[i]*pow(10,3)-myfreq);
      if (freqdiff1<freqdiff2){
  	freqdiff2=freqdiff1;
  	abvfreq=freq[i]*pow(10,3);
  	abvfreqbin=i;
      }
    }
    
    if(myfreq<abvfreq){
      abvfreq=abvfreq;
      belfreq=freq[abvfreqbin-1]*pow(10,3);
      abvfreqbin=abvfreqbin;
      belfreqbin=abvfreqbin-1;
    }
    
    if(myfreq>abvfreq){
      belfreq=abvfreq;
      abvfreq=freq[abvfreqbin+1]*pow(10,3);
      belfreqbin=abvfreqbin;
      abvfreqbin=abvfreqbin+1;
    }
    //cout<<abvfreq<<" "<<belfreq<<" "<<abvfreqbin<<" "<<belfreqbin<<endl; 

    for(Int_t iptR=0;iptR<2;iptR++){
      for(Int_t ipt=0;ipt<2;ipt++){
	entry = AntTree->GetEntryNumberWithIndex(dou_theta[iptR],dou_phi[ipt]);
	AntTree->GetEntry(entry);
	//cout<<"entry is "<<entry<<" "<<dou_theta[iptR]<<" "<<dou_phi[ipt]<<endl;

	//calculate heff for abv freq.
	V_theta[0] = TComplex(v_re_theta[abvfreqbin], v_im_theta[abvfreqbin]);
	V_phi[0] = TComplex(v_re_phi[abvfreqbin], v_im_phi[abvfreqbin]);
	
	wavelength = 0.299792/((Double_t) freq[abvfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
	norm = TComplex(0., -2 * wavelength/z0/v0);
	
	Z[abvfreqbin] = TComplex(z_re[abvfreqbin], z_im[abvfreqbin]);
	h_thetaC[0] = norm * Z[abvfreqbin] * V_theta[0];
	h_phiC[0]   = norm * Z[abvfreqbin] * V_phi[0];

	A_eff_theta = pow(TComplex::Abs(h_thetaC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
	A_eff_phi   = pow(TComplex::Abs(h_phiC[0]),2) * 119.9169*pi /(4*z_re[abvfreqbin]);
	G_theta[0] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
	G_phi[0] = 4. * pi * A_eff_phi / (wavelength * wavelength);
	
	//calculate heff for low freq.
	V_theta[1] = TComplex(v_re_theta[belfreqbin], v_im_theta[belfreqbin]);
	V_phi[1] = TComplex(v_re_phi[belfreqbin], v_im_phi[belfreqbin]);
	
	wavelength = 0.299792/((Double_t) freq[belfreqbin]); // speed of light is 3e8 m/s or 0.3 m/ns
	norm = TComplex(0., -2 * wavelength/z0/v0);
	
	Z[belfreqbin] = TComplex(z_re[belfreqbin], z_im[belfreqbin]);
	h_thetaC[1] = norm * Z[belfreqbin] * V_theta[1];
	h_phiC[1]   = norm * Z[belfreqbin] * V_phi[1];

	A_eff_theta = pow(TComplex::Abs(h_thetaC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
	A_eff_phi   = pow(TComplex::Abs(h_phiC[1]),2) * 119.9169*pi /(4*z_re[belfreqbin]);
	G_theta[1] = 4. * pi * A_eff_theta /(wavelength * wavelength); 
	G_phi[1] = 4. * pi * A_eff_phi / (wavelength * wavelength);
	
	//calculate heff individual components for up and low
	for(Int_t iud=0;iud<2;iud++){
	  h_theta[iud][0] = (Double_t) h_thetaC[iud].Re();
	  h_theta[iud][1] = (Double_t) h_thetaC[iud].Im();
	  h_theta[iud][2] = (Double_t) TComplex::Abs(h_thetaC[iud]);
	  h_theta[iud][3] = (Double_t) h_thetaC[iud].Theta();//180./TMath::Pi();
	  
	  h_phi[iud][0] = (Double_t) h_phiC[iud].Re();
	  h_phi[iud][1] = (Double_t) h_phiC[iud].Im();
	  h_phi[iud][2] = (Double_t) TComplex::Abs(h_phiC[iud]);
	  h_phi[iud][3] = (Double_t) h_phiC[iud].Theta();//180./TMath::Pi();
	  
	  // cout<<" "<<h_theta[iud][0]<<" "<<h_theta[iud][1]<<" "<<h_theta[iud][2]<<" "<<h_theta[iud][3]<<endl;
	  // cout<<" "<<h_phi[iud][0]<<" "<<h_phi[iud][1]<<" "<<h_phi[iud][2]<<" "<<h_phi[iud][3]<<endl;
	}
	
	//calculate the weighted average
	for(Int_t ig=0;ig<4;ig++){
	  h_thetaF[ipt][ig]=h_theta[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_theta[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
	  h_phiF[ipt][ig]=h_phi[0][ig]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+h_phi[1][ig]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
	}

	G_thetaF[ipt]=G_theta[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_theta[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
	G_phiF[ipt]=G_phi[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phi[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));
	
	//cout<<" "<<h_thetaF[ipt][0]<<" "<<h_thetaF[ipt][1]<<" "<<h_thetaF[ipt][2]<<" "<<h_thetaF[ipt][3]<<endl;
	//cout<<" "<<h_phiF[ipt][0]<<" "<<h_phiF[ipt][1]<<" "<<h_phiF[ipt][2]<<" "<<h_phiF[ipt][3]<<endl;
	//cout<<" "<<endl;
      }///ipt loop

      for(Int_t ig=0;ig<4;ig++){
	PHh_thetaF[iptR][ig]=h_thetaF[0][ig]*(1-fabs(myphi-dou_phi[0])/(dou_phi[0]-dou_phi[1]))+h_thetaF[1][ig]*(1-fabs(myphi-dou_phi[1])/(dou_phi[0]-dou_phi[1]));
	PHh_phiF[iptR][ig]=h_phiF[0][ig]*(1-fabs(myphi-dou_phi[0])/(dou_phi[0]-dou_phi[1]))+h_phiF[1][ig]*(1-fabs(myphi-dou_phi[1])/(dou_phi[0]-dou_phi[1]));
      }

      PHG_thetaF[iptR]=G_thetaF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_thetaF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
      PHG_phiF[iptR]=G_phiF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+G_phiF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

      
      //cout<<" "<<PHh_thetaF[iptR][0]<<" "<<PHh_thetaF[iptR][1]<<" "<<PHh_thetaF[iptR][2]<<" "<<PHh_thetaF[iptR][3]<<endl;
      //cout<<" "<<PHh_phiF[iptR][0]<<" "<<PHh_phiF[iptR][1]<<" "<<PHh_phiF[iptR][2]<<" "<<PHh_phiF[iptR][3]<<endl;
      
    }///iptR loop
       
    for(Int_t ig=0;ig<4;ig++){
      TPh_thetaF[ig]=PHh_thetaF[0][ig]*(1-fabs(mytheta-dou_theta[0])/(dou_theta[0]-dou_theta[1]))+PHh_thetaF[1][ig]*(1-fabs(mytheta-dou_theta[1])/(dou_theta[0]-dou_theta[1]));
      TPh_phiF[ig]=PHh_phiF[0][ig]*(1-fabs(mytheta-dou_theta[0])/(dou_theta[0]-dou_theta[1]))+PHh_phiF[1][ig]*(1-fabs(mytheta-dou_theta[1])/(dou_theta[0]-dou_theta[1]));
      if(ig<4){
	output[ig]=TPh_thetaF[ig];
      }
      if(ig>3){
	output[ig]=TPh_phiF[ig-4];
      }
 
    }

     TPG_thetaF=PHG_thetaF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+PHG_thetaF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq)); 
     TPG_phiF=PHG_phiF[0]*(1-fabs(myfreq-abvfreq)/(abvfreq-belfreq))+PHG_phiF[1]*(1-fabs(myfreq-belfreq)/(abvfreq-belfreq));

     output[8]=THG_thetaF[0];
     output[9]=THG_phiF[0];     
    
    //cout<<" "<<TPh_thetaF[0]<<" "<<TPh_thetaF[1]<<" "<<TPh_thetaF[2]<<" "<<TPh_thetaF[3]<<endl;
    //cout<<" "<<TPh_phiF[0]<<" "<<TPh_phiF[1]<<" "<<TPh_phiF[2]<<" "<<TPh_phiF[3]<<endl;

  }

  // delete AntTree;
  // delete ZTree;
  // delete ariannaAntennaFile;
  
  return output;  
}
