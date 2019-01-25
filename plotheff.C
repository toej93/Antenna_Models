#include "readgainVH.C"
#include "WIPLD_heff.C"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "TMath.h"


void plotheff(){
  
  Double_t heff=0;
  TNtuple *ntheff=new TNtuple("","","ifre:hff");
  TNtuple *ntPH=new TNtuple("","","ifre:PH");

  TFile *ariannaFile = new TFile("WIPLD_antennamodel_firn.root");

  TTree *ATree = (TTree*)ariannaFile->Get("AntTree");
  TTree *zTree = (TTree*)ariannaFile->Get("ZTree");
 
  for(Double_t ifr=10;ifr<1000;ifr=ifr+5){

    Double_t *getres=WIPLD_heff(90.0,90.0,ifr,ariannaFile,ATree,zTree);
    heff=sqrt(getres[2]*getres[2]+getres[6]*getres[6]);
    ntheff->Fill(ifr,heff);
    ntPH->Fill(ifr,getres[3]);
    cout<<ifr<<" "<<heff<<endl;
    cout<<"G_theta "<<getres[8]<<" ,G_phi "<<getres[9]<<endl;

    delete getres;
  }
  TCanvas * c1=new TCanvas("c1","c1");
  c1->Divide(1,2);
  c1->cd(1);
  ntPH->Draw("PH:ifre","","PL");
  c1->cd(2);
  ntheff->Draw("hff:ifre","","PL");
  
  // Double_t *getres=WIPLD_heff(25.0,35.0,400);
  // heff=sqrt(getres[2]*getres[2]+getres[6]*getres[6]);
  // cout<<heff<<endl;
  
}
