//#include "readgainVH.C"
#include "WIPLD_heff.C"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "TMath.h"


void writer(){
  FILE *fout = fopen("ARIANNA_LPDA_firn.txt", "w");
  for( Double_t ifr=83.333;ifr<1066.70;ifr+=16.667){	
    fprintf(fout, "freq : %3.2f MHz\n", ifr);
    fprintf(fout, "SWR : %f\n", 2);
    //		fprintf(fout, "Impedance : %f\n", TComplex::Abs(Z[i]));
    fprintf(fout, " Theta \t Phi \t Gain(dB)     \t   Gain     \t    Phase(deg)\n");
    for(int p=-90;p<=270; p+=5){//Loop over phi
      for(int t=-90;t<=90; t+=5){//Loop over theta

	Double_t *getres=WIPLD_heff(t,p,ifr);
	Double_t heff=0;
	heff=sqrt(getres[2]*getres[2]+getres[6]*getres[6]);
	fprintf(fout, "%d \t %d \t %2.2f     \t   %2.2f     \t    %3.2f\n", t+90, p+90, heff, heff, heff ); 
	//printf("%d\t%d\t%d\t%2.2f\t%2.2f\t%3.2f\n",en, t, p, 20.*TMath::Log10(gains[i]), gains[i],  h_theta_phase[i] );
      }
    }
    std::cout << ifr  << std::endl;
  }

}
