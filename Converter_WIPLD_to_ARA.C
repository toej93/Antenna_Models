{
	float theta = 0;
	float phi = 0;
	int N = 199;
	TFile *ariannaAntennaFile = new TFile("WIPLD_ARA_bicone_v1.root");
	TTree *AntTree = (TTree*)ariannaAntennaFile->Get("AntTree");
	TTree *ZTree = (TTree*)ariannaAntennaFile->Get("ZTree");
	
	float freq[239];
	double gains[239];
	double v_re_phi[239];
	double v_im_phi[239];
	double v_re_theta[239];
	double v_im_theta[239];
	double thetas[239];
	double phis[239];
	double z_re[239];
	double z_im[239];
	
	TComplex V_theta[239];
	TComplex V_phi[239];
	TComplex Z[239];
	TComplex h_theta[239];
	TComplex h_phi[239];	

	float h_phi_re[239];
	float h_phi_im[239];
	float h_theta_re[239];
	float h_theta_im[239];

	float h_phi_mag[239];
	float h_phi_phase[239];
	float h_theta_mag[239];
	float h_theta_phase[239];

	float Z_mag[239];
	float V_theta_mag[239];
	float V_theta_phase[239];
	float V_phi_mag[239];
	float V_phi_phase[239];

	float A_eff_theta[239];
	float A_eff_phi[239];
	float G_theta[239];
	float G_phi[239];

	float pi = TMath::Pi();

	AntTree->SetBranchAddress("frequencies", &freq);
	AntTree->SetBranchAddress("gains", &gains);
	AntTree->SetBranchAddress("Re_phi", &v_re_phi);
	AntTree->SetBranchAddress("Im_phi", &v_im_phi);
	AntTree->SetBranchAddress("Re_theta", &v_re_theta);
	AntTree->SetBranchAddress("Im_theta", &v_im_theta);
	AntTree->BuildIndex("thetas","phis");

	// fill the antenna impedance. 
	// This tree only has one entry.
	ZTree->SetBranchAddress("Im_Z", &z_im);
	ZTree->SetBranchAddress("Re_Z", &z_re);
	ZTree->GetEntry(0);
	for(int i=0; i<N; i++){
		Z[i] = TComplex(z_re[i], z_im[i]);
		//		std:: cout << i << " " << z_re[i] << " " << z_im[i] << " " << Z[i] << std::endl;
	}
	TGraph *gV_phi_re = new TGraph(N);
	TGraph *gV_phi_im = new TGraph(N);
	TGraph *gV_theta_re = new TGraph(N);
	TGraph *gV_theta_im = new TGraph(N);
	
	Int_t entry = AntTree->GetEntryNumberWithIndex(theta,phi);
	AntTree->GetEntry(entry);
	
	

	for(int i=0; i<N; i++){
		h_phi_re[i] = 0.;
		h_phi_im[i] = 0.;
		h_theta_re[i] = 0;
		h_theta_im[i] = 0;
		
		V_theta[i] = TComplex(v_re_theta[i], v_im_theta[i]);
		V_phi[i] = TComplex(v_re_phi[i], v_im_phi[i]);

		double z0 = 377.;
		double v0 = 1.0;
		double wavelength = 0.299792/((double) freq[i]); // speed of light is 3e8 m/s or 0.3 m/ns
		TComplex norm = TComplex(0., -2 * wavelength/z0/v0);

		//	TComplex norm = TComplex(0., -2 * wavelength/z0/v0);
		h_theta[i] = norm * Z[i] * V_theta[i];
		//	TComplex norm = TComplex(0., -2 * wavelength/z0/v0);
		//	TComplex norm = TComplex(0., -2 * wavelength/z0/v0);
		h_phi[i]   = norm * Z[i] * V_phi[i];
		//	std::cout << norm << std::endl;
				
		
		gV_phi_re->SetPoint(i, freq[i], (float) v_re_phi[i]);
		gV_theta_re->SetPoint(i, freq[i], (float) v_re_theta[i]);
		gV_phi_im->SetPoint(i, freq[i], (float) v_im_phi[i]);
		gV_theta_im->SetPoint(i, freq[i], (float) v_im_theta[i]);


		h_theta_re[i] = (float) h_theta[i].Re();
		h_theta_im[i] = (float) h_theta[i].Im();
		h_theta_mag[i] = (float) TComplex::Abs(h_theta[i]);
		h_theta_phase[i] = (float) h_theta[i].Theta();//*180./TMath::Pi();

		h_phi_re[i] = (float) h_phi[i].Re();
		h_phi_im[i] = (float) h_phi[i].Im();
		h_phi_mag[i] = (float) TComplex::Abs(h_phi[i]);
		h_phi_phase[i] = (float) h_phi[i].Theta();//*180./TMath::Pi();

		Z_mag[i] = TComplex::Abs(Z[i]);
		V_phi_mag[i] = TComplex::Abs(V_phi[i]);
		V_theta_mag[i] = TComplex::Abs(V_theta[i]);
		V_phi_phase[i] = V_phi[i].Theta();
		V_theta_phase[i] = V_theta[i].Theta();

		// from gain=4*pi*A_eff/lambda^2
    		// and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
    		// gain is unitless value
    		//return 2*sqrt(gain/4/PI*CLIGHT*CLIGHT/(freq*freq*n_medium*n_medium)*Zr/(Z0/n_medium));  // n_medium parts are changed from icemc(I believe this is correct one; E. Hong)
		
		A_eff_theta[i] = pow(h_theta_mag[i],2) * 119.9169*pi /(4*z_re[i]);
		A_eff_phi[i]   = pow(h_phi_mag[i],2) * 119.9169*pi /(4*z_re[i]);
		G_theta[i] = 4. * pi * A_eff_theta[i] /(wavelength * wavelength); 
		G_phi[i] = 4. * pi * A_eff_phi[i] / (wavelength * wavelength);
    
		//		std::cout << freq[i] << " " << norm << " "<< Z[i].Re() << " " << z_re[i] << " " << Z[i].Im() << " " << z_im[i] << " "  << V_theta[i].Re() << " " << V_theta[i].Im() << " " << h_theta[i].Re() << " " << h_theta[i].Im()<< std::endl;
	}

	// voltages from the AntTree
		/*		
	TCanvas *cvolts_re = new TCanvas();
	gV_phi_re->SetTitle(Form("V_{#phi,real}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gV_theta_re->SetTitle(Form("V_{#theta,real}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gV_phi_re->SetLineColor(kRed);
	gV_phi_re->SetFillColor(0);
	gV_theta_re->SetLineColor(kBlack);
	gV_theta_re->SetFillColor(0);
	gV_phi_re->Draw("AL");
	gV_theta_re->Draw("LSAME");

	
	TCanvas *cvolts_im = new TCanvas();
	gV_phi_im->SetTitle(Form("V_{#phi,im}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gV_theta_im->SetTitle(Form("V_{#theta,im}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gV_phi_im->SetLineColor(kRed);
	gV_phi_im->SetFillColor(0);
	gV_theta_im->SetLineColor(kBlack);
	gV_theta_im->SetFillColor(0);
	gV_phi_im->Draw("AL");
	gV_theta_im->Draw("LSAME");
		*/

	// Effective height graphs
	TGraph *gheff_phi_re = new TGraph(N, freq, h_phi_re);
	TGraph *gheff_theta_re = new TGraph(N, freq, h_theta_re);

	TGraph *gheff_phi_im = new TGraph(N, freq, h_phi_im);
	TGraph *gheff_theta_im = new TGraph(N, freq, h_theta_im);


	TGraph *gheff_phi_mag = new TGraph(N, freq, h_phi_mag);
	TGraph *gheff_theta_mag = new TGraph(N, freq, h_theta_mag);

	TGraph *gheff_phi_phase = new TGraph(N, freq, h_phi_phase);
	TGraph *gheff_theta_phase = new TGraph(N, freq, h_theta_phase);

	TGraph *gA_eff_phi = new TGraph(N, freq, A_eff_phi);
	TGraph *gA_eff_theta = new TGraph(N, freq, A_eff_theta);

	TGraph *gG_phi = new TGraph(N, freq, G_phi);
	TGraph *gG_theta = new TGraph(N, freq, G_theta);
	/*
	TCanvas *creal = new TCanvas();
	gheff_phi_re->SetTitle(Form("h_{#phi,real}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gheff_theta_re->SetTitle(Form("h_{#theta,real}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gheff_phi_re->SetLineColor(kRed);
	gheff_phi_re->SetFillColor(0);
	gheff_theta_re->SetLineColor(kBlack);
	gheff_theta_re->SetFillColor(0);
	gheff_phi_re->Draw("AL");
	gheff_theta_re->Draw("LSAME");

	TCanvas *cim = new TCanvas();
	gheff_phi_im->SetTitle(Form("h_{#phi, imag}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gheff_theta_im->SetTitle(Form("h_{#theta, imag}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gheff_phi_im->SetLineColor(kRed);
	gheff_phi_im->SetFillColor(0);
	gheff_theta_im->SetLineColor(kBlack);
	gheff_theta_im->SetFillColor(0);
	gheff_phi_im->Draw("AL");
	gheff_theta_im->Draw("LSAME");
	*/
	TCanvas *cmag = new TCanvas();
	gheff_phi_mag->SetTitle(Form("h_{#phi,mag}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gheff_theta_mag->SetTitle(Form("h_{#theta,mag}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gheff_phi_mag->SetLineColor(kRed);
	gheff_phi_mag->SetFillColor(0);
	gheff_theta_mag->SetLineColor(kBlue);
	gheff_theta_mag->SetFillColor(0);
	gheff_theta_mag->GetXaxis()->SetTitle("Frequency [GHz]");
	gheff_theta_mag->GetYaxis()->SetTitle("Effective height [m]");
	gheff_theta_mag->SetLineWidth(2);
	//	gheff_phi_mag->Draw("AL");
	gheff_theta_mag->Draw("AL");
	/*
	TCanvas *cphase = new TCanvas();
	gheff_phi_phase->SetTitle(Form("h_{#phi,phase}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gheff_theta_phase->SetTitle(Form("h_{#theta,phase}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gheff_phi_phase->SetLineColor(kRed);
	gheff_phi_phase->SetFillColor(0);
	gheff_theta_phase->SetLineColor(kBlack);
	gheff_theta_phase->SetFillColor(0);
	gheff_phi_phase->Draw("AL");
	gheff_theta_phase->Draw("LSAME");
	
	//	TGraph *gG_phi = new TGraph(N, freq, G_phi);
	//	TGraph *gG_theta = new TGraph(N, freq, G_theta);
	TCanvas *cgain = new TCanvas();
	gG_phi->SetTitle(Form("G_{#phi}(#theta=%d, #phi=%d)", (int)theta, (int)phi));
	gG_theta->SetTitle(Form("G_{#theta}(#theta=%d, #phi=%d)",(int) theta,(int) phi));
	gG_phi->SetLineColor(kBlue);
	gG_phi->SetFillColor(0);
	gG_theta->SetLineColor(kRed);
	gG_theta->SetFillColor(0);
	gG_phi->Draw("AL");
	gG_phi->GetXaxis()->SetTitle("Frequency [GHz]");
	gG_phi->GetYaxis()->SetTitle("Gain [dB]");
	gG_theta->Draw("LSAME");
	AntTree->Draw("gains:frequencies", Form("thetas==%d && phis==%d", (int)theta, (int)phi), "*SAME");
        

	leg = new TLegend(0.2,0.2,0.4,0.4);
	//	leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	leg->AddEntry(gG_phi,"G_{#phi}");
	leg->AddEntry(gG_theta,"G_{#theta}");
	//	leg->AddEntry(AntTree,"Data from tree","*P");
	TLegendEntry *le = leg->AddEntry(AntTree,"Data from tree","*P");
	le->SetMarkerStyle(21);
	le->SetMarkerSize(10);
	leg->Draw("SAME");
		*/


	
	// WIPL-D Model graphs
	TGraph *gzant = new TGraph(N,freq, Z_mag);
	TGraph *gvphi_mag = new TGraph(N,freq, V_phi_mag);
	TGraph *gvtheta_mag = new TGraph(N,freq, V_theta_mag);
	TGraph *gvphi_phase = new TGraph(N,freq, V_phi_phase);
	TGraph *gvtheta_phase = new TGraph(N,freq, V_theta_phase);

	gvphi_mag->SetTitle("|I_{#phi, mag}|");
	gvtheta_mag->SetTitle("|I_{#theta, mag}|");
	gvphi_phase->SetTitle("|I_{#phi, phase}|");
	gvtheta_phase->SetTitle("|I_{#theta, phase}|");

	TCanvas *czant = new TCanvas();
	gzant->SetTitle("Z_{bicone (WIPL-D)}");
	gzant->GetXaxis()->SetTitle("Frequency [GHz]");
	gzant->GetYaxis()->SetTitle("Antenna Impedance [#Omega]");
	gzant->SetLineColor(kBlue);
	gzant->SetFillColor(0);
	gzant->SetLineWidth(2);
	gzant->Draw("AL");
	/*
	TCanvas *cVmag= new TCanvas();
	gvphi_mag->SetLineColor(kRed);
	gvphi_mag->SetFillColor(0);
	gvtheta_mag->SetLineColor(kBlack);
	gvtheta_mag->SetFillColor(0);
	gvphi_mag->Draw("AL");
	gvtheta_mag->Draw("LSAME");

	TCanvas *cVphase= new TCanvas();
	gvphi_phase->SetLineColor(kRed);
	gvphi_phase->SetFillColor(0);
	gvtheta_phase->SetLineColor(kBlack);
	gvtheta_phase->SetFillColor(0);
	gvphi_phase->Draw("AL");
	gvtheta_phase->Draw("LSAME");
	std::cout << "done" << std::endl;
	*/	

	/*for(int p=0; p<=360; p+=5){
		for(int t=-90; t<=90; t+=5){
			entry = AntTree->GetEntryNumberWithIndex(t,p);
			if( n>=-1){
				AntTree->GetEntry(entry);
			}
		}
	}*/
	
	/*	
       	// Now that everything works, step through the frequencies 
	// and print out the gain and phase at each phi and theta angle
	TComplex norm = TComplex(0., -2. * wavelength/z0/v0);
	FILE *fout = fopen("ARA_bicone.txt", "w");
	for( int i=1; i<N; i++){
		
		fprintf(fout, "freq : %3.2f MHz\n", freq[i]*1000.);
		fprintf(fout, "SWR : %f\n", 2-0.007*i);
		fprintf(fout, " Theta \t Phi \t Gain(dB)     \t   Gain     \t    Phase(deg)\n");
	
		for(int p=-90;p<=265; p+=5){
		  double z0 = 377.;
		  double v0 = 1.0;
		  double wavelength = 0.299792/((double) freq[i]); // speed of light is 3e8 m/s or 0.3 m/ns
		  for(int t=-90;t<=90; t+=5){
		    Int_t en = AntTree->GetEntryNumberWithIndex(t,p);
		    AntTree->GetEntry(en);
		    
		    h_phi_re[i] = 0.;
		    h_phi_im[i] = 0.;
		    h_theta_re[i] = 0.;
		    h_theta_im[i] = 0.;
		    
		    V_theta[i] = TComplex(v_re_theta[i], v_im_theta[i]);
		    V_phi[i] = TComplex(v_re_phi[i], v_im_phi[i]);
		    
		    h_theta[i] = norm * Z[i] * V_theta[i];
		    //  TComplex norm = TComplex(0., -2. * wavelength/z0/v0);
		    h_phi[i]   = norm * Z[i] * V_phi[i];
		  //    TComplex norm = TComplex(0., -2. * wavelength/z0/v0);
		  								
		    
		    h_theta_re[i] = (float) h_theta[i].Re();
		    h_theta_im[i] = (float) h_theta[i].Im();
		    h_theta_mag[i] = (float) TComplex::Abs(h_theta[i]);
		    h_theta_phase[i] = (float) h_theta[i].Theta();//*180./TMath::Pi();
		    
		    h_phi_re[i] = (float) h_phi[i].Re();
		    h_phi_im[i] = (float) h_phi[i].Im();
		    h_phi_mag[i] = (float) TComplex::Abs(h_phi[i]);
		    h_phi_phase[i] = (float) h_phi[i].Theta();//*180./TMath::Pi();
		    
		    
		    fprintf(fout, "%d \t %d \t %2.2f     \t   %2.2f     \t    %3.2f\n", t+90, p+90, 10.*TMath::Log10(gains[i]), gains[i],  h_theta_phase[i]*TMath::RadToDeg() ); 
		    //printf("%d\t%d\t%d\t%2.2f\t%2.2f\t%3.2f\n",en, t, p, 20.*TMath::Log10(gains[i]), gains[i],  h_theta_phase[i] );
		  }
 		}
		std::cout << i  << std::endl;
	} 
	fclose(fout);
	std::cout << "done" << std::endl;
	}*/	

}


//setw
