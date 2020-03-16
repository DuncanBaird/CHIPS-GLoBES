#include <TROOT.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TTree.h>
#include <TMath.h>
#include <THStack.h>
#include <fstream>

using namespace std;

/*void readfiles_compare(TString filename1, TString filename2, TH1D * h_chi_1, TH1D * h_chi_2, TH1D * h_sig){
	ifstream infile;
	infile.open(filename1);


}*/

void readfile_makeChiplot(TString filename,TH1D * h_chi_1, TH1D * h_chi_2, TH1D * h_sig){
	ifstream infile;
	infile.open(filename);
	if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl;  // no point continuing if the file didn't open...
    } 

	int i = 0;
	int bin;
	double d_bin, chiSq1, chiSq2; 
  // int numBin = 1000;
	while(!infile.eof()){
		infile >> bin >> chiSq1 >> chiSq2;

		h_chi_1->SetBinContent(bin, chiSq1);
		h_chi_2->SetBinContent(bin, chiSq2);

		if(chiSq1 != 0. || chiSq2 != 0.){
			if(chiSq1 > chiSq2){
				h_sig->SetBinContent(bin, TMath::Sqrt(TMath::Sqrt(TMath::Power(chiSq2,2))));
			} else{
				h_sig->SetBinContent(bin, TMath::Sqrt(TMath::Sqrt(TMath::Power(chiSq1,2))));
			}
		}
	else{ //if((chiSq1 != 0. && chiSq2 != 0.)){
		h_sig->SetBinContent(bin, 0.0);
	}

	++i;

}

infile.close();
}


void readfile_makeSmallErrorplot(TString filename,  TH1D * h_pi0, TH1D * h_pi1, TH1D * h_pi2, TH1D * h_zero0, TH1D * h_zero1, TH1D * h_zero2, bool opt){
	ifstream infile;
	infile.open(filename);

	int i = 0;
	int bin;
	double pi_0_a_0, pi_0_error_0, pi_0_a_1, pi_0_error_1,pi_0_a_2, pi_0_error_2;
	double zero_0_a_0 , zero_0_error_0, zero_0_a_1, zero_0_error_1,zero_0_a_2, zero_0_error_2; 
	double pi_1_a_0, pi_1_error_0, pi_1_a_1, pi_1_error_1,pi_1_a_2, pi_1_error_2;
	double zero_1_a_0 , zero_1_error_0, zero_1_a_1, zero_1_error_1,zero_1_a_2, zero_1_error_2;


	while(!infile.eof()){
		//infile >> bin >> pi_a_0 >> pi_error_0 >> pi_a_1 >> pi_error_1 >> pi_a_2 >> pi_error_2 >> zero_a_0 >> zero_error_0 >>
		//		zero_a_1 >> zero_error_1 >> zero_a_2 >> zero_error_2;

		
		infile >>  bin >> pi_0_a_0 >> pi_0_error_0 >> pi_0_a_1 >> pi_0_error_1 >>pi_0_a_2 >> pi_0_error_2 >>
		pi_1_a_0 >> pi_1_error_0 >> pi_1_a_1 >> pi_1_error_1 >>pi_1_a_2 >> pi_1_error_2 >>
		zero_0_a_0  >> zero_0_error_0 >> zero_0_a_1 >> zero_0_error_1 >>zero_0_a_2 >> zero_0_error_2 >>
		zero_1_a_0  >> zero_1_error_0 >> zero_1_a_1 >> zero_1_error_1 >>zero_1_a_2 >> zero_1_error_2;


	 	if(opt == true){  // for exp near
	 		h_pi0->SetBinContent(bin, TMath::Power(pi_0_a_0 / pi_0_error_0, 2));
	 		h_pi1->SetBinContent(bin, TMath::Power(pi_0_a_1 / pi_0_error_1, 2));
	 		h_pi2->SetBinContent(bin, TMath::Power(pi_0_a_2 / pi_0_error_2, 2));

	 		h_zero0->SetBinContent(bin, TMath::Power(zero_0_a_0 / zero_0_error_0, 2));
	 		h_zero1->SetBinContent(bin, TMath::Power(zero_0_a_1 / zero_0_error_1, 2));
	 		h_zero2->SetBinContent(bin, TMath::Power(zero_0_a_2 / zero_0_error_2, 2));
	 	}
		else if(opt != true){ // for exp far
			h_pi0->SetBinContent(bin, TMath::Power(pi_1_a_0 / pi_1_error_0, 2));
			h_pi1->SetBinContent(bin, TMath::Power(pi_1_a_1 / pi_1_error_1, 2));
			h_pi2->SetBinContent(bin, TMath::Power(pi_1_a_2 / pi_1_error_2, 2));

			h_zero0->SetBinContent(bin, TMath::Power(zero_1_a_0 / zero_1_error_0, 2));
			h_zero1->SetBinContent(bin, TMath::Power(zero_1_a_1 / zero_1_error_1, 2));
			h_zero2->SetBinContent(bin, TMath::Power(zero_1_a_2 / zero_1_error_2, 2));
		}
/*
		h_pi0->SetBinContent(bin, pi_a_0);
		h_pi1->SetBinContent(bin, pi_a_1);
		h_zero0->SetBinContent(bin, zero_a_0);
		h_zero1->SetBinContent(bin, zero_a_1);
*/
		++i;
	}
	infile.close();
}

void readfile_makeLargeErrorplot(TString filename,  TH1D * h_pi0, TH1D * h_pi1, TH1D * h_pi2, TH1D * h_zero0, TH1D * h_zero1, TH1D * h_zero2){
	ifstream infile;
	infile.open(filename);

	int i = 0;
	int bin;
	double pi_a_0, pi_error_0, pi_a_1, pi_error_1, pi_a_2, pi_error_2;
	double zero_a_0, zero_error_0, zero_a_1, zero_error_1, zero_a_2, zero_error_2;


	while(!infile.eof()){
		infile >> bin >> pi_a_0 >> pi_error_0 >> pi_a_1 >> pi_error_1 >> pi_a_2 >> pi_error_2 >> zero_a_0 >> zero_error_0 >>
		zero_a_1 >> zero_error_1 >> zero_a_2 >> zero_error_2;

		



		h_pi0->SetBinContent(bin, TMath::Power(pi_a_0 / pi_error_0, 2));
		h_pi1->SetBinContent(bin, TMath::Power(pi_a_1 / pi_error_1, 2));
		h_pi2->SetBinContent(bin, TMath::Power(pi_a_2 / pi_error_2, 2));

		h_zero0->SetBinContent(bin, TMath::Power(zero_a_0 / zero_error_0, 2));
		h_zero1->SetBinContent(bin, TMath::Power(zero_a_1 / zero_error_1, 2));
		h_zero2->SetBinContent(bin, TMath::Power(zero_a_2 / zero_error_2, 2));
		

		++i;
	}
	infile.close();
}


int plotter(){
	printf("debug");

	TFile * output = new TFile("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/TomStuff/2S_100kt_1_1y_ES2_FV2_BN2.root", "RECREATE");

	TString chifile =  "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/TomStuff/2S_100kt_1_1y_ES2_FV2_BN2.dat";
	TString errorfile = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/TomStuff/2S_GLoBES_errors_1_1y_ES2_FV2_BN2.txt";


	TH1D * hist_chiPi = new TH1D("hist_chiPi","Significance of #delta_{CP}; #delta_{CP} value; #chi^{2}", 360, 0 , 360);
	TH1D * hist_chiZero = new TH1D("hist_chiZero","", 360, 0 , 360);
	TH1D * hist_signif = new TH1D("hist_signif", "Significance of #delta_{CP} from two 10kt detectors Same Baseline apart 7mrad NuMI 5+5y; value of #delta_{CP} / ^{o}; #sqrt{#chi^{2}_{min}}", 360, 0, 360);




	TH1D * hist_pi0 = new TH1D("hist_pi0","", 360, 0 , 360);
	TH1D * hist_pi1 = new TH1D("hist_pi1","", 360, 0 , 360);
	TH1D * hist_pi2 = new TH1D("hist_pi2","", 360, 0 , 360);
	TH1D * hist_zero0 = new TH1D("hist_zero0","", 360, 0 , 360);
	TH1D * hist_zero1 = new TH1D("hist_zero1","", 360, 0 , 360);
	TH1D * hist_zero2 = new TH1D("hist_zero2","", 360, 0 , 360);

	TH1D * hist_pi0_far = new TH1D("hist_pi0_far","", 360, 0 , 360);
	TH1D * hist_pi1_far = new TH1D("hist_pi1_far","", 360, 0 , 360);
	TH1D * hist_pi2_far = new TH1D("hist_pi2_far","", 360, 0 , 360);
	TH1D * hist_zero0_far = new TH1D("hist_zero0_far","", 360, 0 , 360);
	TH1D * hist_zero1_far = new TH1D("hist_zero1_far","", 360, 0 , 360);
	TH1D * hist_zero2_far = new TH1D("hist_zero2_far","", 360, 0 , 360);

	TH1D * totpi = new TH1D("totpi","",360,0,360);
	TH1D * totzero = new TH1D("totzero","",360,0,360);

	THStack * Total_Cont_pi = new THStack("Total_Cont_pi", "");
	THStack * Total_Cont_zero = new THStack("Total_Cont_zero", "");

	TH1D * pi_likelihoodOnly = new TH1D("pi_likelihoodOnly","",360,0,360);
	TH1D * zero_likelihoodOnly = new TH1D("zero_likelihoodOnly","",360,0,360);

	hist_pi0->SetStats(0);
	hist_pi1->SetStats(0);
	hist_pi2->SetStats(0);
	hist_zero0->SetStats(0);
	hist_zero1->SetStats(0);
	hist_zero2->SetStats(0);

	hist_chiPi->SetStats(0);
	hist_chiZero->SetStats(0);
	hist_signif->SetStats(0);

	hist_signif->SetLineColor(kMagenta-9);
	hist_signif->SetLineWidth(2);
	hist_signif->SetFillColor(kMagenta-10);

	hist_chiPi->SetLineColor(kGreen+2);
	hist_chiZero->SetLineColor(kBlue);

	hist_pi0->SetLineColor(kGreen+2);
	hist_pi1->SetLineColor(kBlue);
	hist_pi2->SetLineColor(kRed);

	hist_zero0->SetLineColor(kGreen+2);
	hist_zero1->SetLineColor(kBlue);
	hist_zero2->SetLineColor(kRed);


	readfile_makeChiplot(chifile, hist_chiPi, hist_chiZero, hist_signif);

	if(chifile.BeginsWith("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/TomStuff/2")){
		readfile_makeSmallErrorplot(errorfile, hist_pi0, hist_pi1,hist_pi2, hist_zero0, hist_zero1, hist_zero2, true);
		readfile_makeSmallErrorplot(errorfile, hist_pi0_far, hist_pi1_far,hist_pi2_far, hist_zero0_far, hist_zero1_far, hist_zero2_far, false);


		hist_pi0_far->SetStats(0);
		hist_pi1_far->SetStats(0);
		hist_pi2_far->SetStats(0);
		hist_zero0_far->SetStats(0);
		hist_zero1_far->SetStats(0);
		hist_zero2_far->SetStats(0);


		hist_pi0_far->SetLineColor(kGreen+2);
		hist_pi1_far->SetLineColor(kBlue);
		hist_pi2_far->SetLineColor(kRed);

		hist_zero0_far->SetLineColor(kGreen+2);
		hist_zero1_far->SetLineColor(kBlue);
		hist_zero2_far->SetLineColor(kRed);
	}else{
		readfile_makeLargeErrorplot(errorfile, hist_pi0, hist_pi1,hist_pi2, hist_zero0, hist_zero1, hist_zero2);
	}

	Total_Cont_pi->Add(hist_pi0);
	Total_Cont_pi->Add(hist_pi1);
	Total_Cont_pi->Add(hist_pi2);

	Total_Cont_zero->Add(hist_zero0);
	Total_Cont_zero->Add(hist_zero1);
	Total_Cont_zero->Add(hist_zero2);


	totpi->Clone("hist_pi0");
	totpi->Add(hist_pi1);
	totpi->Add(hist_pi2);
	totzero->Clone("hist_zero0");
	totzero->Add(hist_zero1);
	totzero->Add(hist_zero2);

	if(chifile.BeginsWith("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/TomStuff/2")){
		totpi->Add(hist_pi0_far);
		totpi->Add(hist_pi1_far);
		totpi->Add(hist_pi2_far);

		totzero->Add(hist_zero0_far);
		totzero->Add(hist_zero1_far);
		totzero->Add(hist_zero2_far);
	}
	

	pi_likelihoodOnly->Add(hist_chiPi, totpi ,1, -1);
	zero_likelihoodOnly->Add(hist_chiZero,totzero,1, -1);

	pi_likelihoodOnly->SetName("pi_likelihoodOnly");
	pi_likelihoodOnly->SetTitle("True CP = pi - Likelihood only");		

	zero_likelihoodOnly->SetName("zero_likelihoodOnly");
	zero_likelihoodOnly->SetTitle("True CP = zero - Likelihood only");	



	output->Write();
	
	return 0;

}
