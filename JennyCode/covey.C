#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

#include <TH2.h>
#include <TFile.h>
#include <TRandom.h>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TApplication.h>



int userConfirm(){
  int kz;
  cout << "Please enter binary confirmation for plotting: ";
  cin >> kz;
  cout << "The value you entered is " << kz << ".\n";
  return kz;
}

/* 
Try to make a correlation matrix

*/

void histInterpolate(TH1F *Euniverse1){
  //Interpolate gaps
  for(int i = 0; i<200; i++){
     double bin_content =  Euniverse1->GetBinContent(i);
     
     if(bin_content == 0){
       
       double left_bin = Euniverse1->GetBinContent(i+1);
       double right_bin = Euniverse1->GetBinContent(i-1);
       if((left_bin != 0) && (right_bin!=0) ){
         Euniverse1->SetBinContent(i,(left_bin+right_bin)/2);
       }
     }
   }
}
/*
* Apply errors to energy bins, specify errors as decimals
* option = 1 for both, 0 for shift only, -1 for neither
*/
void errorApply(double &energy, double calibration, double shift,int option){
  if(option == 1){
    energy = (1-calibration)*energy + shift*energy;
  }else if (option == 0)
  {
    energy = (1)*energy + shift*energy;
  }else if (option == -1)
  {
    energy = energy;
  } 
}

/**
 * Function to create covariance matrix from spectra.
 * Errors; eoption = 1 for both calibration and shift, 0 for shift only, -1 for neither.
 * */
TH1F *createCovariance(TMatrixD &covariance_matrix, TMatrixD &correlation_matrix, int eoption){
  double calibration_e = 0.04;
  double shift_e = 0.05;

  const int mat_len = 200;
  // TMatrixD covariance_matrix(200,200);
  // TMatrixD correlation_matrix(200,200);
  double_t matrix_data[200*200];

  //two detectors, close together, energy from 0-5

   TH2D *mat = new TH2D("mat","mat",200,0.0,10.0,200,0.0,10.0);
   int iuniverse=100;
   char no[5];

   TH1F *Euniverse1[100];
   TRandom *g = new TRandom();
   for (int i=0; i<iuniverse; i++)
     {
       sprintf (no, "no%d",i);
       Euniverse1[i] = new TH1F(no,no,200,0.0,10.0);
     }

   //generate the event spectra


   /*MINOS*/

   TFile *ff3 = new TFile("TrueE.root");
   TH1F *Espec1 = (TH1F*)ff3->Get("trueE_0to5");
   Espec1->SetName("Espec1");

   
   // Adds some floor noise to bottom of MINOS spectrum
   for(int i = 1; i < Espec1->GetXaxis()->FindBin(1);++i){
     if(Espec1->GetBinContent(i)==0){
       Espec1->SetBinContent(i,3.0);
     }
   }


   double c = 2.99792458e8;
   double pi = 3.1415926535;

   double t23 = 45.0/360.0*2.0*pi;
   double sint23 = sin(t23);
   double sin2t23 = sin(2.0*t23);
   double cost23 = cos(t23);
   double sinsqt23 = sint23*sint23;
   double cossqt23 = cost23*cost23;

   //   double t13 = 8.73/360.0*2.0*pi; //Daya Bay
   double t13 = 9.22/360.0*2.0*pi; //Average
   //double t13 = 15.0/360.0*2.0*pi; //code testing

   double sint13 = sin(t13);
   double sin2t13 = sin(t13*2.0);
   double sinsq2t13 = sin2t13*sin2t13;
   double cost13 = cos(t13);
   double cossqt13 = cost13*cost13;

   double t12 = 33.96/360.0*2.0*pi;
   double sin2t12 = sin(2.0*t12);
   double sinsq2t12 = sin2t12*sin2t12;
   double sint12 = sin(t12);
   double sinsqt12 = sint12*sint12;

   double Dmsq31NH = 2.45e-3;  //Sushant's values
   double Dmsq31IH = -2.31e-3;
   double Dmsq21 = 7.6e-5;

   double countchips, countminos, countt2k, countnova;
   double countminoschips, countminosnova, countt2knova, countall;
   double a = 1./4000.0;
   double ab = -1./4000.0;
   double AhatNHn = a/Dmsq31NH;
   double AhatIHn = a/Dmsq31IH;

   double AhatNHb = ab/Dmsq31NH;
   double AhatIHb = ab/Dmsq31IH;

   double alphaNH = Dmsq21/Dmsq31NH;
   double alphaIH = Dmsq21/Dmsq31IH;

   double deltacp=0.0;

   //make the NOVA spectrum
   TH1D *Espec2 = new TH1D("Espec2","Espec2",100,0.0,5.0);
   for (int i=0;i<1E4;i++)
     {
       Espec2->Fill(g->Gaus(2.0,0.3));
     }

    // FLoor noise addition
    for(int i = 1; i < Espec2->GetXaxis()->FindBin(1);++i){
     if(Espec2->GetBinContent(i)==0){
       Espec2->SetBinContent(i,3.0);
     }
    }
    for(int i = Espec2->GetNbinsX(); i > Espec2->GetXaxis()->FindBin(3);--i){
     if(Espec2->GetBinContent(i)==0){
       Espec2->SetBinContent(i,3.0);
     }
    }

   //make the universes based on the real event spectrum
   //std::cout<<" making the universes "<<std::endl;
   for (int k=0; k<100; k++) 
      {
	//use Fill instead of SetBinContent to allow for resolution smearing later
	//clear and reset the spectra histograms for each fake experiment

        for (int i=0; i<100; i++) //energy bins
          {
	    //E1 is MINOS or CHIPS, E2 is NOVA
	     double E1 = Espec1->GetBinCenter(i+1);
	     double Ewt1 = Espec1->GetBinContent(i+1);
	     double E2 = Espec2->GetBinCenter(i+1);
	     double Ewt2 = Espec2->GetBinContent(i+1);
	     //for each energy bin, throw a possion based on no entries
	     float i1 = g->Poisson(Ewt1);
	     float i2 = g->Poisson(Ewt2);
	     //now lets put in a constant E shift

       errorApply(E1,calibration_e,shift_e,eoption);
       errorApply(E2,calibration_e,shift_e,eoption);
	    //  E1 = E1+0.05*E1;
	    //  E2 = E2+0.05*E2;
	     //float rfake = g->Gaus(E,res[0]/sqrt(E));
	     //float Eres = E;//+rfake;
	     //fake11->Fill(Eres,ifake);

	     Euniverse1[k]->Fill(E1,i1);
	     Euniverse1[k]->Fill(E2+5,i2); //E2+5.0


	  } ///finished spectrum generation for this k universe
    histInterpolate(Euniverse1[k]);
//std::cout<<" finished this universe "<<k<<std::endl;
	for (int i=0; i<200; i++)
	  {
	    for (int j=0; j<200; j++)
	      {
		double alan, mary;
		if(i<100)
		  {
		    alan = Espec1->GetBinContent(i+1) - Euniverse1[k]->GetBinContent(i+1);
		  }
		if(j<100)
		  {
		    mary = Espec1->GetBinContent(j+1) - Euniverse1[k]->GetBinContent(j+1);
		  }
		if(i>101)
		  {
		    alan = Espec2->GetBinContent(i+1-100) - Euniverse1[k]->GetBinContent(i+1);
		  }
		if(j>101)
		  {
		    mary = Espec2->GetBinContent(j+1-100) - Euniverse1[k]->GetBinContent(j+1);
		  }
		mat->SetBinContent(i+1,j+1,alan*mary/float(iuniverse));

    // storing covariance matrix data
    matrix_data[i*(mat_len-1)+j] = alan*mary/float(iuniverse);
    covariance_matrix[i][j] += alan*mary/float(iuniverse);
    correlation_matrix[i][j] += alan*mary/float(iuniverse);

    
		if(k==4&&i==120&&j==120){
        //std::cout<<" i "<<i<<" j "<<j<<"alan "<<alan<<" mary"<<mary<<" j universe bin "<<Euniverse1[k]->GetBinContent(j+100)<<std::endl;
      }
	      }
	  }
      }

    double singular_adjust = 0.001;
    for(int w=0;w<mat_len;w++){
      //printf("diagonal %d value is %f \n",w,correlation_matrix[w][w]);
      correlation_matrix[w][w] = correlation_matrix[w][w] + singular_adjust;
      //printf("diagonal %d value is now %f \n",w,correlation_matrix[w][w]);
    }
   //adding covariance values into matrix object
    // covariance_matrix.SetMatrixArray(matrix_data);
    // correlation_matrix.SetMatrixArray(matrix_data);

    double_t det;
    correlation_matrix.Invert(&det);

    TH2D *cov_hist = new TH2D(covariance_matrix);
    TH2D *cor_hist = new TH2D(correlation_matrix);
    cov_hist->SetName("covariance1");
    cor_hist->SetName("correlation1");
    printf("The determinant is %f \n",det);

   //Output to file
   std::cout<<" made the covey matrix \n";
   TFile *newfile = new TFile("coveySIM.root","RECREATE");
   newfile->cd();
   mat->Write();
   Euniverse1[4]->Write();
   cov_hist->Write();
   cor_hist->Write();
   newfile->Close();
   
   mat->Delete();

   for (int k=0; k<100; k++){
     if (k != 4)
     {
       Euniverse1[k]->Delete();
     }
     
   }

   return Euniverse1[4];

}

/**
 * Function to create relative version of matrix.
 * M_sys/M_nosys = M_rel.
 **/
void makeRelative(TMatrixD &ResultMatrix,TMatrixD &SysMatrix,TMatrixD &OrigMatrix, TH1F *spectrum){
  for(int i = 0; i<200; ++i){
    for(int j = 0; j<200; ++j){
      if(OrigMatrix[i][j]!= 0.0){
      ResultMatrix[i][j] = SysMatrix[i][j] / spectrum->Integral();
      }else if (OrigMatrix[i][j]== 0.0)
      {
        ResultMatrix[i][j] = 0.0;//SysMatrix[i][j] / (OrigMatrix[i][j]+0.001);
      }
      
    }
  }
}

void applyStats(TMatrixD &Matrix){
  for(int i = 0; i<200; ++i){
    Matrix[i][i] = Matrix[i][i]*1;
  }
}

void applyScaling(TMatrixD &Matrix,TH1F *Spectrum){
  for(int i = 0; i<200; ++i){
    Matrix[i][i] = Matrix[i][i]*1;
  }
}

void doPlotMatrices(int save_option, TMatrixD Cov1, TMatrixD Cov2, const char* canvas_name, const char* title1, const char* title2, TH1F *spectrum){

    TH2D *cov_hist_plot = new TH2D("covariance1",title1,200,0.0,10.0,200,0.0,10.0);
    TH2D *cor_hist_plot = new TH2D("correlation1",title2,200,0.0,10.0,200,0.0,10.0);

    TH2D *cov_hist = new TH2D(Cov1);
    TH2D *inv_cov_hist = new TH2D(Cov2);
    
    // const char* canvas_name = "TestingRelative";


    auto myCanvas = new TCanvas(canvas_name,canvas_name);

    myCanvas->SetGrid();
    myCanvas->Divide(2,2);


    for(int i = 0; i<200; i++){
      for(int j = 0; j<200; j++){
        cov_hist_plot->SetBinContent(i+1,j+1,Cov1[i][j]);
        cor_hist_plot->SetBinContent(i+1,j+1,Cov2[i][j]);
      }
    }
    
    cov_hist_plot->GetXaxis()->SetTitle("E GeV");
    cor_hist_plot->GetXaxis()->SetTitle("E GeV");
    cov_hist_plot->GetYaxis()->SetTitle("E GeV");
    cor_hist_plot->GetYaxis()->SetTitle("E GeV");

    
    
    //Relabelling Spectra
  //   for (int i = 0; i < 10-1; ++i) {
  //   int bin{Euniverse1[4]->GetXaxis()->FindBin(i)};
  //   printf("bin number: %d \n",bin);
  //   std::string bin_label{std::to_string(i % 10)};
  //   Euniverse1[4]->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
  //   //bin = Euniverse1[4]->GetYaxis()->FindBin(i);
  //   //Euniverse1[4]->GetYaxis()->SetBinLabel(bin, bin_label.c_str());
  // }

  int universe_choice = 3;

  //Relabel axes
  for(int i=0; i<10; i++){
    int bin = spectrum->GetXaxis()->FindBin(i);
    std::string bin_label;

    if(i<=5){
      bin_label = std::to_string(i);
    }else if (i>5)
    {
      bin_label = std::to_string(i%5);
    }
    
    //std::string bin_label{std::to_string(i % 10)};
    //printf("bin number: %d \n",bin);
    //cout << bin_label << "\n";
    spectrum->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
    inv_cov_hist->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
    cov_hist_plot->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
  }


    //// Plotting Matrices
    myCanvas->cd(3);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    cov_hist_plot->Draw("COLZ");
    cov_hist_plot->LabelsOption("h","X");
    

    myCanvas->cd(4);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    cor_hist_plot->Draw("COLZ");
    cor_hist_plot->LabelsOption("h","X");

    
    
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)cov_hist_plot->FindObject("stats");
    TPaveStats *st2 = (TPaveStats*)cor_hist_plot->FindObject("stats");

    

    //stat box position
    st1->SetX2NDC(0.2); 
    st1->SetY2NDC(0.65); 
    st1->SetX1NDC(0.4); 
    st1->SetY1NDC(0.85); 

    //myCanvas->Update();

    st1->SetOptStat(111110110);

    st2->SetX2NDC(0.2); 
    st2->SetY2NDC(0.65); 
    st2->SetX1NDC(0.4); 
    st2->SetY1NDC(0.85); 

    st2->SetOptStat(111110110);

    myCanvas->Update();
  
   spectrum->LabelsOption("h","X");
   spectrum->SetTitle("Energy Spectrum");

   spectrum->GetXaxis()->SetTitle("E GeV");
   spectrum->SetName("Merged");   
   spectrum->SetStats(false);
  
    ///Plotting Spectra
    myCanvas->cd(1);
    gPad->SetRightMargin(0.15);
    spectrum->Draw("hist");
    //spectrum->LabelsOption("v","X");
  

    myCanvas->cd(2);
    gPad->SetRightMargin(0.15);
    spectrum->Draw("hist");
    //spectrum->LabelsOption("v","X");

    // TPaveStats *st3 = (TPaveStats*)Euniverse1[universe_choice]->FindObject("stats");
    // TPaveStats *st4 = (TPaveStats*)Euniverse1[universe_choice]->FindObject("stats");

    string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/JennyCode/Plotting/";
    string file_svg = ".svg";
    string file_pdf = ".pdf";


    string filestring1 = filepath + (string)canvas_name + file_svg;
    char filename1[filestring1.length() + 1];
    strcpy(filename1,filestring1.c_str());

    string filestring2 = filepath + (string)canvas_name + file_pdf;
    char filename2[filestring2.length() + 1];
    strcpy(filename2,filestring2.c_str());

    if(save_option == 1){
    myCanvas->SaveAs(filename1);
    myCanvas->SaveAs(filename2);
    }
}

void doPlotMatrices_single(int save_option, TMatrixD Cov1, const char* canvas_name, const char* title1, TH1F *spectrum){

    TH2D *cov_hist_plot = new TH2D("covariance1",title1,200,0.0,10.0,200,0.0,10.0);
    // TH2D *cor_hist_plot = new TH2D("correlation1",title2,200,0.0,10.0,200,0.0,10.0);

    TH2D *cov_hist = new TH2D(Cov1);
    
    // const char* canvas_name = "TestingRelative";


    auto myCanvas = new TCanvas(canvas_name,canvas_name,200,10,300,500);

    myCanvas->SetGrid();

    TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0,0);
    TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.3,0);
    pad1->Draw();
    pad2->Draw();
    //myCanvas->Divide(1,2);


    for(int i = 0; i<200; i++){
      for(int j = 0; j<200; j++){
        cov_hist_plot->SetBinContent(i+1,j+1,Cov1[i][j]);
        // cor_hist_plot->SetBinContent(i+1,j+1,Cov2[i][j]);
      }
    }
    
    cov_hist_plot->GetXaxis()->SetTitle("E GeV");
    // cor_hist_plot->GetXaxis()->SetTitle("E GeV");
    cov_hist_plot->GetYaxis()->SetTitle("E GeV");
    // cor_hist_plot->GetYaxis()->SetTitle("E GeV");

    
    
    //Relabelling Spectra
  //   for (int i = 0; i < 10-1; ++i) {
  //   int bin{Euniverse1[4]->GetXaxis()->FindBin(i)};
  //   printf("bin number: %d \n",bin);
  //   std::string bin_label{std::to_string(i % 10)};
  //   Euniverse1[4]->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
  //   //bin = Euniverse1[4]->GetYaxis()->FindBin(i);
  //   //Euniverse1[4]->GetYaxis()->SetBinLabel(bin, bin_label.c_str());
  // }

  int universe_choice = 3;

  //Relabel axes
  for(int i=0; i<10; i++){
    int bin = spectrum->GetXaxis()->FindBin(i);
    std::string bin_label;

    if(i<=5){
      bin_label = std::to_string(i);
    }else if (i>5)
    {
      bin_label = std::to_string(i%5);
    }
    
    //std::string bin_label{std::to_string(i % 10)};
    //printf("bin number: %d \n",bin);
    //cout << bin_label << "\n";
    spectrum->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
    // inv_cov_hist->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
    cov_hist_plot->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
  }


    //// Plotting Matrices
    pad1->cd();
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    cov_hist_plot->Draw("COLZ");
    cov_hist_plot->LabelsOption("h","X");
    

    // myCanvas->cd(4);
    // gPad->SetLogz();
    // gPad->SetRightMargin(0.15);
    // cor_hist_plot->Draw("COLZ");
    // cor_hist_plot->LabelsOption("h","X");

    
    
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)cov_hist_plot->FindObject("stats");
    // TPaveStats *st2 = (TPaveStats*)cor_hist_plot->FindObject("stats");

    

    //stat box position
    st1->SetX2NDC(0.2); 
    st1->SetY2NDC(0.65); 
    st1->SetX1NDC(0.4); 
    st1->SetY1NDC(0.85); 

    //myCanvas->Update();

    st1->SetOptStat(111110110);

    // st2->SetX2NDC(0.2); 
    // st2->SetY2NDC(0.65); 
    // st2->SetX1NDC(0.4); 
    // st2->SetY1NDC(0.85); 

    // st2->SetOptStat(111110110);

    myCanvas->Update();
  
   spectrum->LabelsOption("h","X");
   spectrum->SetTitle("Energy Spectrum");

   spectrum->GetXaxis()->SetTitle("E GeV");
   spectrum->SetName("Merged");   
   spectrum->SetStats(false);
  
    ///Plotting Spectra
    pad2->cd();
    gPad->SetRightMargin(0.15);
    spectrum->Draw("hist");
    //spectrum->LabelsOption("v","X");
  

    // myCanvas->cd(2);
    // gPad->SetRightMargin(0.15);
    // spectrum->Draw("hist");
    //spectrum->LabelsOption("v","X");

    // TPaveStats *st3 = (TPaveStats*)Euniverse1[universe_choice]->FindObject("stats");
    // TPaveStats *st4 = (TPaveStats*)Euniverse1[universe_choice]->FindObject("stats");

    string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/JennyCode/Plotting/Report/";
    string file_svg = ".svg";
    string file_pdf = ".pdf";


    string filestring1 = filepath + (string)canvas_name + file_svg;
    char filename1[filestring1.length() + 1];
    strcpy(filename1,filestring1.c_str());

    string filestring2 = filepath + (string)canvas_name + file_pdf;
    char filename2[filestring2.length() + 1];
    strcpy(filename2,filestring2.c_str());

    if(save_option == 1){
    myCanvas->SaveAs(filename1);
    myCanvas->SaveAs(filename2);
    }
}
int main(int argc, char* argv[])
{ //argv[0] is name of programme 

  TApplication *app = new TApplication("app", &argc, argv);

  TMatrixD covariance_matrix_noe(200,200);
  TMatrixD covariance_matrix_bothe(200,200);
  TMatrixD covariance_matrix_shifte(200,200);

  TMatrixD covariance_matrix(200,200);

  TMatrixD inv_covariance_matrix(200,200);
  TMatrixD inv_covariance_matrix2(200,200);
  TMatrixD inv_covariance_matrix3(200,200);

  TH1F *spectrum = createCovariance(covariance_matrix_bothe,inv_covariance_matrix,1);
  TH1F *spectrum2 = createCovariance(covariance_matrix_noe,inv_covariance_matrix2,-1);
  TH1F *spectrum3 = createCovariance(covariance_matrix_shifte,inv_covariance_matrix3,0);

  TMatrixD covariance_matrix_rel(200,200);
  TMatrixD inv_covariance_matrix_rel(200,200);
  TMatrixD covariance_matrix_rel_shift(200,200);
  TMatrixD inv_covariance_matrix_rel_shift(200,200);

  makeRelative(covariance_matrix_rel,covariance_matrix_bothe,inv_covariance_matrix2,spectrum);
  makeRelative(inv_covariance_matrix_rel,inv_covariance_matrix,inv_covariance_matrix2,spectrum);
  applyStats(covariance_matrix_rel);
  applyScaling(covariance_matrix_rel,spectrum);

  makeRelative(covariance_matrix_rel_shift,covariance_matrix_shifte,inv_covariance_matrix3,spectrum3);
  makeRelative(inv_covariance_matrix_rel_shift,inv_covariance_matrix3,inv_covariance_matrix3,spectrum3);
  applyStats(covariance_matrix_rel_shift);
  applyScaling(covariance_matrix_rel_shift,spectrum3);


   if (userConfirm() == 1){
//Some Plotting of matrices and other stuff
    std::cout<<"Starting plotting \n ";

    // doPlotMatrices(1,covariance_matrix_bothe,inv_covariance_matrix,"Normal Both Errors","Covariance Both Errors","Inverse Covariance Both Errors",spectrum);
    // doPlotMatrices(1,covariance_matrix_rel,inv_covariance_matrix_rel,"Relative Both Errors no scaling","Covariance Both Errors Relative","Inverse Covariance Both Errors Relative",spectrum);
    
    // doPlotMatrices(1,covariance_matrix_shifte,inv_covariance_matrix3,"Normal Shift Error","Covariance Shift Error","Inverse Covariance Shift Error",spectrum3);
    // doPlotMatrices(1,covariance_matrix_rel_shift,inv_covariance_matrix_rel_shift,"Relative Shift Error no scaling","Covariance Shift Error Relative","Inverse Covariance Shift Error Relative",spectrum);

    // doPlotMatrices(1,covariance_matrix_bothe,covariance_matrix_noe,"Normal Both and None Errors","Covariance Both Errors","Covariance No Errors",spectrum);
    // doPlotMatrices(1,covariance_matrix_rel,covariance_matrix_noe,"Relative Both Errors","Covariance Both Errors Relative","Covariance No Errors",spectrum);
  

    doPlotMatrices_single(1,covariance_matrix_noe,"NoErrors","Covariance with No Errors",spectrum2);
    doPlotMatrices_single(1,inv_covariance_matrix,"InverseNoErrors","Inverse with No Errors",spectrum);
    doPlotMatrices_single(1,covariance_matrix_shifte,"ShiftError","Covariance with Shift Error",spectrum3);
    doPlotMatrices_single(1,covariance_matrix_bothe,"BothErrors","Covariance with Shift and Calibration Errors",spectrum);
    doPlotMatrices_single(1,covariance_matrix_rel,"RelativeBoth Errors","Relative Covariance with Noth Errors",spectrum);
    doPlotMatrices_single(1,inv_covariance_matrix_rel,"InverseRelative Both Errors","Inverse Relative Covariance with both Errors",spectrum);

    app->Run();
	  
}
}


