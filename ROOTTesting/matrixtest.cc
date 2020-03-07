#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

#include <TMatrixD.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TFile.h>

int setIdentitiy(TMatrixD &dummy1,int len){
  for(int i = 0;i<len;i++){
    for(int j = 0;j<len;j++){
        dummy1[i][j] = 0;
      }
    }
  for(int i = 0;i<len;i++){

    dummy1[i][i]=1;

  }
  return 0;
}

int setValue(TMatrixD &dummy1,int len_r,int len_c, int value){
  for(int i = 0;i<len_r;i++){
    for(int j = 0;j<len_c;j++){
        dummy1[i][j] = value;
      }
    }

  return 0;
}

void generate_covariance(TMatrixD covariance_matrix){
  
 const int mat_len = 200;
  //TMatrixD covariance_matrix(200,200);
  TMatrixD correlation_matrix(200,200);
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

   TFile *ff3 = new TFile("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/ROOTTesting/TrueE.root");
   TH1F *Espec1 = (TH1F*)ff3->Get("trueE_0to5");
   Espec1->SetName("Espec1");


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
   for (int i=0;i<10000;i++)
     {
       Espec2->Fill(g->Gaus(2.0,0.3));
     }

   //make the universes based on the real event spectrum
   std::cout<<" making the universes "<<std::endl;
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
	     E1 = E1+0.05*E1;
	     E2 = E2+0.05*E2;
	     //float rfake = g->Gaus(E,res[0]/sqrt(E));
	     //float Eres = E;//+rfake;
	     //fake11->Fill(Eres,ifake);

	     Euniverse1[k]->Fill(E1,i1);
	     Euniverse1[k]->Fill(E2+5,i2); //E2+5.0

	  } ///finished spectrum generation for this k universe
	std::cout<<" finished this universe "<<k<<std::endl;
	for (int i=0; i<200; i++)
	  {
	    for (int j=0; j<200; j++)
	      {
		double alan, mary;
		if(i<100)
		  {
		    alan = Espec1->GetBinContent(i+1) - Euniverse1[k]->GetBinContent(i);
		  }
		if(j<100)
		  {
		    mary = Espec1->GetBinContent(j+1) - Euniverse1[k]->GetBinContent(j);
		  }
		if(i>101)
		  {
		    alan = Espec2->GetBinContent(i+1) - Euniverse1[k]->GetBinContent(i);
		  }
		if(j>101)
		  {
		    mary = Espec2->GetBinContent(j+1) - Euniverse1[k]->GetBinContent(j);
		  }
		mat->SetBinContent(i+1,j+1,alan*mary/float(iuniverse));

    // storing covariance matrix data
    matrix_data[i*(mat_len-1)+j] = alan*mary/float(iuniverse);
    covariance_matrix[i][j] += alan*mary/float(iuniverse);
    correlation_matrix[i][j] += alan*mary/float(iuniverse);

    
		if(k==4&&i==120&&j==120)std::cout<<" i "<<i<<" j "<<j<<"alan "<<alan<<" mary"<<mary<<" j universe bin "<<Euniverse1[k]->GetBinContent(j+100)<<std::endl;
	      }
	  }
      }

    double singular_adjust = 0.001;
    for(int w=0;w<mat_len;w++){
      printf("diagonal %d value is %f \n",w,correlation_matrix[w][w]);
      correlation_matrix[w][w] = correlation_matrix[w][w] + singular_adjust;
      printf("diagonal %d value is now %f \n",w,correlation_matrix[w][w]);
    }
   //adding covariance values into matrix object
    // covariance_matrix.SetMatrixArray(matrix_data);
    // correlation_matrix.SetMatrixArray(matrix_data);

    // double_t det;
    // correlation_matrix.Invert(&det);

}


void matrixtest(){

  const int len = 4;

  TMatrixD test1(len,len);
  TMatrixD test2(len,len);

  TMatrixD vertical(len,1);
  TMatrixD vertical2(len,1);
  TMatrixD horizontal(1,len);

  // for(int i = 0;i<len;i++){
  //   for(int j = 0;j<len;j++){
  //       test1[i][j] = 0;
  //     }
  //   }
  
  // for(int i = 0;i<len;i++){

  // test1[i][i]=1;

  // }

  setIdentitiy(test1,len);

  setValue(test2,len,len,4);

  setValue(horizontal,1,len,4);
  setValue(vertical,len,1,4);
  

  // test1.Print();
  // test2.Print();

  // TMatrixD test3 = test2 * test2;

  // test3.Print();

  // vertical.Print();
  // horizontal.Print();

  vertical2.T().Print();

  test2[0][0] = 12.0;
  test2.Print();


  double det1;

  test2.Invert(&det1);
  test2.Print();

  
  TMatrixD decomp = horizontal * test2 * vertical;

  decomp.Print();

  double x = decomp[0][0];

  cout << "chi is " << x << "\n";

  TMatrixD delta1(200,1);
  TMatrixD delta2(200,1);

  TMatrixD covariance(200,200);
  generate_covariance(covariance);
  // setValue(covariance,200,200,0.1);
  // covariance[0][0] = 100.0;

  for(int i=0;i<200;i++){
    delta1[i][0] = 5E10;
    delta2[i][0] = 5E10;
  }
  delta1.T();
  (delta1*covariance.Invert(&det1)*delta2).Print();
}