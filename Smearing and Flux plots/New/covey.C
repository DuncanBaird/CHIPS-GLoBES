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

void plotSmear(string filename){

  ifstream file;
  string title = filename;
  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/New/Data/Smear/";
  filepath.append(filename.append(".dat"));
file.open(filepath.c_str());
printf(filepath.c_str());
if (!file.is_open()) {printf("Error file could not open\n");}

const int number_rows = 86;
const int number_columns = 88;
double my_array [number_rows][number_columns]{};
double value;
for (int i{}; i != number_rows; ++i) {
    for (int j{}; j != number_columns; ++j) {
      file >> value;
      if(value==85){
        my_array[i][j] = value;
      }else{
        my_array[i][j] = value;
      }
    }
}



  double bins_array[number_rows+1] = {0,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4};

  for (int edge=0;edge<=number_rows;edge++){
    bins_array[edge]=bins_array[edge]+bins_array[edge-1];
    //printf("Bin edges: %f\n",bins_array[edge]);
  }
  char dist_title[30];
  strcpy(dist_title,title.c_str());

  TH2D *hist = new TH2D(dist_title,"2D Histrogram of Smearing Matrix",number_rows,bins_array,number_columns-2,bins_array);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");


  for (int i=1;i<number_rows;i++){
     for(int j=3;j<number_columns;j++){
         
         //hist->Fill(i,j-2,my_array[i][j]);
         hist->SetBinContent(i,j-2,my_array[i][j]);
         //printf("Bin array: %f\n",my_array[i][j]);
     }
  }
  hist->Draw("COLZ");

  gPad->Update();
  TPaveStats *st = (TPaveStats*)hist->FindObject("stats");

  //stat box position
  st->SetX2NDC(0.6); 
  st->SetY2NDC(0.15); 
  st->SetX1NDC(0.8); 
  st->SetY1NDC(0.35); 

  //st->SetOptStat(111110110);
}

void plotCross(){

}


int main(int argc, char* argv[])
{ //argv[0] is name of programme 

  TApplication *app = new TApplication("app", &argc, argv);

  //  //Output to file
  //  std::cout<<" made the covey matrix \n";
  //  TFile *newfile = new TFile("coveySIM.root","RECREATE");
  //  newfile->cd();
  //  mat->Write();
  //  Euniverse1[4]->Write();
  //  cov_hist->Write();
  //  cor_hist->Write();
  //  newfile->Close();

   if (userConfirm() == 1){
//Some Plotting of matrices and other stuff
    
    const char* canvas_name = "Smearing Matrices";


    auto myCanvas = new TCanvas(canvas_name,canvas_name);

    myCanvas->SetGrid();
    myCanvas->Divide(3,3);

    string smear_filenames[7] = {"smear_anu_mucc_sk2","smear_anu_nqe_sk2","smear_anu_qe_sk2","smear_nc_sk2","smear_nu_mucc_sk2","smear_nu_nqe_sk2","smear_nu_qe_sk2"};

    for(int i=0;i<7;i++){
      myCanvas->cd(i+1);
      plotSmear(smear_filenames[i]);
    }

    // for(int i = 0; i<200; i++){
    //   for(int j = 0; j<200; j++){
    //     cov_hist_plot->SetBinContent(i+1,j+1,covariance_matrix[i][j]);
    //     cor_hist_plot->SetBinContent(i+1,j+1,correlation_matrix[i][j]);
    //   }
    // }
    
    // cov_hist_plot->GetXaxis()->SetTitle("E GeV");
    // cor_hist_plot->GetXaxis()->SetTitle("E GeV");
    // cov_hist_plot->GetYaxis()->SetTitle("E GeV");
    // cor_hist_plot->GetYaxis()->SetTitle("E GeV");

    // //// Plotting Matrices
    // myCanvas->cd(3);
    // gPad->SetLogz();
    // gPad->SetRightMargin(0.15);
    // cov_hist_plot->Draw("COLZ");
    

    // myCanvas->cd(4);
    // gPad->SetLogz();
    // gPad->SetRightMargin(0.15);
    // cor_hist_plot->Draw("COLZ");

    
    
    // gPad->Update();
    // TPaveStats *st1 = (TPaveStats*)cov_hist_plot->FindObject("stats");
    // TPaveStats *st2 = (TPaveStats*)cor_hist_plot->FindObject("stats");

    

    // //stat box position
    // st1->SetX2NDC(0.2); 
    // st1->SetY2NDC(0.65); 
    // st1->SetX1NDC(0.4); 
    // st1->SetY1NDC(0.85); 

    // //myCanvas->Update();

    // st1->SetOptStat(111110110);

    // st2->SetX2NDC(0.2); 
    // st2->SetY2NDC(0.65); 
    // st2->SetX1NDC(0.4); 
    // st2->SetY1NDC(0.85); 

    // st2->SetOptStat(111110110);

    // myCanvas->Update();
    


 
    // string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux Plots/New/";
    // string file_svg = ".svg";
    // string file_pdf = ".pdf";


    // string filestring1 = filepath + (string)canvas_name + file_svg;
    // char filename1[filestring1.length() + 1];
    // strcpy(filename1,filestring1.c_str());

    // string filestring2 = filepath + (string)canvas_name + file_pdf;
    // char filename2[filestring2.length() + 1];
    // strcpy(filename2,filestring2.c_str());

    // if(userConfirm()==1){
    // myCanvas->SaveAs(filename1);
    // myCanvas->SaveAs(filename2);
    }
   
   
    

    app->Run();
	  
}


