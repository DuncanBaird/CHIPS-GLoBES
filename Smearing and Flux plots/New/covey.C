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
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>



int userConfirm(){
  int kz;
  cout << "Please enter binary confirmation for plotting: ";
  cin >> kz;
  cout << "The value you entered is " << kz << ".\n";
  return kz;
}

void plotSmear(string filename, string title){

  string title_prefix = "Smearing Matrix of ";
  title_prefix.append(title);

  ifstream file;
  //string title = filename;
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

  TH2D *hist = new TH2D(dist_title,title_prefix.c_str(),number_rows,bins_array,number_columns-2,bins_array);//"distribution name","Name of plot",86,0,86,86,0,86);
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
  //TPaveStats *st = (TPaveStats*)hist->FindObject("stats");

  // //stat box position
  // st->SetX2NDC(0.6); 
  // st->SetY2NDC(0.15); 
  // st->SetX1NDC(0.8); 
  // st->SetY1NDC(0.35); 

  //st->SetOptStat(111110110);
  
}

void plotCross(string filename, string title){

  //wc_XCC.dat
  
  auto graph = new TMultiGraph();

  
  graph->GetXaxis()->SetTitle("E_{#nu} GeV");
  graph->GetYaxis()->SetTitle("#frac{#sigma_{L}}{E_{#nu_{L}}} #frac{10^{-14} b}{GeV}");

  // string title = filename;
  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/New/Data/Cross/";
  filepath.append(filename.append(".dat"));

  graph->SetTitle(title.c_str());

  cout << filepath;

  const int file_length = 1001;
  double e[file_length];
  double loge[file_length];
  double sig_mu[file_length];
  double sig_e[file_length];
  double sig_tau[file_length];
  double sig_mu_a[file_length];
  double sig_e_a[file_length];
  double sig_tau_a[file_length];

  ifstream infile;   

  infile.open(filepath.c_str());// file containing numbers in 3 columns 
     if (!infile.is_open()) {printf("Error file could not open\n");}
       
       int counter = 0;
       while(
         infile >> loge[counter] >> sig_mu[counter] >> sig_e[counter] >> sig_tau[counter] >> sig_mu_a[counter] >> sig_e_a[counter] >> sig_tau_a[counter])
         {
            //Looped area
            ++counter;
       } 

         printf("Values for 1000 are %f,%f,%f,%f\n",loge[1000],sig_mu[1000],sig_e[1000]);
        for(int i=0;i<file_length;++i){
          e[i] = pow(10.0,loge[i]);
          sig_mu[i] = sig_mu[i]/e[i];
          sig_e[i] = sig_e[i]/e[i];
          sig_tau[i] = sig_tau[i]/e[i];
          sig_mu_a[i] = sig_mu_a[i]/e[i];
          sig_e_a[i] = sig_e_a[i]/e[i];
          sig_tau_a[i] = sig_tau_a[i]/e[i];
          //printf("Sig mu value at point %d,energy %f, is %f \n",i,e[i],sig_mu[i]);
        }
        auto spa = new TGraph(file_length,e,sig_mu);
        auto spb = new TGraph(file_length,e,sig_e);
        auto spc = new TGraph(file_length,e,sig_tau);
        auto spd = new TGraph(file_length,e,sig_mu_a);
        auto spe = new TGraph(file_length,e,sig_e_a);
        auto spf = new TGraph(file_length,e,sig_tau_a);
        
        // spa->SetMarkerStyle(2);
        // spb->SetMarkerStyle(4);
        // spc->SetMarkerStyle(26);
        // spd->SetMarkerStyle(37);
        // spe->SetMarkerStyle(40);
        // spf->SetMarkerStyle(42);
        
        // spa->SetMarkerColor(2);
        // spb->SetMarkerColor(3);
        // spc->SetMarkerColor(4);
        // spd->SetMarkerColor(1);
        // spe->SetMarkerColor(6);
        // spf->SetMarkerColor(7);

        int marker_style = 5;
        int marker_style2 = 4;

        spa->SetMarkerStyle(marker_style);
        spb->SetMarkerStyle(marker_style2);
        spc->SetMarkerStyle(marker_style);
        spd->SetMarkerStyle(marker_style2);
        spe->SetMarkerStyle(marker_style);
        spf->SetMarkerStyle(marker_style);

        double marker_scale = 0.1;

        spa->SetMarkerSize(6*marker_scale);
        spb->SetMarkerSize(2*marker_scale);
        spc->SetMarkerSize(3*marker_scale);
        spd->SetMarkerSize(4*marker_scale);
        spe->SetMarkerSize(2*marker_scale);
        spf->SetMarkerSize(2*marker_scale);
        
        spa->SetMarkerColor(2);
        spb->SetMarkerColor(3);
        spc->SetMarkerColor(4);
        spd->SetMarkerColor(1);
        spe->SetMarkerColor(6);
        spf->SetMarkerColor(7);

        graph->Add(spa);
        graph->Add(spb);
        graph->Add(spc);
        graph->Add(spd);
        graph->Add(spe);
        graph->Add(spf);

        

        graph->Draw("APL");

        TLegend* legend = new TLegend();
        legend->SetHeader("Neutrino Flavour: L");
        legend->AddEntry(spa,"#mu","p");
        legend->AddEntry(spb,"e","p");
        legend->AddEntry(spc,"#tau","p");
        legend->AddEntry(spd,"#bar{#mu}","p");
        legend->AddEntry(spe,"#bar{e}","p");
        legend->AddEntry(spf,"#bar{#tau}","p");
        legend->Draw("SAME");

        gPad->SetLogx();
        gPad->SetLogy();
}

void doPlotSmear(int option){
  
  const char* canvas_name = "Smearing Matrices";


    auto myCanvas = new TCanvas(canvas_name,canvas_name);

    myCanvas->SetGrid();
    myCanvas->Divide(3,3);

    string smear_filenames[7] = {"smear_anu_mucc_sk2","smear_anu_nqe_sk2","smear_anu_qe_sk2","smear_nc_sk2","smear_nu_mucc_sk2","smear_nu_nqe_sk2","smear_nu_qe_sk2"};
    string smear_titles[7] = {"#bar{#nu}_{#mu} CC","#bar{#nu}_{e} NQE","#bar{#nu}_{e} QE","#nu_{#mu} NC","#nu_{#mu} CC","#nu_{e} CC NQE","#nu_{e} CC QE"};


    for(int i=0;i<7;i++){
      myCanvas->cd(i+1);
      plotSmear(smear_filenames[i],smear_titles[i]);
    }
  
  
  
  if (option == 1){

  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/New/";
  string file_svg = ".svg";
  string file_pdf = ".pdf";

  string filestring1 = filepath + (string)canvas_name + file_svg;
  char filename1[filestring1.length() + 1];
  strcpy(filename1,filestring1.c_str());

  string filestring2 = filepath + (string)canvas_name + file_pdf;
  char filename2[filestring2.length() + 1];
  strcpy(filename2,filestring2.c_str());

  myCanvas->SaveAs(filename1);
  myCanvas->SaveAs(filename2);
  }
}

void doPlotCross(int option){

   const char* canvas_name = "Cross Sections";


    auto myCanvas = new TCanvas(canvas_name,canvas_name);

    myCanvas->SetGrid();
    myCanvas->Divide(2,2);

  string cross_filenames[4] = {"wc_XCC","wc_XCCNonQE","wc_XNC","wc_XQE"};
  string cross_titles[4] = {"Charged Current","Charged Current \n Non-Quasielastic","Neutral Current","Charged Current Quasielastic"};

  
for(int i=0;i<4;i++){
      myCanvas->cd(i+1);
      gPad->SetLeftMargin(0.15);
      plotCross(cross_filenames[i],cross_titles[i]);
    }
  
  
  
  if (option == 1){

  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/New/";
  string file_svg = ".svg";
  string file_pdf = ".pdf";

  string filestring1 = filepath + (string)canvas_name + file_svg;
  char filename1[filestring1.length() + 1];
  strcpy(filename1,filestring1.c_str());

  string filestring2 = filepath + (string)canvas_name + file_pdf;
  char filename2[filestring2.length() + 1];
  strcpy(filename2,filestring2.c_str());

  myCanvas->SaveAs(filename1);
  myCanvas->SaveAs(filename2);
  }
  


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
  /* Sets width of marker lines*/
  gStyle->SetLineScalePS(0.25);
  
  gStyle->SetOptStat(0);
  


   if (userConfirm() == 1){
//Some Plotting of matrices and other stuff
    
    doPlotSmear(1);
    
    }
    userConfirm();

    doPlotCross(1);
   
   
    

    app->Run();
	  
}


