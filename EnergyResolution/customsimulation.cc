/*
 * Duncan Baird; Simulation for doing several plots.
 * Example: Create sensitivity plot for reactor experiments with near
 * and far detectors.
 * Compile with ``make customsimulation''
 */

#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include "initialsimulation.h"
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TF1.h>
#include <TError.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>
using namespace std;

/***************************************************************************
 *                   P L O T T I N G   F U N C T I O N                     *
 ***************************************************************************/
const int tSteps = 1;
void doPlotROOT(double plot_data_a[][tSteps],double plot_data_b[][tSteps],int option,const char* canvas_name){

for (int g = 0;g<tSteps; g++){
    printf("Data a is: index: %f and value: %f \n",plot_data_a[0][g],plot_data_a[1][g]);
    printf("Data b is: index: %f and value: %f \n",plot_data_b[0][g],plot_data_b[1][g]);
  }
  plot_data_a[1][0] = plot_data_a[1][1];
  plot_data_b[1][0] = plot_data_b[1][1];

  auto myCanvas = new TCanvas(canvas_name,canvas_name);
  myCanvas->SetGrid();
  auto spa = new TGraph(tSteps,plot_data_a[0],plot_data_a[1]);
  auto spb = new TGraph(tSteps,plot_data_b[0],plot_data_b[1]);

  spa->SetMarkerStyle(4);
  spa->SetMarkerColor(4);
  spb->SetMarkerStyle(3);
  spb->SetMarkerColor(3);
  
  auto mg = new TMultiGraph();

  mg->SetTitle("AB: Log Plot of sensitivity curve systematics vs stats only.");
  mg->Add(spa);
  mg->Add(spb);
  mg->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  mg->GetYaxis()->SetTitle("sin(2*theta_13)^2 sensitiivity");
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);

  gPad->SetLogx();
  gPad->Update();
  mg->Draw("APL");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  legend->AddEntry(spa,"series A: With Systematics","lp");
  legend->AddEntry(spb,"series B: Statistics only","lp");
  legend->Draw();

  myCanvas->Update();
  
  myCanvas->Modified();
  myCanvas->Update();

  if (option == 1){

  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewSystematicsROOTv2/Plotting/Output/";
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

void userConfirm(){
  int kz;
  cout << "Please enter an integer to proceed: ";
  cin >> kz;
  cout << "The value you entered is " << kz;
  cout << " and its double is " << kz*2 << ".\n";
  //return 0;
}

///READ CSV DATA
std::vector<std::pair<std::string, std::vector<float>>> read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<float>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    float val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<float> {}});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
        }
    }

    // Close file
    myFile.close();

    return result;
}  

void doResolutionHist(THStack* hs,float mean,float sigma,int j){
  TRandom3 rndgen;
  char no[5];
  sprintf (no, "no%d",j);
  TH1* h1 = new TH1D(no, no, 1000, 0.0, 5000.0);

  TF1  *f1 = new TF1(no,"gaus",0,5000);
  f1->FixParameter(1,1.0);
  f1->FixParameter(2,0.5);

  


  for (int i=0;i < 1E3; ++i){

    h1->Fill(rndgen.Gaus(mean,sigma));
  }

  h1->Fit(f1);

  // for (int i = 0;i<h1->GetNbinsX();i++){

  //   h1->SetBinContent(i,0.0);
  // }

  h1->SetFillColorAlpha(0,0);
  h1->SetMarkerColorAlpha(0,0);

  h1->Draw("C");
  f1->Draw("SAME");

  hs->Add(h1);

  
  //delete h1;
  //delete f1;

}

float getSigma(float r_val, float mean){
  return (r_val * mean)/(2*sqrt(2*log(2)));
}

void doPlotResolutions(std::vector<std::pair<std::string, std::vector<float>>> plot_data){
  const char* canvas_name = "Energy Resolution Peak";
  const char* canvas_name2 = "Energy Resolution Peak Stack";

  THStack* hs = new THStack("hs",canvas_name2);


  
  auto myCanvas = new TCanvas(canvas_name,canvas_name);
  myCanvas->SetGrid();

  

  gPad->SetLogy();

  std::vector<float> energy_val = plot_data[0].second;
  std::vector<float> e_val = plot_data[1].second;
  std::vector<float> mu_val = plot_data[2].second;
  int count = 0;
  for(auto i = energy_val.begin(); i!= energy_val.end();++i ){

    doResolutionHist(hs,energy_val[count],getSigma(mu_val[count],energy_val[count]),count);
    printf("Number is %d: \n",count);
    ++count;
    
  }
  
  hs->Draw("nostack");
}

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
  //TApplication app{"customsimulation", &argc, argv};
  TApplication *app = new TApplication("app", &argc, argv);
  //gErrorIgnoreLevel = kFatal;

  TFile *file = new TFile("test.root", "RECREATE");
  file->Write();
  file->Close();
  delete file;

  std::vector<std::pair<std::string, std::vector<float>>> data = read_csv("resolution.csv");
  
  // std::pair<std::string, std::vector<double>> test = (data[0]);
  // std::vector<double> test2 = test.second;
  // double test3 = test2[0];
  // for (int i =0;i<7;i++){
  //   printf("extracted value is %f \n",test2[i]);
  // }
/***************************************************************************
 *                                P L O T T I N G                          *
 ***************************************************************************/

  doPlotResolutions(data);

  app->Run();

  

  return 0;
}