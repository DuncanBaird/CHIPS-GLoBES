#include <fstream>

TH1* chi2Hist = NULL;

void graphinghist(){

  ifstream file;
  file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/histsim.dat");
  if (!file.is_open()) {printf("Error\n");}else{

  const int len = 86;
  double x[85]{};
  double w[85]{};
  double value;
  double weight;
  for (int i = 0; i != 86; ++i) {
      file >> value >> weight;
      x[i] = value;
      w[i] = weight;
      printf("bin %f",x[i]);
  }
file.close();


  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  chi2Hist = new TH1D("Interaction channel 0","Histogram of Neutrino Energies",86,0,86); //last numbers: number of bins, lower bound of bins, upper bound of bins
  chi2Hist->GetXaxis()->SetTitle("Energy GeV");
  chi2Hist->GetYaxis()->SetTitle("Frequency");
  for(int i = 0; i != 86; ++i) {
  chi2Hist->Fill(x[i],w[i]);}

  chi2Hist->Draw("hist");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)chi2Hist->FindObject("stats");
  myCanvas->Update();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/CHIPShist-First.pdf");

  }
}