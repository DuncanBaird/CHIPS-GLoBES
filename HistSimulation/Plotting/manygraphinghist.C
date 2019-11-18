#include <fstream>

TH1* chi2Hist = NULL;

void manygraphinghist(){

  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  myCanvas->DivideSquare(18);

  const int number_histograms = 18;

  char buf[120];
  char str[120];

  for (int k = 0; k<number_histograms;k++){
    myCanvas->cd(k+1);
    gPad->Update();

  sprintf(buf,"/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/histsim_channel%d.dat",k);
  strcpy(str,buf);
  
  ifstream file;
  file.open(str);
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
  }
file.close();

  char buf2[25];
  char str2[25];
  sprintf(buf2,"Interaction channel %d",k);
  strcpy(str2,buf2);

  chi2Hist = new TH1D(str2,str2,86,0,86); //last numbers: number of bins, lower bound of bins, upper bound of bins
  chi2Hist->GetXaxis()->SetTitle("Energy GeV");
  chi2Hist->GetYaxis()->SetTitle("Frequency");
  for(int i = 0; i != 86; ++i) {
  chi2Hist->Fill(x[i],w[i]);}

  chi2Hist->Draw("hist");
  
  TPaveStats *st = (TPaveStats*)chi2Hist->FindObject("stats");


  }
  
  }

  myCanvas->Update();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/CHIPShist-many.pdf");


}