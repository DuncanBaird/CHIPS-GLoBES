#include <fstream>

TH1* chi2Hist = NULL;

void printer3(int w){
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  myCanvas->DivideSquare(18);

  const int number_histograms = 18;

  char buf[120];
  char str[120];

  for (int k = 0; k<number_histograms;k++){
    myCanvas->cd(k+1);
    gPad->Update();

  sprintf(buf,"/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/%dtime2histsim_channel%d.dat",w,k);
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
  sprintf(buf2,"Time: %d, Interaction channel %d",w,k);
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

  char buf4[120];
  char str4[120];

  sprintf(buf4,"/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/CHIPShist-time-%d.pdf",w);
  strcpy(str4,buf4);

  myCanvas->SaveAs();

  }

void manygraphinghistmanytime(){

int tSteps = 30;
  
for (int j=0; j <= tSteps; j++){

  printer3(j);

}
 

}