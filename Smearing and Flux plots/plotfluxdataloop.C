#include <fstream>
void plotfluxdataloop(){



ifstream file;
file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/flux-numi-me-7mrad-minus.dat");
//file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/flux-numi-me-7mrad-plus.dat");
if (!file.is_open()) {printf("Error\n");}

const int number_entries = 500;

double energy [number_entries]{};
double flux_mu [number_entries]{};
double flux_e [number_entries]{};
double flux_tau [number_entries]{};
double flux_mu_anti [number_entries]{};
double flux_e_anti [number_entries]{};
double flux_tau_anti [number_entries]{};



for (int i{}; i != number_entries; ++i) {
      file >> energy[i] >> flux_mu[i] >> flux_e[i] >> flux_tau[i] >> flux_mu_anti[i] >> flux_e_anti[i] >> flux_tau_anti[i];
}



TH1D *hist_mu = new TH1D("mu","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_mu->GetXaxis()->SetTitle("Energy GeV");
  hist_mu->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_mu->Fill(energy[i],flux_mu[i]);
  }
TCanvas *myCanvas_mu = new TCanvas();
myCanvas_mu->SetGrid();
hist_mu->Draw("HIST");
myCanvas_mu->Update();
myCanvas_mu->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/mufluxminus.svg");
myCanvas_mu->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/mufluxminus.pdf");


TH1D *hist_e = new TH1D("e","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_e->GetXaxis()->SetTitle("Energy GeV");
  hist_e->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_e->Fill(energy[i],flux_e[i]);
  }
TCanvas *myCanvas_e = new TCanvas();
myCanvas_e->SetGrid();
hist_e->Draw("HIST");
myCanvas_e->Update();
myCanvas_e->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/efluxminus.svg");
myCanvas_e->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/efluxminus.pdf");

TH1D *hist_tau = new TH1D("tau","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_tau->GetXaxis()->SetTitle("Energy GeV");
  hist_tau->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_tau->Fill(energy[i],flux_tau[i]);
  }
TCanvas *myCanvas_tau = new TCanvas();
myCanvas_tau->SetGrid();
hist_tau->Draw("HIST");
myCanvas_tau->Update();
myCanvas_tau->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/taufluxminus.svg");
myCanvas_tau->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/taufluxminus.pdf");


TH1D *hist_mu_anti = new TH1D("anti mu","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_mu_anti->GetXaxis()->SetTitle("Energy GeV");
  hist_mu_anti->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_mu_anti->Fill(energy[i],flux_mu_anti[i]);
  }
TCanvas *myCanvas_mu_anti = new TCanvas();
myCanvas_mu_anti->SetGrid();
hist_mu_anti->Draw("HIST");
myCanvas_mu_anti->Update();
myCanvas_mu_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antimufluxminus.svg");
myCanvas_mu_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antimufluxminus.pdf");


  TH1D *hist_e_anti = new TH1D("anti e","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_e_anti->GetXaxis()->SetTitle("Energy GeV");
  hist_e_anti->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_e_anti->Fill(energy[i],flux_e_anti[i]);
  }
TCanvas *myCanvas_e_anti = new TCanvas();
myCanvas_e_anti->SetGrid();
hist_e_anti->Draw("HIST");
myCanvas_e_anti->Update();
myCanvas_e_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antiefluxminus.svg");
myCanvas_e_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antiefluxminus.pdf");


  TH1D *hist_tau_anti = new TH1D("anti tau","Histogram of Flux Minus",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist_tau_anti->GetXaxis()->SetTitle("Energy GeV");
  hist_tau_anti->GetYaxis()->SetTitle("Flux");
  for (int i=0;i<number_entries;i++){
         
         hist_tau_anti->Fill(energy[i],flux_tau_anti[i]);
  }
TCanvas *myCanvas_tau_anti = new TCanvas();
myCanvas_tau_anti->SetGrid();
hist_tau_anti->Draw("HIST");
myCanvas_tau_anti->Update();
myCanvas_tau_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antitaufluxminus.svg");
myCanvas_tau_anti->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/antitaufluxminus.pdf");


  gPad->Update();
  //TPaveStats *st = (TPaveStats*)hist_mu->FindObject("stats");

  //stat box position
  //st->SetX2NDC(0.6); 
  //st->SetY2NDC(0.15); 
  //st->SetX1NDC(0.8); 
  //st->SetY1NDC(0.35); 

  //st->SetOptStat(111110110);

  //myCanvas->Update();
  //myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/GraphSmearMatrixCorrectedNewAxis.svg");
}