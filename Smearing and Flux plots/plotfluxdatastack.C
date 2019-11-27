#include <fstream>
void plotfluxdatastack(){

   TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

ifstream file;
file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/flux-numi-me-7mrad-plus.dat");
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

// Create Histogram Stack
THStack *hs = new THStack("hs","test stacked histograms");


TH1D *hist_mu = new TH1D("mu","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_mu->Fill(energy[i],flux_mu[i]);
  }

TH1D *hist_e = new TH1D("e","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_e->Fill(energy[i],flux_e[i]);
  }

TH1D *hist_tau = new TH1D("tau","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_tau->Fill(energy[i],flux_tau[i]);
  }

TH1D *hist_mu_anti = new TH1D("anti mu","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_mu_anti->Fill(energy[i],flux_mu_anti[i]);
  }

  TH1D *hist_e_anti = new TH1D("anti e","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_e_anti->Fill(energy[i],flux_e_anti[i]);
  }

  TH1D *hist_tau_anti = new TH1D("anti tau","Histogram of Flux",100,0,30);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  //hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  //hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");
  for (int i=0;i<number_entries;i++){
         
         hist_tau_anti->Fill(energy[i],flux_tau_anti[i]);
  }

  gStyle->SetPalette(kOcean);

  hs->Add(hist_mu);
  hs->Add(hist_e);
  hs->Add(hist_tau);
  hs->Add(hist_mu_anti);
  hs->Add(hist_e_anti);
  hs->Add(hist_tau_anti);

  hs->Draw();

  gPad->Update();
  //TPaveStats *st = (TPaveStats*)hist_mu->FindObject("stats");

  //stat box position
  //st->SetX2NDC(0.6); 
  //st->SetY2NDC(0.15); 
  //st->SetX1NDC(0.8); 
  //st->SetY1NDC(0.35); 

  //st->SetOptStat(111110110);

  myCanvas->Update();
  //myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/GraphSmearMatrixCorrectedNewAxis.svg");
}