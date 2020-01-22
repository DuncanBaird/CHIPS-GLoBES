// Code to read flux data from .root file
// and export to .dat file.

void ReadFluxROOT(){

// Opening .root file
TFile f("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/New Flux Files/FluxReader_actualchipslocation_0.root");

//f.cd("enufullfine/CHIPSoffAXIS");
//f.cd("CHIPSoffAXIS");

// Opening Histograms (subdirectory specified here)
TH1D *h_nue = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_anue = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_anue_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_numu = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_numu_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_anumu = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_anumu_allpar_NoXSec_CHIPSoffAXIS");

h_nue->SetDirectory(nullptr);
h_anue->SetDirectory(nullptr);
h_numu->SetDirectory(nullptr);
h_anumu->SetDirectory(nullptr);
f.Close();

// Getting max number of bins
int bins_nue = h_nue->GetNbinsX();
int bins_anue = h_anue->GetNbinsX();
int bins_numu = h_numu->GetNbinsX();
int bins_anumu = h_anumu->GetNbinsX();

printf("Bins are %d\n",bins_nue);

// Arrays to hold values in memory
double energy [500];
double flux_nue [500];
double flux_anue [500];
double flux_numu [500];
double flux_anumu [500];
double flux_nutau [500];
double flux_anutau [500];

// Loop to get data from histograms
for(int i=1; i<502;i++){
  energy[i] = h_nue->GetBinLowEdge(i);
  flux_nue[i] = h_nue->GetBinContent(i);
  flux_anue[i] = h_anue->GetBinContent(i);
  flux_numu[i] = h_numu->GetBinContent(i);
  flux_anumu[i] = h_anumu->GetBinContent(i);

  flux_nutau[i] = 0;
  flux_anutau[i] = 0;

  printf("Energy is %f\n",energy[i]);


}

// Setup file output
ofstream outdata;

outdata.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/New Flux Files/example.dat");

// Loop to wrte each line of file
for(int i=1; i<502;i++){
  outdata << energy[i];
  outdata << " " << flux_nue[i];
  outdata << " " << flux_numu[i];
  outdata << " " << flux_nutau[i];
  outdata << " " << flux_anue[i];
  outdata << " " << flux_anumu[i];
  outdata << " " << flux_anutau[i] << endl;


}

outdata.close();

// TCanvas *myCanvas = new TCanvas();
// myCanvas->SetGrid();


// h_nue->Draw();
// h_anue->Draw("SAME");
// h_numu->Draw("SAME");
// h_anumu->Draw("SAME");


// enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS
// enufullfine_anue_allpar_NoXSec_CHIPSoffAXIS
// enufullfine_numu_allpar_NoXSec_CHIPSoffAXIS
// enufullfine_anumu_allpar_NoXSec_CHIPSoffAXIS


////OLD
  // TFile* f = TFile::Open("/home/duncan/Documents/Globes Sample Data Analysis/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  // TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  // TCanvas *myCanvas = new TCanvas();
  // myCanvas->SetGrid();
  // myCanvas->Divide(2,2);
  // myCanvas->cd(1);
  // t1->Draw("ToT:Channel");
  // myCanvas->cd(2);
  // t1->Draw("ToT:Channel", "Channel == 1.0");
  // myCanvas->cd(3);
  // t1->Draw("ToT:TimeStamp_ns", "Channel == 1.0");
  // myCanvas->cd(4);
  // t1->Draw("ToT:TimeStamp_s");
  // myCanvas->SaveAs("/home/duncan/Documents/Globes Sample Data Analysis/Graphs/GraphGamma3.pdf");

}
