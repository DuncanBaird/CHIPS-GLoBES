// Code to read flux data from .root file
// and export to .dat file.

void ReadFluxROOT(){

TFile f("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/New Flux Files/FluxReader_actualchipslocation_0.root");

//f.cd("enufullfine/CHIPSoffAXIS");
//f.cd("CHIPSoffAXIS");


TH1D *h_nue = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_nue_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_anue = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_anue_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_numu = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_numu_allpar_NoXSec_CHIPSoffAXIS");
TH1D *h_anumu = (TH1D*)f.Get("enufullfine/CHIPSoffAXIS/enufullfine_anumu_allpar_NoXSec_CHIPSoffAXIS");

h_nue->SetDirectory(nullptr);
h_anue->SetDirectory(nullptr);
h_numu->SetDirectory(nullptr);
h_anumu->SetDirectory(nullptr);
f.Close();

TCanvas *myCanvas = new TCanvas();
myCanvas->SetGrid();


h_nue->Draw();
h_anue->Draw("SAME");
h_numu->Draw("SAME");
h_anumu->Draw("SAME");


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
