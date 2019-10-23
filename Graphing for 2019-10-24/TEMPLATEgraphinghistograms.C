TH1* chi2Hist = NULL;

void graphinghistograms(){

  TFile* f = TFile::Open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  chi2Hist = new TH1D("distribution name","Name of plot",100,4,20); //last numbers: number of bins, lower bound of bins, upper bound of bins
  chi2Hist->GetXaxis()->SetTitle("X axis name (bins)");
  chi2Hist->GetYaxis()->SetTitle("Yaxis Name (frequency of bin)");

  chi2Hist->Draw("test1");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)chi2Hist->FindObject("stats");
}