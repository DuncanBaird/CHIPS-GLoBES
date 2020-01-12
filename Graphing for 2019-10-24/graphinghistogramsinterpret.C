
void graphinghistogramsinterpret(){

  TFile* f = TFile::Open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TBranch *b1 = (TBranch*)t1->GetBranch("TimeStamp_ns");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

  UInt_t TimeStamp_ns;
  t1->SetBranchAddress("TimeStamp_ns",&TimeStamp_ns);

  t1->Draw("This->GetReadEntry():TimeStamp_ns", "Channel == 1.0");
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->SetTitle("Plot of events against timestamp ns for Channel = 1");
  htemp->GetXaxis()->SetTitle("Time stamp ns");
  htemp->GetYaxis()->SetTitle("Event number");
  
  gPad->Update();
  myCanvas->Update();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/GraphModuloInterpret.pdf");
}