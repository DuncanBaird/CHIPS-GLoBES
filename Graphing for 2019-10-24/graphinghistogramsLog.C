TH1* chi2Hist = NULL;

void graphinghistogramsLog(){

  TFile* f = TFile::Open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TBranch *b1 = (TBranch*)t1->GetBranch("TimeStamp_ns");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  chi2Hist = new TH1D("Time events modulo 50k","Histrogram of events over time",100,0,500); //last numbers: number of bins, lower bound of bins, upper bound of bins
  chi2Hist->GetXaxis()->SetTitle("modulo(Timestamp ns)");
  chi2Hist->GetYaxis()->SetTitle("Log(Frequency of event timestamp)");

  UInt_t TimeStamp_ns;
  t1->SetBranchAddress("TimeStamp_ns",&TimeStamp_ns);


  int entries = t1->GetEntries();
  for (int i=0;i<entries;i++){
    t1->GetEntry(i);
    chi2Hist->Fill(TimeStamp_ns%50000); //take modulo 50000 of each timestamp
  }

  chi2Hist->Draw();
  gPad->SetLogy();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)chi2Hist->FindObject("stats");

  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/GraphModulofiftythouLog.pdf");

  //printf("%d",3 % 2);
}