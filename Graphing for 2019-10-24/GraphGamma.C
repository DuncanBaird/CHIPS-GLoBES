void GraphGamma(){

  TFile* f = TFile::Open("/home/duncan/Documents/Globes Sample Data Analysis/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  t1->Draw("ToT:Channel");
  myCanvas->cd(2);
  t1->Draw("ToT:Channel", "Channel == 1.0");
  myCanvas->cd(3);
  t1->Draw("ToT:TimeStamp_ns", "Channel == 1.0");
  myCanvas->cd(4);
  t1->Draw("ToT:TimeStamp_s");
  myCanvas->SaveAs("/home/duncan/Documents/Globes Sample Data Analysis/Graphs/GraphGamma3.pdf");
}
