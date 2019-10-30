void GraphSensitivity(){

  TFile* f = TFile::Open("/home/duncan/Documents/Globes Sample Data Analysis/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  TGraphErrors* sp = new TGraphErrors("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/InitialSimulation/insimb.dat","%lg %lg");
  
  sp->SetTitle("B: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  sp->SetMarkerStyle(4);
  sp->SetMarkerColor(4);

  gPad->SetLogx();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)sp->FindObject("stats");
  
  
  sp->Draw();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/InitialSimulation/Plotting/SensitivityplotB.pdf");
}
