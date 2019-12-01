void MultiGraphSensitivity(){

  TFile* f = TFile::Open("/home/duncan/Documents/Globes Sample Data Analysis/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  TGraphErrors* spa = new TGraphErrors("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/InitialSimulation/insima.dat","%lg %lg");
  
  TGraphErrors* spb = new TGraphErrors("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/InitialSimulation/insimb.dat","%lg %lg");

  TMultiGraph* mg = new TMultiGraph();

  spa->SetMarkerStyle(4);
  spa->SetMarkerColor(4);
  spb->SetMarkerStyle(3);
  spb->SetMarkerColor(3);

  mg->SetTitle("AB: Log Plot of sensitivity curve for Initial Simulation with systematic corrections For CHIPS .glb");
  mg->Add(spa);
  mg->Add(spb);
  mg->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  mg->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);



  gPad->SetLogx();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)spa->FindObject("stats");
  
  mg->Draw("APL");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  legend->AddEntry(spa,"series A: No Systematics","lp");
  legend->AddEntry(spb,"series B: All systematics (no spectra shift)","lp");
  legend->Draw();



  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/InitialSimulation/Plotting/SensitivityplotABnew.pdf");
}
