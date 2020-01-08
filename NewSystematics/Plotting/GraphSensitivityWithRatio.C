#include <fstream>
#include<cmath>
int GraphSensitivityWithRatio(){

  double x[31];// array that can hold 100 numbers for 1st column 
  double y[31];// array that can hold 100 numbers for 2nd column 

  ifstream infile;   

  infile.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewSystematics/Plotting/Input/allaCorr2.dat");// file containing numbers in 3 columns 
     if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
       infile.ignore(10000,'\n');
       for(int k=0;k<31;++k) // reads file to end of *file*, not line
      { 
         infile >> x[k]; // read first column number
         infile >> y[k]; // read second column number
         printf("%f\n",x[k]);


         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
double y_new[31];
for(int u=0;u<31;++u){
  y_new[u] = pow(10.0,y[u]);

}

double x2[31];// array that can hold 100 numbers for 1st column 
  double y2[31];// array that can hold 100 numbers for 2nd column 

  ifstream infile2;   

  infile2.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewSystematics/Plotting/Input/allaCorrDouble.dat");// file containing numbers in 3 columns 
     if(infile2.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
       infile2.ignore(10000,'\n');
       for(int k=0;k<31;++k) // reads file to end of *file*, not line
      { 
         infile2 >> x2[k]; // read first column number
         infile2 >> y2[k]; // read second column number
         printf("%f\n",x2[k]);


         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
double y2_new[31];
for(int u=0;u<31;++u){
  y2_new[u] = pow(10.0,y2[u]);

}

  //TCanvas *c1 = new TCanvas("c1","3 Graphs",700,900);

  
  auto myCanvas = new TCanvas();
  // myCanvas->SetGrid();

  auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3); p2->Draw();
  p2->SetTopMargin(0.001);
  p2->SetBottomMargin(0.3);
  p2->SetGrid();
  p2->SetLogx ();
  p2->SetLogy();

  auto *p1 = new TPad("p1","p1",0.,0.3,1.,0.9);  p1->Draw();
  p1->SetBottomMargin(0.001);
  p1->cd();
  p1->SetGrid();
  p1->SetLogx ();
  //p1->SetLogy();

  auto spa = new TGraph(31,x,y_new);
  auto spb = new TGraph(31,x2,y2_new);

  //TGraphErrors* sp = new TGraphErrors(g->GetN(),pow(g->GetY(),10.0),g->GetX());

  spa->SetMarkerStyle(4);
  spa->SetMarkerColor(4);
  spb->SetMarkerStyle(3);
  spb->SetMarkerColor(3);
  
  auto mg = new TMultiGraph();

  //gStyle->SetTitleFontSize(0.04);
  mg->SetTitle("Sensitivity Curve Comparison for Detector Configurations With GLoBES Correlatons.");
  mg->Add(spa);
  mg->Add(spb);
  //mg->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  //mg->GetYaxis()->SetTitle("sin(2*theta_13)^2 sensitiivity");
  //mg->GetXaxis()->CenterTitle(true);
  //mg->GetYaxis()->CenterTitle(true);
  // sp->SetTitle("A: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  // sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  // sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  // sp->SetMarkerStyle(4);
  // sp->SetMarkerColor(4);

  //gPad->SetLogx();
  //gPad->Update();
  mg->Draw("APL");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  legend->AddEntry(spa,"series A: Both 700 km baseline","lp");
  legend->AddEntry(spb,"series B: Single detector, double fiducial mass","lp");
  legend->Draw();

  //// Ratio Plot
  double y_ratio[31];
  for(int z=0; z<31; z++){
    y_ratio[z] = y_new[z]/y2_new[z];
  }

  p2->cd();
  auto r2 = new TGraph(31,x,y_ratio);
  auto mg2 = new TMultiGraph();
  mg2->Add(r2);
  mg2->Draw("APL");
  printf("test %f",x[30]);
  //mg2->GetXaxis()->SetLimits(0.,x[30]);
  //mg->GetXaxis()->SetLimits(0.,x[30]);

  mg2->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  mg2->GetXaxis()->CenterTitle(true);
  mg2->GetXaxis()->SetTitleSize(0.1);
  mg->GetYaxis()->SetTitle("Chi Squared Sensitivity");
  mg->GetYaxis()->CenterTitle(true);
  mg2->GetYaxis()->SetTitle("Ratio A/B");
  mg2->GetYaxis()->CenterTitle(true);
  mg2->GetYaxis()->SetTitleSize(0.06);



  //  TGraph*r = new TGraph(31); r->SetTitle("");
  //  //r->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  //  for (int i=0; i<31; i++) r->SetPoint(i, x[i], y_new[i]/y2_new[i]);
  //     r->GetXaxis()->SetLabelSize(0.075);
  //     r->GetYaxis()->SetLabelSize(0.075);
  //  r->Draw("AL");


   


  myCanvas->Update();
  
  myCanvas->Modified();
  myCanvas->Update();


  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewSystematics/Plotting/Output/SensitivityplotAChiComparedRatio.svg");

return 0; // everything went right.
}
