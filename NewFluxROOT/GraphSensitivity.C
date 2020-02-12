#include <fstream>
#include<cmath>
int GraphSensitivity(){

  double x[31];// array that can hold 100 numbers for 1st column 
  double y[31];// array that can hold 100 numbers for 2nd column 

  ifstream infile;   

  infile.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewFluxROOT/insima.dat");// file containing numbers in 3 columns 
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



  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();
  TGraph* sp = new TGraph(31,x,y_new);

  //TGraphErrors* sp = new TGraphErrors(g->GetN(),pow(g->GetY(),10.0),g->GetX());

  sp->SetTitle("A: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  sp->SetMarkerStyle(4);
  sp->SetMarkerColor(4);

  

  gPad->SetLogx();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)sp->FindObject("stats");
  
  
  sp->Draw();

  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewFluxROOT/Plotting/SensitivityplotA.pdf");

return 0; // everything went right.
}
