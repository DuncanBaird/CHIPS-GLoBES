#include <cmath>

double chisquared(double data[],double model[]){
  double chi2 = 0;
  int len = *(&data + 1) - data;
  double mean;
  double x_sum;

  for (int i = 0; i<len;i++){
    x_sum += data[i];
  }
  mean = x_sum/len;

  for (int i = 0; i<len;i++){
    chi2 += pow((data[i]-model[i])/(pow((pow((data[i]-mean),2)),0.5))/(len-1),2);
    
  }


  return chi2;
}


void covariancescratch(){

  //DATA stuff//
  const int n = 5;
  double data_y [n] = {7,10,13,16,19};
  double data_x [n] = {1,2,3,4,5};

  double m = 3;
  double b = 4;

  double model_u[n] = {6,11,14,15,18};

  //Plotting Stuff//
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

  TGraph* gr1 = new TGraph(n,data_x,data_y);
  TGraph* gr2 = new TGraph(n,data_x,model_u);

  gr1->SetMarkerStyle(4);
  gr1->SetMarkerColor(4);
  gr1->SetLineColor(4);
  gr2->SetMarkerStyle(3);
  gr2->SetMarkerColor(3);
  gr2->SetLineColor(3);

  TMultiGraph* mg = new TMultiGraph();

  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("AC*");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  legend->AddEntry(gr1,"Data","lp");
  legend->AddEntry(gr2,"Model","lp");
  legend->Draw();

  printf("Chi is %f\n",chisquared(data_x,model_u));

}