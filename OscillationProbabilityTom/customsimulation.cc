/*
 * Duncan Baird; Make probability plots.
 * 
 * Compile with ''make customsimulation''
 */

#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include "initialsimulation.h"
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <iostream>
#include <stdio.h>
using namespace std;


/* "True" oscillation parameters */
double theta12,theta13,theta23,deltacp,sdm,ldm;

/* GLoBES parameter structures */
glb_params true_values;
glb_params test_values;
glb_params input_errors;

#define EXP_FAR  0
#define EXP_NEAR 1


const int baseline_steps = 6000;

void getOscillation(double osc_data[][baseline_steps],double min_base,double max_base,double energy,int l,int m,int option){
  double baseline = ((max_base-min_base)/baseline_steps)*1000;
  if(option ==1){
  for(int i = 0;i< baseline_steps;i++){
    osc_data[0][i] = (i * baseline)/(energy);
    osc_data[1][i] = glbVacuumProbability(l, m, 1,energy,i*baseline);
    // osc_data[0][i] = baseline * i;
    // osc_data[1][i] = i*energy;
    }
  }else if(option==0){
    double energy_step = energy/baseline_steps;
    for(int i = 0;i< baseline_steps;i++){
    osc_data[0][i] = i * energy_step;
    osc_data[1][i] = glbProfileProbability(EXP_FAR,l, m, 1,i*energy_step);
  }
  }
  
}

/***************************************************************************
 *                   P L O T T I N G   F U N C T I O N                     *
 ***************************************************************************/

/*
* Function to create and plot oscillation probabilities
* l is first lepton number, m is second lepton number, 
* option turns plot saving on and off.
*/
void doPlotROOTProb(int l1, int m1, int l2,int m2,int option,const char* canvas_name,int plot_type){


  double plot_data_prob_1[2][baseline_steps];
  double plot_data_prob_2[2][baseline_steps];

// for (int g = 0;g<baseline_steps; g++){
//     printf("Data a is: index: %f and value: %f \n",plot_data_a[0][g],plot_data_a[1][g]);
//     printf("Data b is: index: %f and value: %f \n",plot_data_b[0][g],plot_data_b[1][g]);
//   }
char *plottitle = NULL;
char *xtitle = NULL;

  if(plot_type == 1){
    plottitle = (char*)"Plot of Vacuum Probabilities for 10 GeV";
    xtitle = (char*)"Baseline L/E km/GeV";


    getOscillation(plot_data_prob_1,10.0,5E2,10,l1,m1,1);
    getOscillation(plot_data_prob_2,10.0,5E2,10,l2,m2,1);

  }else if(plot_type == 0){
    plottitle = (char*)"Plot of Matter Probabilities";
    xtitle = (char*)"E GeV";

  getOscillation(plot_data_prob_1,100.0,4E2,10,l1,m1,0);
  getOscillation(plot_data_prob_2,100.0,4E2,10,l2,m2,0);

  }else{
    plottitle = (char*)"error1";
    xtitle = (char*)"errorx1";
  }
  plot_data_prob_1[1][0] = plot_data_prob_1[1][1];
  plot_data_prob_2[1][0] = plot_data_prob_2[1][1];

  auto myCanvas = new TCanvas(canvas_name,canvas_name);
  myCanvas->SetGrid();
  auto spa = new TGraph(baseline_steps,plot_data_prob_1[0],plot_data_prob_1[1]);
  auto spb = new TGraph(baseline_steps,plot_data_prob_2[0],plot_data_prob_2[1]);

  //TGraphErrors* sp = new TGraphErrors(g->GetN(),pow(g->GetY(),10.0),g->GetX());

  spa->SetMarkerStyle(4);
  spa->SetMarkerColor(4);
  spb->SetMarkerStyle(3);
  spb->SetMarkerColor(3);
  
  auto mg = new TMultiGraph();

  mg->SetTitle(plottitle);
  mg->Add(spa);
  mg->Add(spb);
  mg->GetXaxis()->SetTitle(xtitle);
  mg->GetYaxis()->SetTitle("Probability");
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);
  // sp->SetTitle("A: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  // sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  // sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  // sp->SetMarkerStyle(4);
  // sp->SetMarkerColor(4);

  //gPad->SetLogx();
  gPad->Update();
  mg->Draw("APL");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  char dummy1[30];
  sprintf(dummy1,"series A: %d->%d",l1,m1);
  char dummy2[30];
  sprintf(dummy2,"series A: %d->%d",l2,m2);
  legend->AddEntry(spa,dummy1,"lp");
  legend->AddEntry(spb,dummy2,"lp");
  legend->Draw();

  myCanvas->Update();
  
  myCanvas->Modified();
  myCanvas->Update();

  if (option == 1){

  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/OscillationProbability/Plotting/Output/";
  string file_svg = ".svg";
  string file_pdf = ".pdf";

  string filestring1 = filepath + (string)canvas_name + file_svg;
  char filename1[filestring1.length() + 1];
  strcpy(filename1,filestring1.c_str());

  string filestring2 = filepath + (string)canvas_name + file_pdf;
  char filename2[filestring2.length() + 1];
  strcpy(filename2,filestring2.c_str());

  myCanvas->SaveAs(filename1);
  myCanvas->SaveAs(filename2);
  }

}

void userConfirm(){
  int kz;
  cout << "Please enter binary true to proceed: ";
  cin >> kz;
  cout << "The value you entered is " << kz;
  cout << " and its double is " << kz*2 << ".\n";
  //return 0;
}



/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
  //TApplication app{"customsimulation", &argc, argv};
  TApplication *app = new TApplication("app", &argc, argv);

  TFile *file = new TFile("test.root", "RECREATE");
  file->Write();
  file->Close();
  delete file;

  double *old_sys_errors = NULL;      /* Temp. pointer to systematical error array */
  int sys_dim;                        /* Abbrv. for number of systematical errors */
  int n_bins=86;                      /* Number of bins */
  int i;

  /* Set standard oscillation parameters (cf. hep-ph/0405172v5) */
  theta12 = asin(sqrt(0.3));
  theta13 = 0.0;
  theta23 = M_PI/4;
  deltacp = M_PI/2;
  sdm = 7.9e-5;
  ldm = 2.6e-3;

  glbInit(argv[0]);                    /* Initialize GLoBES and define chi^2 functions */
  

  /* Load 2 experiments: DC far (#0) and near (#1) detectors */
  char *Far_file = (char*)"CHIPS-GLB/glb-CHIPS10-7mrad-ME-FAR.glb";
  char *Near_file = (char*)"CHIPS-GLB/glb-CHIPS10-7mrad-ME-NEAR.glb";


  glbClearExperimentList();
  glbInitExperiment(Far_file, &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment(Near_file, &glb_experiment_list[0], &glb_num_of_exps);
  if (glbGetNumberOfBins(EXP_FAR) != n_bins || glbGetNumberOfBins(EXP_NEAR) != n_bins)
  {
    printf("ERROR: Number of bins changed in AEDL file, but not in C code (or vice-versa).\n");
    printf("nbins %i",n_bins);
    printf("far bins %i",glbGetNumberOfBins(EXP_FAR));
    printf("near bins %i",glbGetNumberOfBins(EXP_NEAR));

    return -1;
  }
  else
    n_bins = glbGetNumberOfBins(EXP_FAR);

  /* Initialize parameter vectors */
  true_values  = glbAllocParams();
  test_values  = glbAllocParams();
  input_errors = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  glbDefineParams(input_errors, 0.1*theta12, 0, 0.15*theta23, 0, 0.05*sdm, 0.05*ldm);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  glbSetOscillationParameters(true_values);
  glbSetInputErrors(input_errors);

/***************************************************************************
 *                                P L O T T I N G                          *
 ***************************************************************************/

doPlotROOTProb(3,3,2,2,0,"OSCProbabilityPlotSurvMUvTAU",1);
doPlotROOTProb(1,1,2,2,0,"OSCProbabilityPlotSurvMUvE",1);

userConfirm();

doPlotROOTProb(2,2,2,1,0,"OSCProbabilityProfilePlotv3",0);

  app->Run();

  /* Clean up */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);

  return 0;
}