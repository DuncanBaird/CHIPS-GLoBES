/*
 * Duncan Baird; Simulation for doing several plots.
 * Example: Create sensitivity plot for reactor experiments with near
 * and far detectors.
 * Compile with ``make customsimulation''
 */

#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include "initialsimulation.h"
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TError.h>
#include <TLatex.h>
#include <iostream>
#include <stdio.h>
using namespace std;


/* "True" oscillation parameters */
double theta12,theta13,theta23,deltacp,sdm,ldm;

/* GLoBES parameter structures */
glb_params true_values;
glb_params test_values;
glb_params input_errors;


    /* Bin-to-bin error */

#define EXP_FAR  0
#define EXP_NEAR 1

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{
  //TApplication app{"customsimulation", &argc, argv};
  TApplication *app = new TApplication("app", &argc, argv);
  gErrorIgnoreLevel = kFatal;
 glbInit(argv[0]); // needed to get libraries


/*
theta12:  0.586168 +- 0.0133455 (s2t12 = 0.306)
theta13: 0.147693 +- 0.0026509 (s2t13 = 0.021655)
theta23:  0.72626 +- 0.041681 (0.872843 +- 0.041681)
  (s2t23 = 0.441 (0.587) )
d21:    7.5e-05 +- 1.77e-06
d31:    0.002524 +- 3.93e-05 (-0.002439 +- 3.73e-05)
*/
  /* Set standard oscillation parameters (cf. hep-ph/0405172v5) */

  //ASSUMING NORMAL HEIRACRCHY

  double theta12 = 0.586168;
  double theta13 = 0.147693;
  double theta23 = 0.72626;
  double deltacp = 0. ;
  double sdm = 7.5e-05;
  double ldm = 0.002524;

  // 1 sigma errors on parameters
  double e_theta12 = 0.0133455;
  double e_theta13 = 0.0026509;
  double e_theta23 = 0.041681;
  double e_sdm = 1.77e-06;
  double e_ldm = 3.93e-05;


  int n_bins=9;   /* Number of bins in AEDL Files*/



  /*   Initiates the Experiments      */
  glbClearExperimentList();
  glbInitExperiment("CHIPS-GLB/TOMCHIPS_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("CHIPS-GLB/TOMCHIPS_far.glb", &glb_experiment_list[0], &glb_num_of_exps);

  /*   Checks for Binning Consistancy     */
  if (glbGetNumberOfBins(EXP_FAR) != n_bins || glbGetNumberOfBins(EXP_NEAR) != n_bins)
  {
    printf("ERROR: Number of bins changed in AEDL file, but not in C code (or vice-versa) %d %d.\n" , glbGetNumberOfBins(EXP_FAR), glbGetNumberOfBins(EXP_NEAR));
    return -1;
  }
  else
    n_bins = glbGetNumberOfBins(EXP_FAR);
  
  printf("\n\nExperiment 0 Loaded : %s\nExperiment 1 Loaded : %s\n",glbGetFilenameOfExperiment(0) , glbGetFilenameOfExperiment(1) );


 true_values  = glbAllocParams();
  glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm);
  glbSetDensityParams(true_values, 1.0, GLB_ALL);

  input_errors = glbAllocParams();
  glbDefineParams(input_errors, e_theta12 ,e_theta13, e_theta23, 0.1, e_sdm, e_ldm);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);

  double deltacp_fit = 0.;
  test_values  = glbAllocParams(); 
  glbDefineParams(test_values, theta12 ,theta13, theta23, deltacp_fit, sdm, ldm);
  glbSetDensityParams(test_values, 1.0, GLB_ALL);

  //glbSetRates();
   glbSetOscParams(true_values, 3.141592654, GLB_DELTA_CP);
  glbSetOscillationParameters(true_values);
  glbSetRates();


cout << "Starting testing stuff \n";

cout << glbGetRunningTime(EXP_FAR,0) << "\n";

glbShowRuleRates(stdout,EXP_FAR,0,GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);

//runChiCurve(0,2*M_PI,100,1,1,0,"testing");

cout << "Finished testing stuff \n";

 printf("\n\nExperiment 0 Loaded : %s\nExperiment 1 Loaded : %s\n",glbGetFilenameOfExperiment(0) , glbGetFilenameOfExperiment(1) );

  app->Run();

  /* Clean up */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);


  return 0;
}