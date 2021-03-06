/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * Duncan Baird; Initial Simultion; based on example5
 * Example: Create sensitivity plot for reactor experiments with near
 * and far detectors.
 * Compile with ``make example5''
 */

#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include "initialsimulation.h"


/* If filenames given, write to file; for empty filenames, write to screen */
char MYFILE1[]="Plotting/Input/allaCorr100.dat";
char MYFILE2[]="Plotting/Input/allbCorr100.dat";
// char MYFILE3[]="test5c.dat";
// char MYFILE4[]="test5d.dat";
// char MYFILE5[]="test5e.dat";

/* "True" oscillation parameters */
double theta12,theta13,theta23,deltacp,sdm,ldm;

/* GLoBES parameter structures */
glb_params true_values;
glb_params test_values;
glb_params input_errors;

/* GSL stuff */
const gsl_root_fsolver_type *T;   /* GSL algorithm for root finding */
gsl_root_fsolver *s;              /* GSL Solver object */
gsl_function gsl_func;            /* Function to minimize */

/* Experiment and analysis parameters */
const double min_runtime = 1e-2;  /* Minimum running time to consider [years] */
const double max_runtime = 1e4;   /* Maximum running time to consider [years] */
const double chi2_goal   =  2.7;  /* Desired chi^2 value in root finder (2.7 = 90% C.L.) */
const int tSteps = 30;            /* Number of data points for each curve */
const double logs22th13_precision = 0.005; /* Desired precision of log(sin[th13]^2) in root finder */

#define MAX_SYS 200
double sys_errors[MAX_SYS];       /* Uncertainties of systematics params */
double sys_startval[MAX_SYS];     /* Starting values for systematics minimizer */
double sigma_binbin = 0.0;        /* Bin-to-bin error */

#define EXP_FAR  0
#define EXP_NEAR 1


/***************************************************************************
 *                        H E L P E R   F U N C T I O N S                  *
 ***************************************************************************/

/* Minimum of two numbers */
inline double min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

/* Square of real number */
inline double square(double x)
{
  return x*x;
}

/* Gauss likelihood (this is sufficient here due to the large event numbers
 * in a reactor experiment; for other setups, one should use Poisson statistics) */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}



/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/

/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors:                                                                 *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 ***************************************************************************/
double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double signal_fit_rates_N[n_bins];
  double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  /* Apply energy calibration error */
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
                      signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
                      signal_fit_rates_N, n_bins, emin, emax);

  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical part of chi^2 for far detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[1] */
    fit_rate  = signal_norm_F * signal_fit_rates_F[i];
    chi2 += likelihood(true_rates_F[i], fit_rate, true_rates_F[i]);

    /* Statistical part of chi^2 for near detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[2] */
    fit_rate  = signal_norm_N * signal_fit_rates_N[i];
    chi2 += likelihood(true_rates_N[i], fit_rate, true_rates_N[i]);
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    sys_startval[i] = x[i];

  return chi2;
}


/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors:                                                                 *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 *   x[5]...x[5+nBins-1]: Spectral error                                   *
 * and the bin-to-bin error sigma_binbin. The calculation is based on      *
 * Eq. (9) of hep-ph/0303232.                                              *
 ***************************************************************************/
double chiDCSpectral(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double signal_fit_rates_N[n_bins];
  double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  /* Apply energy calibration error */
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
                      signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
                      signal_fit_rates_N, n_bins, emin, emax);

  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical part of chi^2 for far detector
     * Normalization is affected by flux error x[0], fiducial mass error x[1],
     * spectral perturbations x[5+i], and bin-to-bin error sigma_binbin */
    fit_rate  = (signal_norm_F + x[5+i]) * signal_fit_rates_F[i];
    chi2 += likelihood( true_rates_F[i], fit_rate,
                    true_rates_F[i] * (1.0 + true_rates_F[i]*square(sigma_binbin)) );

    /* Statistical part of chi^2 for near detector
     * Normalization is affected by flux error x[0], fiducial mass error x[2],
     * spectral perturbations x[5+i], and bin-to-bin error sigma_binbin */
    fit_rate  = (signal_norm_N + x[5+i]) * signal_fit_rates_N[i];
    chi2 += likelihood( true_rates_N[i], fit_rate,
                    true_rates_N[i] * (1.0 + true_rates_N[i]*square(sigma_binbin)) );
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    sys_startval[i] = x[i];

  return chi2;
}


/***************************************************************************
 *          W R A P P E R   F O R   G S L   R O O T   F I N D E R          *
 ***************************************************************************/

/* Callback function for GSL root finder */
double DoChiSquare(double x, void *dummy)
{
  double thetheta13, chi2;

  /* Set vector of test values */
  thetheta13 = asin(sqrt(pow(10,x)))/2;
  glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
  glbSetOscParams(test_values, 200.0/2 * (x+4)*M_PI/180, GLB_DELTA_CP);

  /* Set starting values for systematics minimiyer to the coordinates of
   * minimum in the last iteration. This accelerates the minimization and
   * prevents convergence problems. */
  glbSetSysStartingValuesList(EXP_FAR, 0, GLB_ON, sys_startval);

  /* Compute Chi^2 for all loaded experiments and all rules
   * Correlations are unimportant in reactor experiments, so glbChiSys is sufficient */
  chi2 = glbChiTheta13(test_values,NULL, GLB_ALL);

  return chi2 - chi2_goal;
}


/* Error handler for GSL errors */
void gslError(const char *reason, const char *file, int line, int gsl_errno)
{
  static int n_errors=0;

  printf("GSL Error in file %s:%d : %s\n", file, line, gsl_strerror(gsl_errno));
  if (++n_errors > 1000)
  {
    printf("Too many GSL errors. Program aborted.\n");
    abort();
  }
}


/* Compute sensitivity to sin^2(2*th13) for different running times */
void ComputeSensitivityCurve()
{
  double t_factor = pow(max_runtime/min_runtime, 1.0/tSteps);
  double t;
  double x;               /* Current value of log[sin^2(2*th13)] */
  double x_lo, x_hi;      /* Bracketing interval for root finder */
  double x_sens[tSteps+1];/* Calculated sensitivities to log[sin^2(2*th13)] */
  int gsl_status;         /* GSL error code */
  int iter;               /* Iteration counter for root finder */
  const int max_iter=100; /* Maximum number of iterations allowed in root finder */
  int j;

  for (j=0; j <= tSteps; j++)      /* Initialize output vector */
    x_sens[j] = 1;                 /* As x is always <= 0 this signals "no value" */

  /* Loop over all data points in this dataset, using log stepping */
  for (j=0; j <= tSteps; j++)
  {
    /* Set running time */
    t = min_runtime * pow(t_factor,j);
    glbSetRunningTime(EXP_FAR, 0, t);
    glbSetRunningTime(EXP_NEAR, 0, t);

    /* Calculate "true" event rates */
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbDefineParams(input_errors, 0.1*theta12, 0, 0.15*theta23, 0, 0.05*sdm, 0.05*ldm);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetCentralValues(true_values);
    glbSetInputErrors(input_errors);
    glbSetRates();

    /* Determine sensitivity to sin^2(2*th13) by  using the GSL Brent-Dekker
     * algorithm to find the th13, for which chi^2 = 2.7 (90%) */
    x_lo = -3.0;
    x_hi = -0.1;
    iter = 0;

    /* Start root finder. The initial search interval is guessed based on the
     * results from a neighboring grid point */
    double deviation=0.1;
    do {
      if (j > 0) {
        x_lo = x_sens[j-1] - deviation;
        x_hi = min(x_sens[j-1] + deviation, -0.001);
      }
      gsl_status = gsl_root_fsolver_set(s, &gsl_func, x_lo, x_hi);
      deviation *= 2;
    } while (gsl_status != GSL_SUCCESS);

    /* Iterate root finder */
    do {
      gsl_status = gsl_root_fsolver_iterate(s);
      x          = gsl_root_fsolver_root(s);
      x_lo       = gsl_root_fsolver_x_lower(s);
      x_hi       = gsl_root_fsolver_x_upper(s);
      gsl_status = gsl_root_test_interval (x_lo, x_hi, logs22th13_precision, 0);
    } while (gsl_status == GSL_CONTINUE && iter < max_iter);

    /* Save results */
    x_sens[j] = x;
    AddToOutput2(t, x);
  }
}

const int max_index = tSteps*8;

void ComputeSensitivityCurve2(double plot_data[][max_index])
{
  double t_factor = pow(max_runtime/min_runtime, 1.0/tSteps);
  double t;
  double x;               /* Current value of log[sin^2(2*th13)] */
  double x_lo, x_hi;      /* Bracketing interval for root finder */
  double x_sens[tSteps+1];/* Calculated sensitivities to log[sin^2(2*th13)] */
  int gsl_status;         /* GSL error code */
  int iter;               /* Iteration counter for root finder */
  const int max_iter=100; /* Maximum number of iterations allowed in root finder */
  int j;

 
  for(int w = 0; w < 8; w++){
    printf("w value: %d \n",w);
    double baseline = (w*100.0) + 200.0;
    printf("Starting baseline: %f \n",baseline);

    glbSetBaselineInExperiment(EXP_NEAR,baseline);
    glbSetBaselineInExperiment(EXP_FAR,baseline);

    printf("Baselines set \n");
  


  for (j=0; j <= tSteps; j++)      /* Initialize output vector */
    x_sens[j] = 1;                 /* As x is always <= 0 this signals "no value" */

  /* Loop over all data points in this dataset, using log stepping */
  for (j=0; j <= tSteps; j++)
  {
    /* Set running time */
    t = min_runtime * pow(t_factor,j);
    glbSetRunningTime(EXP_FAR, 0, t);
    glbSetRunningTime(EXP_NEAR, 0, t);

    /* Calculate "true" event rates */
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbDefineParams(input_errors, 0.1*theta12, 0, 0.15*theta23, 0, 0.05*sdm, 0.05*ldm);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetCentralValues(true_values);
    glbSetInputErrors(input_errors);
    glbSetRates();

    /* Determine sensitivity to sin^2(2*th13) by  using the GSL Brent-Dekker
     * algorithm to find the th13, for which chi^2 = 2.7 (90%) */
    x_lo = -3.0;
    x_hi = -0.1;
    iter = 0;

    /* Start root finder. The initial search interval is guessed based on the
     * results from a neighboring grid point */
    double deviation=0.1;
    do {
      if (j > 0) {
        x_lo = x_sens[j-1] - deviation;
        x_hi = min(x_sens[j-1] + deviation, -0.001);
      }
      gsl_status = gsl_root_fsolver_set(s, &gsl_func, x_lo, x_hi);
      deviation *= 2;
    } while (gsl_status != GSL_SUCCESS);

    /* Iterate root finder */
    do {
      gsl_status = gsl_root_fsolver_iterate(s);
      x          = gsl_root_fsolver_root(s);
      x_lo       = gsl_root_fsolver_x_lower(s);
      x_hi       = gsl_root_fsolver_x_upper(s);
      gsl_status = gsl_root_test_interval (x_lo, x_hi, logs22th13_precision, 0);
    } while (gsl_status == GSL_CONTINUE && iter < max_iter);

    /* Save results */
    x_sens[j] = x;

    //Baseline data
    plot_data[0][(w*tSteps)+j] = baseline;
    //Exposure data
    plot_data[1][(w*tSteps)+j] = t;
    //Chi data
    plot_data[2][(w*tSteps)+j] = x;
  }

  }
  // return 0;
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

  /* Initialization */
  for (i=0; i < MAX_SYS; i++)
    sys_startval[i] = 0.0;

  /* Set standard oscillation parameters (cf. hep-ph/0405172v5) */
  theta12 = asin(sqrt(0.3));
  theta13 = 0.0;
  theta23 = M_PI/4;
  deltacp = M_PI/2;
  sdm = 7.9e-5;
  ldm = 2.6e-3;

  glbInit(argv[0]);                    /* Initialize GLoBES and define chi^2 functions */
  glbDefineChiFunction(&chiDCNorm,     5,        "chiDCNorm",     NULL);
  glbDefineChiFunction(&chiDCSpectral, 5+n_bins, "chiDCSpectral", NULL);

  gsl_set_error_handler(&gslError);    /* Initialize GSL root finder */
  gsl_func.function = &DoChiSquare;
  gsl_func.params   = NULL;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

  /* Load 2 experiments: DC far (#0) and near (#1) detectors */
  glbClearExperimentList();
  glbInitExperiment("CHIPS-GLB/glb-CHIPS10-7mrad-ME-FAR.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("CHIPS-GLB/glb-CHIPS10-7mrad-ME-NEAR.glb", &glb_experiment_list[0], &glb_num_of_exps);
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

  // /* Calculate sensitivity curve without systematics */
  // InitOutput(MYFILE1,"Format: Running time   Log(10,s22th13) sens. \n");
  glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
  // ComputeSensitivityCurve();
  // glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);

  // /* Calculate sensitivity curve with normalization and energy calibration errors,
  //  * as defined in the AEDL files */
  // InitOutput(MYFILE2,"Format: Running time   Log(10,s22th13) sens. \n");
  // ComputeSensitivityCurve();


  printf("Before curve \n");
  double plot_data[3][max_index];
  ComputeSensitivityCurve2(plot_data);
  printf("After curve");


  /* PLOTTING */
  auto myCanvas = new TCanvas();

  TH2* h2 = new TH2D("Experiment viability", "Surface plot of Chi Sensitivity Surface", 8, 200, 900, tSteps, plot_data[1][0], plot_data[1][29]);

  myCanvas->SetGrid();
  h2->GetXaxis()->SetTitle("Baseline km");
  h2->GetYaxis()->SetTitle("Integrated detector luminosity GW t years");
  h2->GetZaxis()->SetTitle("Chi Squared Sensitivity of 23 mixing angle");
  
  
  for(int q = 1; q < tSteps*8; ++q) {
    //printf("Data inspection: index: %d, pow chi: %f, baseline: %f, exposure: %f \n",q,pow(10.0,plot_data[2][q]),plot_data[0][q],plot_data[1][q]);
  
  
  h2->Fill(plot_data[0][q],plot_data[1][q],pow(10.0,plot_data[2][q]));
  }

  h2->Draw("SURF");
  gPad->SetLogz();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
  myCanvas->Update();
  //myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/CHIPShist-First.pdf");

  app->Run();




  /* Clean up */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  gsl_root_fsolver_free(s);

  return 0;
}