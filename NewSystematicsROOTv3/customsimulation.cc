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
#include <iostream>
#include <stdio.h>
using namespace std;


/* If filenames given, write to file; for empty filenames, write to screen */
char MYFILE1[]="Plotting/Input/v2testa.dat";
char MYFILE2[]="Plotting/Input/v2testb.dat";
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

/************************
 * Interpolate
 * 
 */

void createMergedRate(TMatrixD data){
  int bins = glbGetNumberOfSamplingPoints(EXP_FAR);

  int far_rules = glbGetNumberOfRules(EXP_FAR);
  int near_rules = glbGetNumberOfRules(EXP_FAR);
  double *true_rate_0 = glbGetRuleRatePtr(EXP_FAR, 0);
  double *true_rate_1 = glbGetRuleRatePtr(EXP_FAR, 1);

for(int i = 0; i<bins;++i){
  data[i][0] = true_rate_0[i] + true_rate_1[i];
}

}


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

// Taken from glb_sys.c in GLoBES Library
double glb_likelihood(double true_rate, double fit_rate)
{
  double res;
  res = fit_rate - true_rate;
  if (true_rate > 0)
  {
    if (fit_rate <= 0.0)
      res = 1e100;
    else
      res += true_rate * log(true_rate/fit_rate);
  }
  else
    res = fabs(res);

  return 2.0 * res;

}

// Taken from glb_sys.c in GLoBES Library
double glb_prior(double x, double center, double sigma)
{
  double tmp = (x - center)/sigma;
  return tmp*tmp;
}

/***************************************************
 * Function to create covariance matrix
 * Pulled from covey.C
 * 
 ***************************************************/
void generate_covariance(TMatrixD covariance_matrix){
  
 const int mat_len = 200;
  //TMatrixD covariance_matrix(200,200);
  TMatrixD correlation_matrix(200,200);
  double_t matrix_data[200*200];

  //two detectors, close together, energy from 0-5

   TH2D *mat = new TH2D("mat","mat",200,0.0,10.0,200,0.0,10.0);
   int iuniverse=100;
   char no[5];

   TH1F *Euniverse1[100];
   TRandom *g = new TRandom();
   for (int i=0; i<iuniverse; i++)
     {
       sprintf (no, "no%d",i);
       Euniverse1[i] = new TH1F(no,no,200,0.0,10.0);
     }

   //generate the event spectra


   /*MINOS*/

   TFile *ff3 = new TFile("TrueE.root");
   TH1F *Espec1 = (TH1F*)ff3->Get("trueE_0to5");
   Espec1->SetName("Espec1");


   double c = 2.99792458e8;
   double pi = 3.1415926535;

   double t23 = 45.0/360.0*2.0*pi;
   double sint23 = sin(t23);
   double sin2t23 = sin(2.0*t23);
   double cost23 = cos(t23);
   double sinsqt23 = sint23*sint23;
   double cossqt23 = cost23*cost23;

   //   double t13 = 8.73/360.0*2.0*pi; //Daya Bay
   double t13 = 9.22/360.0*2.0*pi; //Average
   //double t13 = 15.0/360.0*2.0*pi; //code testing

   double sint13 = sin(t13);
   double sin2t13 = sin(t13*2.0);
   double sinsq2t13 = sin2t13*sin2t13;
   double cost13 = cos(t13);
   double cossqt13 = cost13*cost13;

   double t12 = 33.96/360.0*2.0*pi;
   double sin2t12 = sin(2.0*t12);
   double sinsq2t12 = sin2t12*sin2t12;
   double sint12 = sin(t12);
   double sinsqt12 = sint12*sint12;

   double Dmsq31NH = 2.45e-3;  //Sushant's values
   double Dmsq31IH = -2.31e-3;
   double Dmsq21 = 7.6e-5;

   double countchips, countminos, countt2k, countnova;
   double countminoschips, countminosnova, countt2knova, countall;
   double a = 1./4000.0;
   double ab = -1./4000.0;
   double AhatNHn = a/Dmsq31NH;
   double AhatIHn = a/Dmsq31IH;

   double AhatNHb = ab/Dmsq31NH;
   double AhatIHb = ab/Dmsq31IH;

   double alphaNH = Dmsq21/Dmsq31NH;
   double alphaIH = Dmsq21/Dmsq31IH;

   double deltacp=0.0;

   //make the NOVA spectrum
   TH1D *Espec2 = new TH1D("Espec2","Espec2",100,0.0,5.0);
   for (int i=0;i<10000;i++)
     {
       Espec2->Fill(g->Gaus(2.0,0.3));
     }

   //make the universes based on the real event spectrum
   std::cout<<" making the universes "<<std::endl;
   for (int k=0; k<100; k++) 
      {
	//use Fill instead of SetBinContent to allow for resolution smearing later
	//clear and reset the spectra histograms for each fake experiment

        for (int i=0; i<100; i++) //energy bins
          {
	    //E1 is MINOS or CHIPS, E2 is NOVA
	     double E1 = Espec1->GetBinCenter(i+1);
	     double Ewt1 = Espec1->GetBinContent(i+1);
	     double E2 = Espec2->GetBinCenter(i+1);
	     double Ewt2 = Espec2->GetBinContent(i+1);
	     //for each energy bin, throw a possion based on no entries
	     float i1 = g->Poisson(Ewt1);
	     float i2 = g->Poisson(Ewt2);
	     //now lets put in a constant E shift
	     E1 = E1+0.05*E1;
	     E2 = E2+0.05*E2;
	     //float rfake = g->Gaus(E,res[0]/sqrt(E));
	     //float Eres = E;//+rfake;
	     //fake11->Fill(Eres,ifake);

	     Euniverse1[k]->Fill(E1,i1);
	     Euniverse1[k]->Fill(E2+5,i2); //E2+5.0

	  } ///finished spectrum generation for this k universe
	std::cout<<" finished this universe "<<k<<std::endl;
	for (int i=0; i<200; i++)
	  {
	    for (int j=0; j<200; j++)
	      {
		double alan, mary;
		if(i<100)
		  {
		    alan = Espec1->GetBinContent(i+1) - Euniverse1[k]->GetBinContent(i);
		  }
		if(j<100)
		  {
		    mary = Espec1->GetBinContent(j+1) - Euniverse1[k]->GetBinContent(j);
		  }
		if(i>101)
		  {
		    alan = Espec2->GetBinContent(i+1) - Euniverse1[k]->GetBinContent(i);
		  }
		if(j>101)
		  {
		    mary = Espec2->GetBinContent(j+1) - Euniverse1[k]->GetBinContent(j);
		  }
		mat->SetBinContent(i+1,j+1,alan*mary/float(iuniverse));

    // storing covariance matrix data
    matrix_data[i*(mat_len-1)+j] = alan*mary/float(iuniverse);
    covariance_matrix[i][j] += alan*mary/float(iuniverse);
    correlation_matrix[i][j] += alan*mary/float(iuniverse);

    
		if(k==4&&i==120&&j==120)std::cout<<" i "<<i<<" j "<<j<<"alan "<<alan<<" mary"<<mary<<" j universe bin "<<Euniverse1[k]->GetBinContent(j+100)<<std::endl;
	      }
	  }
      }

    double singular_adjust = 0.001;
    for(int w=0;w<mat_len;w++){
      printf("diagonal %d value is %f \n",w,correlation_matrix[w][w]);
      correlation_matrix[w][w] = correlation_matrix[w][w] + singular_adjust;
      printf("diagonal %d value is now %f \n",w,correlation_matrix[w][w]);
    }
   //adding covariance values into matrix object
    // covariance_matrix.SetMatrixArray(matrix_data);
    // correlation_matrix.SetMatrixArray(matrix_data);

    // double_t det;
    // correlation_matrix.Invert(&det);

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

/********************************************************
 * My Custom chi squared function
 * Uses covariance matrix in computation
 ********************************************************/

TMatrixD covariance_matrix_1(200,200);

  


double chiCOV(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{
  double *true_rates       = glbGetRuleRatePtr(exp, rule);
  double *signal_fit_rates = glbGetSignalFitRatePtr(exp, rule);
  double *bg_fit_rates     = glbGetBGFitRatePtr(exp, rule);
  double bg_norm_center, bg_tilt_center;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  glbGetBGCenters(exp, rule, &bg_norm_center, &bg_tilt_center);
  for (i=ew_low; i <= ew_high; i++)
  {
    fit_rate = signal_fit_rates[i] + bg_norm_center * bg_fit_rates[i];
    chi2 += glb_likelihood(true_rates[i], fit_rate);
  }

  //NEW Contribution
  //define TMatrixD delta = fit - true;
  TMatrixD delta1(200,1);
  TMatrixD delta2(200,1);
  TMatrixD delta3(200,1);

  for(int i=0;i<200;i++){
    delta1[i][0] = 5E10;
    delta2[i][0] = 5E10;
  }

  createMergedRate(delta3);
  cout << "debug1";
  // delta.Transpose()
  delta1.T();
  cout << "debug2";
  //covariance
  double_t det1;
  //chi2+= delta_transpose * inverse_covariance * delta
  TMatrixD matrix_result = delta1 * covariance_matrix_1.Invert(&det1) * delta2;
  // TMatrixD dummy1;
  // TMatrixD dummy2;
  // cout << "debug4";
  // dummy1.Mult(delta1,covariance_matrix_1);
  // cout << "debug5";
  // dummy2.Mult(dummy1,delta2);
   cout << "debug6";
  double test_result = matrix_result[0][0];
  cout << "testing output: "<< test_result << "\n";
  //cout << "hello world\n";

  return chi2 +test_result;
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
  chi2 = glbChiSys(test_values,GLB_ALL,GLB_ALL);//glbChiTheta13(test_values,NULL, GLB_ALL);
  //cout << chi2<<"\n";
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

void ComputeSensitivityCurve3(double plot_data[][tSteps],int option)
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

  if(option == 1){
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
  } else if (option == 0) {
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
  } else {
    printf("No option available \n");
  }

 
  for(int w = 5; w < 6; w++){
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

    //Expsoure data
    plot_data[0][j] = t;
    //Chi data
    plot_data[1][j] = pow(10.0,x);
  }

  }
  // return 0;
}

/***************************************************************************
 *                   P L O T T I N G   F U N C T I O N                     *
 ***************************************************************************/

void doPlotROOT(double plot_data_a[][tSteps],double plot_data_b[][tSteps],int option,const char* canvas_name){

for (int g = 0;g<tSteps; g++){
    printf("Data a is: index: %f and value: %f \n",plot_data_a[0][g],plot_data_a[1][g]);
    printf("Data b is: index: %f and value: %f \n",plot_data_b[0][g],plot_data_b[1][g]);
  }
  plot_data_a[1][0] = plot_data_a[1][1];
  plot_data_b[1][0] = plot_data_b[1][1];

  auto myCanvas = new TCanvas(canvas_name,canvas_name);
  myCanvas->SetGrid();
  auto spa = new TGraph(tSteps,plot_data_a[0],plot_data_a[1]);
  auto spb = new TGraph(tSteps,plot_data_b[0],plot_data_b[1]);

  //TGraphErrors* sp = new TGraphErrors(g->GetN(),pow(g->GetY(),10.0),g->GetX());

  spa->SetMarkerStyle(4);
  spa->SetMarkerColor(4);
  spb->SetMarkerStyle(3);
  spb->SetMarkerColor(3);
  
  auto mg = new TMultiGraph();

  mg->SetTitle("AB: Log Plot of sensitivity curve systematics vs stats only.");
  mg->Add(spa);
  mg->Add(spb);
  mg->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  mg->GetYaxis()->SetTitle("sin(2*theta_13)^2 sensitiivity");
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);
  // sp->SetTitle("A: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  // sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  // sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  // sp->SetMarkerStyle(4);
  // sp->SetMarkerColor(4);

  gPad->SetLogx();
  gPad->Update();
  mg->Draw("APL");

  TLegend* legend = new TLegend();
  legend->SetHeader("Legend Title");
  legend->AddEntry(spa,"series A: With Systematics","lp");
  legend->AddEntry(spb,"series B: Statistics only","lp");
  legend->Draw();

  myCanvas->Update();
  
  myCanvas->Modified();
  myCanvas->Update();

  if (option == 1){

  string filepath = "/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/NewSystematicsROOTv2/Plotting/Output/";
  string file_svg = ".svg";
  string file_pdf = ".pdf";

  // char filename1[256];
  // strncat(filename1,filepath,sizeof(filename1));
  // printf(filename1);
  // strncat(filename1,canvas_name,sizeof(filename1));
  // strncat(filename1,file_svg,sizeof(filename1));
  // printf(filename1);
  
  // char filename2[256];
  // strncat(filename2,filepath,sizeof(filename2));
  // strncat(filename2,canvas_name,sizeof(filename2));
  // strncat(filename2,file_pdf,sizeof(filename2));

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
  cout << "Please enter an integer to proceed: ";
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
  gErrorIgnoreLevel = kFatal;

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

  // Defining my chi function in GLoBES
  glbDefineChiFunction(&chiCOV,     5,        "chiCOV",     NULL);

  gsl_set_error_handler(&gslError);    /* Initialize GSL root finder */
  gsl_func.function = &DoChiSquare;
  gsl_func.params   = NULL;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

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

  // /* Calculate sensitivity curve without systematics */
  // InitOutput(MYFILE1,"Format: Running time   Log(10,s22th13) sens. \n");
  glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
  // ComputeSensitivityCurve();
  // glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);

  // /* Calculate sensitivity curve with normalization and energy calibration errors,
  //  * as defined in the AEDL files */
  // InitOutput(MYFILE2,"Format: Running time   Log(10,s22th13) sens. \n");
  // ComputeSensitivityCurve();

/***************************************************************************
 *                                P L O T T I N G                          *
 ***************************************************************************/

// Notes
// Plot to compare Systematics on and off

cout << "Starting testing stuff \n";



cout << "Finished testing stuff \n";
userConfirm();

double plot_data_statvsys_a[2][tSteps];
double plot_data_statvsys_b[2][tSteps];

double plot_data_statvsysSPECTRAL_a[2][tSteps];
double plot_data_statvsysSPECTRAL_b[2][tSteps];

double plot_data_statvsysCOV_a[2][tSteps];
double plot_data_statvsysCOV_b[2][tSteps];


//ComputeSensitivityCurve();

  // printf("1Before curve \n");
  // double plot_data[3][max_index];
  // ComputeSensitivityCurve2(plot_data);
  // printf("1After curve");


printf("Starting curve calculations \n");
glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiNoSysSpectrum", sys_errors);
glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiNoSysSpectrum", sys_errors);
cout <<"after break\n";
ComputeSensitivityCurve3(plot_data_statvsys_a,1);
ComputeSensitivityCurve3(plot_data_statvsys_b,0);


  glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiZero", sys_errors);
  glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiZero", sys_errors);
  ComputeSensitivityCurve3(plot_data_statvsysSPECTRAL_a,1);
  ComputeSensitivityCurve3(plot_data_statvsysSPECTRAL_b,0);


  //glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiDCNormCOV", sys_errors);
  //glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiDCNormCOV", sys_errors);
  cout << "testing my function \n";

  generate_covariance(covariance_matrix_1);

  glbSetChiFunction(EXP_FAR, 0, GLB_ON, "chiCOV", sys_errors);
  glbSetChiFunction(EXP_NEAR, 0, GLB_ON, "chiCOV", sys_errors);

  
  // Moved to before function call
  // TMatrixD covariance_matrix_1(200,200);
  // generate_covariance(covariance_matrix_1);
  
  ComputeSensitivityCurve3(plot_data_statvsysCOV_a,1);
  ComputeSensitivityCurve3(plot_data_statvsysCOV_b,0);





printf("Curve calculations complete \n");
printf("Commencing potting \n");

  
  userConfirm();
  doPlotROOT(plot_data_statvsys_a,plot_data_statvsys_b,0,"statvsysSpectrum");

  userConfirm();
  doPlotROOT(plot_data_statvsysCOV_a,plot_data_statvsysCOV_b,0,"statvCOV");

  // for (int g = 0;g<tSteps; g++){
  //   printf("Data a is: index: %f and value: %f \n",plot_data_statvsys_a[0][g],plot_data_statvsys_a[1][g]);
  //   printf("Data b is: index: %f and value: %f \n",plot_data_statvsys_b[0][g],plot_data_statvsys_b[1][g]);
  // }
  // plot_data_statvsys_a[1][0] = plot_data_statvsys_a[1][1];
  // plot_data_statvsys_b[1][0] = plot_data_statvsys_b[1][1];

  // auto myCanvas = new TCanvas("test");
  // myCanvas->SetGrid();
  // auto spa = new TGraph(tSteps,plot_data_statvsys_a[0],plot_data_statvsys_a[1]);
  // auto spb = new TGraph(tSteps,plot_data_statvsys_b[0],plot_data_statvsys_b[1]);

  // //TGraphErrors* sp = new TGraphErrors(g->GetN(),pow(g->GetY(),10.0),g->GetX());

  // spa->SetMarkerStyle(4);
  // spa->SetMarkerColor(4);
  // spb->SetMarkerStyle(3);
  // spb->SetMarkerColor(3);
  
  // auto mg = new TMultiGraph();

  // mg->SetTitle("AB: Log Plot of sensitivity curve systematics vs stats only.");
  // mg->Add(spa);
  // mg->Add(spb);
  // mg->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  // mg->GetYaxis()->SetTitle("sin(2*theta_13)^2 sensitiivity");
  // mg->GetXaxis()->CenterTitle(true);
  // mg->GetYaxis()->CenterTitle(true);
  // // sp->SetTitle("A: Log Plot of sensitivity curve for Initial Simulation \n For CHIPS .glb");
  // // sp->GetXaxis()->SetTitle("Integrated detector luminosity GW t years");
  // // sp->GetYaxis()->SetTitle("sin(2*theta_23)^2 sensitiivity");
  // // sp->SetMarkerStyle(4);
  // // sp->SetMarkerColor(4);

  // gPad->SetLogx();
  // gPad->Update();
  // mg->Draw("APL");

  // TLegend* legend = new TLegend();
  // legend->SetHeader("Legend Title");
  // legend->AddEntry(spa,"series A: With Systematics","lp");
  // legend->AddEntry(spb,"series B: Statistics only","lp");
  // legend->Draw();

  // myCanvas->Update();
  
  // myCanvas->Modified();
  // myCanvas->Update();





/***************OLD**/

  // printf("Before curve \n");
  // double plot_data[3][max_index];
  // ComputeSensitivityCurve2(plot_data);
  // printf("After curve");



  // auto myCanvas = new TCanvas();

  // TH2* h2 = new TH2D("Experiment viability", "Surface plot of Chi Sensitivity Surface", 8, 200, 900, tSteps, plot_data[1][0], plot_data[1][29]);

  // myCanvas->SetGrid();
  // h2->GetXaxis()->SetTitle("Baseline km");
  // h2->GetYaxis()->SetTitle("Integrated detector luminosity GW t years");
  // h2->GetZaxis()->SetTitle("Chi Squared Sensitivity of 23 mixing angle");
  
  
  // for(int q = 1; q < tSteps*8; ++q) {
  //   //printf("Data inspection: index: %d, pow chi: %f, baseline: %f, exposure: %f \n",q,pow(10.0,plot_data[2][q]),plot_data[0][q],plot_data[1][q]);
  
  
  // h2->Fill(plot_data[0][q],plot_data[1][q],pow(10.0,plot_data[2][q]));
  // }

  // h2->Draw("SURF");
  // gPad->SetLogz();
  // gPad->Update();
  // TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
  // myCanvas->Update();
  // //myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/HistSimulation/Plotting/CHIPShist-First.pdf");

  app->Run();




  /* Clean up */
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  gsl_root_fsolver_free(s);

  return 0;
}