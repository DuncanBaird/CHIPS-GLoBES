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



void ComputeSensitivityCurve3(double plot_data[][tSteps])
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

    /* Determine sensvoid ComputeSensitivityCurve3(double plot_data[][tSteps])
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

    plot_data[0][j] = t;
    plot_data[1][j] = x; 
  }
}itivity to sin^2(2*th13) by  using the GSL Brent-Dekker
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

    plot_data[0][j] = t;
    plot_data[1][j] = x; 
  }
}