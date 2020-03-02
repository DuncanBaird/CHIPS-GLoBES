//// Custom systematics functions - Duncan Baird

double chiDCNorm(int exp, int rule, int dim, double *x,
double *errors, void* user_data)
{
  double chi2 = 0.0;
  int i;
  ... /* Some code to calculate the chi2 */


  /* Here the systematics priors (penalties) are added: */
  for(i=0;i<dim;i++) chi2 += (x[i]*x[i])/(errors[i]*errors[i]);
return chi2;
}


int main(int argc,char*argv[])
{
  ...
  glbInit(agrv[0]);
  glbDefineChiFunction(&chiDCNorm,5,"chiDCNorm",NULL);
  ...
}