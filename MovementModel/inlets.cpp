#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <inlets.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  I.allocate("I");
  J.allocate("J");
  r.allocate(1,I,"r");
  R.allocate(1,I,1,J,"R");
  dumm.allocate("dumm");
    if (dumm != 999)
    { // Break out if the check is not correct
      cout<<"Error reading data.\n Fix it. \n dumm = " << dumm <<endl;
      ad_exit(1);
    }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logitTheta.allocate(1,I,1,J,"logitTheta");
  logM.allocate(-1,"logM");
  logFbar.allocate(-1,"logFbar");
  alpha.allocate(1,J,-1,"alpha");
  d.allocate(-1,"d");
  Z.allocate(1,J,-1,"Z");
  F.allocate(1,J,-1,"F");
  Theta.allocate(1,I,1,J,"Theta");
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  Fbar.allocate("Fbar");
  #ifndef NO_AD_INITIALIZE
  Fbar.initialize();
  #endif
  Rhat.allocate(1,I,1,J,"Rhat");
  #ifndef NO_AD_INITIALIZE
    Rhat.initialize();
  #endif
  pRhat.allocate(1,I,1,J,"pRhat");
  #ifndef NO_AD_INITIALIZE
    pRhat.initialize();
  #endif
  mnLogLike.allocate(1,I,"mnLogLike");
  #ifndef NO_AD_INITIALIZE
    mnLogLike.initialize();
  #endif
  dirLogPrior.allocate(1,I,"dirLogPrior");
  #ifndef NO_AD_INITIALIZE
    dirLogPrior.initialize();
  #endif
  logPostFun.allocate(1,I,"logPostFun");
  #ifndef NO_AD_INITIALIZE
    logPostFun.initialize();
  #endif
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  pen.allocate(1,I,"pen");
  #ifndef NO_AD_INITIALIZE
    pen.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  nll =0.0;
  //// Convert parameters to natural scale
  // Transition probabilities
  Theta = elem_div(mfexp(logitTheta),(1 + mfexp(logitTheta)));
  // Natural mortality
  M = mfexp(logM);
  // Average fishing mortality
  Fbar = mfexp(logFbar);
  // Compute total mortality from natural and fishing mortality
  // Z = F + M;
  // Compute the matrix of predicted recoveries Rhat
  predRecoveries();
  // Turn Rhat into a matrix of probabilities by rescaling rows
  predRecProb();
  // Compute log prior pdf (Dirichlet) for each row (inlet)
  logPriorVector();
  // Compute log likelihood pdf (multinomial) for each row (inlet)
  logLikelihoodVector();
  // Compute log posterior pdf (Dirichlet...?) for each row (inlet)
  logPostVector();
  // Compute penalty for non-Markovian transition matrix Theta
  for (int i=1; i <= I; i++)
  {
    pen[i] = pow(log(sum(Theta[i])),2.);
  }
  // compute negative log likelihood of joint distribution of rows (inlets)
  nll = -1.*sum(logPostFun)+1000*sum(pen);
}

void model_parameters::predRecoveries(void)
{
  // Loop to predict the number of recoveries in area j from release area i
  dvar_vector fishMortProp=elem_div(F,Z);
  dvar_vector mortAve=1.-mfexp(-Z);
  dvar_vector mort = elem_prod(mortAve,fishMortProp);
  for (int i=1; i<=I; i++)
  {  
    Rhat(i)=elem_prod(Theta(i),mort);
    // cout << "Rhat[" << i << "] = " << Rhat(i) << endl;
  }
}

void model_parameters::predRecProb(void)
{
  // Loop to turn recoveries in Rhat into proportions 
  for (int i=1; i<=I; i++)
  { 
    pRhat(i) = Rhat(i)/sum(Rhat(i));
  }
}

void model_parameters::logPriorVector(void)
{
  // Loop to compute prior for each row of the parameter matrix Theta (inlet)
  dvariable coeffNum = gammln(sum(alpha));
  dvariable coeffDenom = sum(gammln(alpha));
  for (int i=1; i<=I; i++)
  {  
    dvar_vector paramProd = elem_prod(alpha,log(Theta[i]));
    dirLogPrior[i] = coeffNum - coeffDenom + sum(paramProd);
  }
}

void model_parameters::logLikelihoodVector(void)
{
  // Loop to compute the likelihood function for each row of data R (inlet)
  for (int i=1; i<=I; i++)
  {
    dvar_vector logRhat = log(pRhat[i]);
    dvar_vector logPow = elem_prod(R[i],logRhat);
    // cout << "logRhat = " << logRhat << endl;
    // cout << "logPow = " << logPow << endl;
    mnLogLike[i] = sum(logPow);
  }
}

void model_parameters::logPostVector(void)
{
  // Computation of the vector of posterior probabilities for each inlet
  logPostFun = dirLogPrior + mnLogLike;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "## Bayesian pdf values" << endl;
  report << "# priorVector" << endl;
  report << dirLogPrior << "\n" << endl;
  report << "# likelihoodVector" << endl;
  report << mnLogLike << "\n" << endl;
  report << "# postVector" << endl;
  report << logPostFun << "\n" << endl;
  report << "# PostValue" << endl;
  report << nll << "\n" << endl;
  report << "## Maximum Likelihood values" << endl;
  report << "# Theta" << endl;
  report << Theta << "\n" << endl;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
