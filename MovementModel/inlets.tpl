// ----------------------------------------------------------------------------
// Inlet movement model
// 
// Authors: Michelle Jones and Samuel Johnson
// 
// Based on the model of McGarvey and Feenstra (2002)
// ----------------------------------------------------------------------------
// The following is an implementation of the McGarvey and Feenstra (2002) model
// for estimating movement probabilities from tag recovery data. It doesn't
// require tag releases, though these are included in case we wish to use the
// code as a simulator. 
// 
// There is code in here from an aborted attempt to introduce years at liberty
// which will take longer than I have now. That's a job for a flowchart
// and some equations.
// Note that by setting T = 1 the model collapses
// to a time averaged model, and t represents time at liberty rather than a
// time step.
// ----------------------------------------------------------------------------


DATA_SECTION
  // Read in simulated fishery and survey data from default input file
  init_int I;                 // Total Release areas
  init_int J;                 // total recovery areas
  // init_int T;              // Total number of time steps (years)
  init_vector r(1,I);         // tag releases (area x time)
  init_matrix R(1,I,1,J);     // observed tag recoveries (area x area x time)
  // init_matrix E(1,J,1,T);  // raw fishing effort by area, year
  // init_vector A(1,J,1,T);  // area fished in square km, by area, year
  init_number dumm;            // Check for correct data structure

  LOC_CALCS
    if (dumm != 999)
    { // Break out if the check is not correct
      cout<<"Error reading data.\n Fix it. \n dumm = " << dumm <<endl;
      ad_exit(1);
    }


PARAMETER_SECTION
  // Estimated parameters-----------------------------------------------------
  init_matrix logitTheta(1,I,1,J);  // array of movement log probabilities
  init_number logM(-1);               // log of natural mortality
  init_number logFbar(-1);            // log of average fishing mortality
  init_vector alpha(1,J,-1);         // Parameters of the Dirichlet prior
  init_number d(-1);                  // Sample size for mn jumping distribution
  init_vector Z(1,J,-1);             // Vector of area total mortalities
  init_vector F(1,J,-1);             // Fishing mortality by area
  //--------------------------------------------------------------------------

  // Derived Parameters-------------------------------------------------------
  sdreport_matrix Theta(1,I,1,J);          // Movement probabilities
  number M;                       // natural mortality
  number Fbar;                    // fishing mortality
  // vector F(1,J);                  // Fishing mortality by area
  // vector Z(1,J);                  // total mortality by area
  matrix Rhat(1,I,1,J);           // Predicted recoveries by area
  matrix pRhat(1,I,1,J);          // predicted prob of recovery by area
  //--------------------------------------------------------------------------  
  
  // Calculated pdfs----------------------------------------------------------
  vector mnLogLike(1,I);     // multinomial likelihood of each release area
  vector dirLogPrior(1,I);   // Dirichlet prior for each movement prob by rel area
  vector logPostFun(1,I);    // posterior pdf by release area and YAL
  //objective function value - what value to minimize?
  objective_function_value nll;
  //--------------------------------------------------------------------------

  // penalty variable---------------------------------------------------------
  vector pen(1,I);      // Variable to hold penalty for non-Markovian matrix
  //--------------------------------------------------------------------------
  
PROCEDURE_SECTION
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

FUNCTION predRecoveries
  // Loop to predict the number of recoveries in area j from release area i
  dvar_vector fishMortProp=elem_div(F,Z);
  dvar_vector mortAve=1.-mfexp(-Z);
  dvar_vector mort = elem_prod(mortAve,fishMortProp);
  for (int i=1; i<=I; i++)
  {  
    Rhat(i)=elem_prod(Theta(i),mort);
    // cout << "Rhat[" << i << "] = " << Rhat(i) << endl;
  }

FUNCTION predRecProb
  // Loop to turn recoveries in Rhat into proportions 
  for (int i=1; i<=I; i++)
  { 
    pRhat(i) = Rhat(i)/sum(Rhat(i));
  }

FUNCTION logPriorVector
  // Loop to compute prior for each row of the parameter matrix Theta (inlet)
  dvariable coeffNum = gammln(sum(alpha));
  dvariable coeffDenom = sum(gammln(alpha));
  for (int i=1; i<=I; i++)
  {  
    dvar_vector paramProd = elem_prod(alpha,log(Theta[i]));
    dirLogPrior[i] = coeffNum - coeffDenom + sum(paramProd);
  }

FUNCTION logLikelihoodVector
  // Loop to compute the likelihood function for each row of data R (inlet)
  for (int i=1; i<=I; i++)
  {
    dvar_vector logRhat = log(pRhat[i]);
    dvar_vector logPow = elem_prod(R[i],logRhat);
    // cout << "logRhat = " << logRhat << endl;
    // cout << "logPow = " << logPow << endl;
    mnLogLike[i] = sum(logPow);
  }

FUNCTION logPostVector
  // Computation of the vector of posterior probabilities for each inlet
  logPostFun = dirLogPrior + mnLogLike;

REPORT_SECTION
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