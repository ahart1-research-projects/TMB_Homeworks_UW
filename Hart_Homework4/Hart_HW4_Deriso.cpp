// Amanda Hart
// Homework 4 FISH 559

#include <TMB.hpp>

// Optional function which square the value provided
template <class Type> Type square(Type x){return x*x;}

// Objective function which returns the value to minimize (often the negative log likelihood)
template<class Type>
  Type objective_function<Type>::operator() ()
 {
  ///// Data section /////
  DATA_VECTOR(survey_vals); // Survey abundance
  DATA_VECTOR(survey_yr); // Survey years
  DATA_VECTOR(survey_SE); // Survey SE
  DATA_VECTOR(catch_obs); // Catch numbers, shorter than biomass, recruitment, and catch_pred (added 1 year to start of these timeseries) so offset by 2 in equations (for indexing reasons) rather than 1 observed in normal equation form (in reality refer to catches offset by 1 year)
  DATA_VECTOR(catch_yr); // Catch years
  DATA_SCALAR(rho); // Brody growth coefficient
  DATA_SCALAR(w_dat); // weight gain parameter
  DATA_SCALAR(M_dat); // Natural mortality rate
  
  //DATA_INTEGER(data1); 
  //DATA_IVECTOR(data2); // Vector of integers ??
  // DATA_VECTOR(data_object_name); // Vector of data
  
  
  
  ///// Parameter section /////
    PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
    PARAMETER(Bzero); // Bzero = K = biomass at equilibrium
    PARAMETER(h_steep); // Steepness parameter
    PARAMETER_VECTOR(logF_y); // Vector of fishing mortalities for years 2 on (1972 - 2000, the years for which catch > 0), length of 29
  
  // Retransform variables so not in log space 
  // Type local_variable1 = exp(log_variable1); 
  
  // Local variables
     //Type local_variable2; // single value variable which is NOT an integer
      //local_variable2 = 5; 
  
      //vector<Type> local_vector(5); // vector of length 5
      //matrix<Type> local_matrix(3,4); // 3X4col matrix
  vector<Type> biomass(catch_obs.size()); // Biomass storage vector
  vector<Type> recruitment(catch_obs.size()); // Recruitment storage vector
  vector<Type> catch_pred(catch_obs.size()); // Catch observation storage vector
  Type Rzero = 0; // Initial value for Rzero
  Type temp_logith = 0; 
  
  Type obj_fun = 0; // NegativeLogLikelihood initialized at zero
  
  
  ///// Calculate parameters (could also do in R and pass in as DATA) /////
    
  
  ///// Code which contributes to objective function /////
    // Generally make a prediction and compare that prediction to data (likelihood)
  // minimize comparison so you pick parameter values to best predict data
  // obj_fun += likelihood_component1 + likelihood_component2...
  //obj_fun += dummy*dummy; // dummy objective function used to debug (make sure code compiles) without estimating all parameters
  
  // Equilibrium starting conditions
  Rzero = (Bzero - (1.0 + rho)*exp(-1.0*M_dat)*Bzero + rho*exp(-2.0*M_dat)*Bzero)/(1.0 - rho*w_dat*exp(-1.0*M_dat)); // !!!! check to make sure this isn't zero, initialized at 0 but this line isn't being calculated
  //recruitment(iyear) = Rzero;
  //biomass(iyear) = Bzero;
  //catch_pred(iyear) = 0;
  
  // catch_obs in biomass equatio should be catch_pred
  
  for(int iyear=0; iyear<catch_obs.size()-1; iyear++){
    if (iyear == 0){ 
      // 1970
      recruitment(iyear) = Rzero;
      biomass(iyear) = Bzero; // (1.0 + rho)*exp(-M_dat)*Bzero - rho*exp(-2.0*M_dat)*Bzero - rho*w_dat*exp(-M_dat)*Rzero + Rzero;
      catch_pred(iyear) = 0;
      // project 1971
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1.0 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - 0) - rho*exp(-2.0*M_dat)*(1.0 - (0/biomass(iyear)))*(Bzero - 0) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (0/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear == 1){ // project 1972 
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear)-catch_pred(iyear)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear == 2){ // project 1973
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear)-catch_pred(iyear)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear == 3){ // project 1974
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear)-catch_pred(iyear)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else { // project 1975 - 2000 (since catch_obs.size()-1)
      recruitment(iyear+1) = (4*h_steep*Rzero*biomass(iyear-4)/Bzero)/((1-h_steep) + (5*h_steep-1)*biomass(iyear-4)/Bzero); 
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear)-catch_pred(iyear)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*recruitment(iyear) + recruitment(iyear+1); 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    }
    
    
    // Catch likelihood component
    // obj_fun += dnorm(catch_pred(iyear), catch_obs(iyear), Type (0.2)); // Want predicted catch to be very well fit to observed
    Type Catchdiff = 0;
    Catchdiff = catch_obs(iyear) - catch_pred(iyear);
    obj_fun -= dnorm(Catchdiff, Type(0), Type(0.2), TRUE); // TRUE means log, - so NLL
    
    if(iyear == catch_obs.size()-2){ // if last year
      // Catch likelihood component
      // obj_fun += dnorm(catch_pred(iyear), catch_obs(iyear), Type (0.2)); // Want predicted catch to be very well fit to observed
      Type Catchdiff = 0;
      Catchdiff = catch_obs(iyear+1) - catch_pred(iyear+1);
      obj_fun -= dnorm(Catchdiff, Type(0), Type(0.2), TRUE); // TRUE means log, - so NLL
    }
    
    // Biomass likelihood component
    Type Biodiff = 0;
    for (int isurvey=0; isurvey<survey_yr.size(); isurvey++){
      if(catch_yr(iyear) == survey_yr(isurvey)){
        Biodiff = survey_vals(isurvey) - biomass(iyear);
        obj_fun -= dnorm(Biodiff, Type(0), survey_SE(isurvey), TRUE); // second number is mean, maybe this shouldn't be zero??????????/
      }
    } 
    
  }
  
  ///// Prior likelihood components /////
  // logit h (steepness)
  temp_logith = logit((h_steep - 0.2)/0.8);
  obj_fun -= dnorm(temp_logith, Type(0.51), Type(2),TRUE); // True logs tish
  
  
  
  
  
  ///// ADReport reports deviation /////
  // these are in sdreport file
  //ADREPORT(Bzero); // Bzero = K = biomass at equilibrium
  //ADREPORT(h_steep); // Steepness parameter
  //ADREPORT(logF_y); // Fishing mortality
  
  ///// Report /////
  REPORT(recruitment); // Report variable or parameter value
  REPORT(biomass);
  REPORT(catch_pred);
  REPORT(catch_obs);
  REPORT(obj_fun);
  REPORT(rho);
  REPORT(w_dat);
  REPORT(M_dat);
  
  
  ///// Return objective function /////
  return(obj_fun);
  
  ///// Advice if not compiling /////
    // Check that all lines end in ;
  // Check that * between multiplied objects in equations
  // Check indexing on objects in equations
  // Check indexing of storage objects (e.g. StorageVector = DataVector1[i]*DataVector2[i] won't work since this produces a single value but there is no indexing of StorageVector[i])
                                        // Check indexing of for() loops (indexing starts at 0)
                                        // Check spelling of objects in equations
                                        // Try commenting out parts of the line causing the error and gradually add them back in until it causes an error again to ID the source
 }
