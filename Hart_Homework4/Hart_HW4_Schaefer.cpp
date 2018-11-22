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
  DATA_INTEGER(Nproj); // Number of projection years to run, set to 0 if you don't want projections, if Nproj > 0, you need to estimate logF_proj
  DATA_SCALAR(catch_proj); // Log F used in projection period, set to if Nproj = 0 this will not be used
  
  
  //DATA_INTEGER(data1); 
  //DATA_IVECTOR(data2); // Vector of integers ??
  // DATA_VECTOR(data_object_name); // Vector of data
  
  
  ///// Parameter section /////
  PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
  PARAMETER(Bzero); // Bzero = K = biomass at equilibrium
  PARAMETER(logr_growth); // Steepness parameter
  PARAMETER_VECTOR(logF_y); // Vector of fishing mortalities for years 2 on (1972 - 2000, the years for which catch > 0), length of 29
  
  // Retransform variables so not in log space 
  // Type local_variable1 = exp(log_variable1); 
  
  // Local variables
  //Type local_variable2; // single value variable which is NOT an integer
  Type r_growth = exp(logr_growth);
  
  //vector<Type> local_vector(5); // vector of length 5
  //matrix<Type> local_matrix(3,4); // 3X4col matrix
  vector<Type> biomass(catch_obs.size()+Nproj); // Biomass storage vector
  //vector<Type> recruitment(catch_obs.size()); // Recruitment storage vector
  vector<Type> catch_pred(catch_obs.size()+Nproj); // Catch observation storage vector
  Type rNLL = 0; // Initial value for Rzero
  Type BzeroNLL = 0; 
  
  Type obj_fun = 0; // NegativeLogLikelihood initialized at zero
  
  
  ///// Calculate parameters (could also do in R and pass in as DATA) /////
  
  
  ///// Code which contributes to objective function /////
  // Generally make a prediction and compare that prediction to data (likelihood)
  // minimize comparison so you pick parameter values to best predict data
  // obj_fun += likelihood_component1 + likelihood_component2...
  //obj_fun += dummy*dummy; // dummy objective function used to debug (make sure code compiles) without estimating all parameters
  
  // Equilibrium starting conditions
  //Rzero = (Bzero - (1.0 + rho)*exp(-1.0*M_dat)*Bzero + rho*exp(-2.0*M_dat)*Bzero)/(1.0 - rho*w_dat*exp(-1.0*M_dat)); // !!!! check to make sure this isn't zero, initialized at 0 but this line isn't being calculated
  //recruitment(iyear) = Rzero;
  //biomass(iyear) = Bzero;
  //catch_pred(iyear) = 0;
  
  // catch_obs in biomass equatio should be catch_pred
  
  for(int iyear=0; iyear<catch_obs.size()-1+Nproj; iyear++){
    if (iyear == 0){ 
      // 1970
      biomass(iyear) = Bzero; 
      catch_pred(iyear) = 0;
      // project 1971
      biomass(iyear+1) = biomass(iyear) + r_growth*biomass(iyear)*(1.0 - (biomass(iyear)/Bzero)) - 0; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear >= catch_obs.size()-1){ // If you enter the projection period
      biomass(iyear+1) = biomass(iyear) + r_growth*biomass(iyear)*(1.0 - (biomass(iyear)/Bzero)) - catch_pred(iyear); 
      catch_pred(iyear+1) = catch_proj;
    } else { // project 1972 - 2000
      biomass(iyear+1) = biomass(iyear) + r_growth*biomass(iyear)*(1.0 - (biomass(iyear)/Bzero)) - catch_pred(iyear); 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } 
    
    // If you are not in the projection period 
    if(iyear < catch_obs.size()-1){ 
      //// Catch likelihood component
      Type Catchdiff = 0;
      Catchdiff = catch_obs(iyear) - catch_pred(iyear);
      obj_fun -= dnorm(Catchdiff, Type(0), Type(0.2), TRUE); // TRUE means log, - so NLL
      // if last year (non-projection)
      if(iyear == catch_obs.size()-2){ 
        // Catch likelihood component
        // obj_fun += dnorm(catch_pred(iyear), catch_obs(iyear), Type (0.2)); // Want predicted catch to be very well fit to observed
        Type Catchdiff = 0;
        Catchdiff = catch_obs(iyear+1) - catch_pred(iyear+1);
        obj_fun -= dnorm(Catchdiff, Type(0), Type(0.2), TRUE); // TRUE means log, - so NLL
      }
      //// Biomass likelihood component
      Type Biodiff = 0;
      for (int isurvey=0; isurvey<survey_yr.size(); isurvey++){
        if(catch_yr(iyear) == survey_yr(isurvey)){
          Biodiff = survey_vals(isurvey) - biomass(iyear);
          obj_fun -= dnorm(Biodiff, Type(0), survey_SE(isurvey), TRUE); // second number is mean, maybe this shouldn't be zero??????????/
        }
        if(iyear == catch_obs.size()-2 && catch_yr(iyear+1) == survey_yr(isurvey)){ // if you are in the last year and the projected year matches a survey year
          Biodiff = survey_vals(isurvey) - biomass(iyear+1); // compare to projected biomass
          obj_fun -= dnorm(Biodiff, Type(0), survey_SE(isurvey), TRUE); // second number is mean, maybe this shouldn't be zero??????????/
        } // This last part of the likelihood was super wrong for a long time, pay attention to end year conditions
      }
    } 
    
  }
  
  ///// Prior likelihood components /////
  // log r (growth rate)
  rNLL = dnorm(logr_growth, Type(-1.2), Type(0.6), TRUE); 
  obj_fun -= rNLL;
  // uniform Bzero (carrying capacity, K)
  BzeroNLL = 1/(15000 - 500); // Each has equal chance of being chosen since bounded in R   //BzeroNLL = (Bzero - 500)/(15000 - 500); 
  obj_fun += BzeroNLL;
  
  
  ///// ADReport reports deviation /////
  // these are in sdreport file
  //ADREPORT(Bzero); // Bzero = K = biomass at equilibrium
  //ADREPORT(h_steep); // Steepness parameter
  //ADREPORT(logF_y); // Fishing mortality
  
  ///// Report /////
  REPORT(obj_fun);
  REPORT(rNLL);
  REPORT(BzeroNLL);
  REPORT(biomass);
  REPORT(catch_pred);
  
  
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
