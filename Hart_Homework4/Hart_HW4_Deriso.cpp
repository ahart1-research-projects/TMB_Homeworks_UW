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
  DATA_INTEGER(Nproj); // Number of projection years to run, set to 0 if you don't want projections, if Nproj > 0, you need to estimate logF_proj
  DATA_SCALAR(catch_proj); // Log F used in projection period, set to if Nproj = 0 this will not be used
  
  ///// Parameter section /////
    PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
    PARAMETER(Bzero); // Bzero = K = biomass at equilibrium
    PARAMETER(h_steep); // Steepness parameter
    PARAMETER_VECTOR(logF_y); // Vector of fishing mortalities for years 2 on (1972 - 2000, the years for which catch > 0), length of 29
  
  // Local variables

  vector<Type> biomass(catch_obs.size()+Nproj); // Biomass storage vector
  vector<Type> recruitment(catch_obs.size()+Nproj); // Recruitment storage vector
  vector<Type> catch_pred(catch_obs.size()+Nproj); // Catch observation storage vector
  Type Rzero = 0; // Initial value for Rzero
  Type temp_logith = 0; 
  
  Type obj_fun = 0; // NegativeLogLikelihood initialized at zero

  // Equilibrium starting conditions
  Rzero = (Bzero - (1.0 + rho)*exp(-1.0*M_dat)*Bzero + rho*exp(-2.0*M_dat)*Bzero)/(1.0 - rho*w_dat*exp(-1.0*M_dat)); // !!!! check to make sure this isn't zero, initialized at 0 but this line isn't being calculated
  //recruitment(iyear) = Rzero;
  //biomass(iyear) = Bzero;
  //catch_pred(iyear) = 0;
  
  // catch_obs in biomass equatio should be catch_pred
  
  for(int iyear=0; iyear<(catch_obs.size()-1+Nproj); iyear++){ // projection not actually storing biomass/catch/recruitment????
//  for(int iyear=0; iyear<(catch_obs.size()-1); iyear++){
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
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear-1)-catch_pred(iyear-1)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear == 2){ // project 1973
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear-1)-catch_pred(iyear-1)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear == 3){ // project 1974
      recruitment(iyear+1) = Rzero;
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear-1)-catch_pred(iyear-1)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*Rzero + Rzero; 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    } else if (iyear >= catch_obs.size()-1){ // If you enter the projection period
      recruitment(iyear+1) = (4*h_steep*Rzero*biomass(iyear-4)/Bzero)/((1-h_steep) + (5*h_steep-1)*biomass(iyear-4)/Bzero); 
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear-1)-catch_pred(iyear-1)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*recruitment(iyear) + recruitment(iyear+1); 
      catch_pred(iyear+1) = catch_proj;
    } else { // project 1975 - 2000 (since catch_obs.size()-1)
      recruitment(iyear+1) = (4*h_steep*Rzero*biomass(iyear-4)/Bzero)/((1-h_steep) + (5*h_steep-1)*biomass(iyear-4)/Bzero); 
      biomass(iyear+1) = (1 + rho)*exp(-1.0*M_dat)*(biomass(iyear) - catch_pred(iyear)) - rho*exp(-2.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*(biomass(iyear-1)-catch_pred(iyear-1)) - rho*w_dat*exp(-1.0*M_dat)*(1.0 - (catch_pred(iyear)/biomass(iyear)))*recruitment(iyear) + recruitment(iyear+1); 
      catch_pred(iyear+1) = exp(logF_y(iyear))*biomass(iyear+1);
    }
    
    if(iyear < catch_obs.size()-1){ // If you are not in the projection period 
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
        }
      }
    } 
 
  }
  
  ///// Prior likelihood components /////
  // logit h (steepness)
  temp_logith = logit((h_steep - 0.2)/0.8);
  obj_fun -= dnorm(temp_logith, Type(0.51), Type(2),TRUE); // True logs tish
  
  
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
