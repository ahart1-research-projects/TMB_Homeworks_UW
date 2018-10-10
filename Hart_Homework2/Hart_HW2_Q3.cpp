// Amanda Hart
// Homework 2 FISH 559 Question1

#include <TMB.hpp>

// Optional function which square the value provided
template <class Type> Type square(Type x){return x*x;}

// Objective function which returns the value to minimize (often the negative log likelihood)
template<class Type>
  Type objective_function<Type>::operator() ()
{
  ///// Data section ///// (Years = data$Year, Males = data$Sex_M, Females = data$Sex_F, fishCount = fishCount)
  DATA_VECTOR(Years); // Vector of years
  DATA_VECTOR(Males); // Vector of male counts in each year
  DATA_VECTOR(Females); // Vector of female counts in each year
  DATA_VECTOR(fishCount); // Vector of total fish counts (Males+Females) in each year
  // DATA_IVECTOR(data_object_name); // Vector of data
  
  
  ///// Parameter section /////
  PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
  PARAMETER(logitp_ratio); // sex ratio parameter
  
  // Retransform variables
  Type p_ratio = exp(logitp_ratio)/(1+exp(logitp_ratio)); // -0.999 + Type(2)*0.999/(1+exp(LogitRho));
  
  // Local variables
  Type LL = 0.0;
  Type NLL = 0.0;// Always need objective function object
  
  ///// Code which contributes to objective function /////
    // Generally make a prediction and compare that prediction to data (likelihood)
  // minimize comparison so you pick parameter values to best predict data
  // obj_fun += likelihood_component1 + likelihood_component2...
  // LL += dummy*dummy; // dummy objective function used to debug (make sure code compiles) without estimating all parameters
  
  
  for (int iyr=0; iyr<61; iyr++) { // loop from 0 to iloop < 5 (0,1,2,3,4), int makes this a local loop variable
    LL += Males(iyr)*log(p_ratio) + Females(iyr)*log(1-p_ratio);
  }
  
  // Convert to negative log likelihood
  NLL = -1*LL;
  
  ///// ADReport reports deviation /////
    ADREPORT(p_ratio);
  
  ///// Return objective function /////
    return(NLL);
  
  ///// Advice if not compiling /////
    // Check that all lines end in ;
  // Check that * between multiplied objects in equations
  // Check indexing on objects in equations
  // Check indexing of storage objects (e.g. StorageVector = DataVector1[i]*DataVector2[i] won't work since this produces a single value but there is no indexing of StorageVector[i])
                                        // Check indexing of for() loops (indexing starts at 0)
                                        // Check spelling of objects in equations
                                        // Try commenting out parts of the line causing the error and gradually add them back in until it causes an error again to ID the source
  }
                                        