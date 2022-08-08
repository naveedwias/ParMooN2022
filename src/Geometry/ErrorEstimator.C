#include "ErrorEstimator.h"
#include "LocalAssembling.h"
#include "Domain.h"
#include "MooNMD_Io.h"
#include <cmath>

template <int d>
ErrorEstimator<d>::ErrorEstimator(const ParameterDatabase& param_db)
:  db(default_error_estimator_database())
{
  db.merge(param_db, false);
}

template<int d>
ParameterDatabase ErrorEstimator<d>::default_error_estimator_database()
{
  ParameterDatabase db("default algebraic flux correction database");

  db.add("estimator_type", 1u, "Chose which type of estimator to use" 
         "The estimator type depends on the derived class", 0u,10u);
  db.add("estimate_Dirichlet_cells", true,
         "Decide if Dirichlet cells take part in the error estimator.");
  
  auto la_db = LocalAssembling<d>::default_local_assembling_database();
  auto domain_db = TDomain::default_domain_parameters();
  
  db.add(domain_db["conforming_closure"]);
  db.add(la_db["space_discretization_type"]);
  
  
  return db;
}

template <int d>
void ErrorEstimator<d>::add_eta_K(const ErrorEstimator& estimator_inital,
                                  double alpha)
{ 
  if(this->currentCollection != estimator_inital.get_collection())
    ErrThrow("The collections don't match!");      
  maximal_local_error = -1.0;
  estimated_global_error = 0.0;
  for(int i=0; i < (int) eta_K.size(); i++)      
  {
    eta_K[i] = eta_K[i]*eta_K[i] + alpha*estimator_inital.get_eta_K(i)*estimator_inital.get_eta_K(i);
    maximal_local_error = std::max(eta_K[i], maximal_local_error);
    estimated_global_error += eta_K[i];
    eta_K[i] = std::sqrt(eta_K[i]);
  }
  estimated_global_error = std::sqrt(2*estimated_global_error);
  maximal_local_error = std::sqrt(maximal_local_error);
}


#ifdef __3D__
template class ErrorEstimator<3>;
#else
template class ErrorEstimator<2>;
#endif
