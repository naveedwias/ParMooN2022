/**
 * @brief test if derived classes of the system classes are possible
 * 
 * This essentially only tries to compile code which would normally be in a
 * separate repository where ParMooN is used as a library. Here we want to make
 * sure that this functionality is not accidentally broken. This test does not
 * run meaningful code, but if it compiles it can be considered successful.
 * 
 * We create a derived class whose base class part is copy-constructed. The base
 * class' copy constructor is protected, so that the other tests do not test it.
 */
#ifdef __2D__
constexpr int dim = 2;
#else
constexpr int dim = 3;
#endif

#include "Database.h"
#include "ConvectionDiffusion.h"
#include "NavierStokes.h"
#include "TimeNavierStokes.h"
#include "ParMooN.h"

#ifdef _MPI
#include <mpi.h>
#endif

/* ************************************************************************* */
template <int d>
class CD_derived : public ConvectionDiffusion<d>
{
  public:
    CD_derived(const ConvectionDiffusion<d>& cd)
     : ConvectionDiffusion<d>(cd)
    {
      Output::print("created CD_derived<", d, "> object");
    }
  
    void do_all()
    {
      this->assemble();
      this->solve();
      this->output();
    }
};

/* ************************************************************************* */
template <int d>
class NSE_derived : public NavierStokes<d>
{
  public:
    NSE_derived(const NavierStokes<d>& cd)
     : NavierStokes<d>(cd)
    {
      Output::print("created NSE_derived<", d, "> object");
    }
  
    void do_all()
    {
      this->assemble_linear_terms();
      this->solve();
      this->output();
    }
};

/* ************************************************************************* */
template <int d>
class TNSE_derived : public TimeNavierStokes<d>
{
  public:
    TNSE_derived(const TimeNavierStokes<d>& cd)
     : TimeNavierStokes<d>(cd)
    {
      Output::print("created TNSE_derived<", d, "> object");
    }
  
    void do_all()
    {
      this->assemble_initial_time();
      this->assemble_matrices_rhs(0);
      this->solve();
      this->output();
    }
};

// =======================================================================
// main program
// =======================================================================
int main(
#ifdef _MPI
  int argc, char* argv[]
#endif
)
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  Output::setVerbosity(2);
  
  db.add("refinement_n_initial_steps", 1u, "", 0u, 3u);
  TDomain domain(db);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  ConvectionDiffusion<dim> cd(domain, db);
  CD_derived<dim> cd_derived(cd);
  cd_derived.do_all();
  
  NavierStokes<dim> nse(domain, db);
  NSE_derived<dim> nse_derived(nse);
  nse_derived.do_all();
  
  TDatabase::ParamDB->PRESSURE_SPACE = -4711; // reset to a valid value
  TimeNavierStokes<dim> tnse(domain, db);
  TNSE_derived<dim> tnse_derived(tnse);
  tnse_derived.do_all();
  
  parmoon::parmoon_finalize();
}
