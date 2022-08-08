#include <cmath>
#include <string.h>

#include <LocalAssembling.h>

#include <Database.h>
#include <MainUtilities.h>
#include <MooNMD_Io.h>
#include "QuadFormula.h"

#include "BaseCell.h"

#include <ConvDiff.h>
#include "FEDatabase.h"
#ifdef __3D__
#include <FEFunction3D.h>
#include <ConvDiff3D.h>
#include <RefTrans3D.h>
#endif
#ifdef __2D__
#include <FEFunction2D.h>
#include <ConvDiff2D.h>
#include <RefTrans2D.h>
#endif
#include "NSE_local_assembling_routines.h"
#include "CD_local_assembling_routines.h"
#include "Time_NSE_local_assembling_routines.h"
#include <DarcyMixed.h>

#include <NonNewtonianViscosity.h>

#include <numeric> // std::iota
#include <algorithm>

// definitions - some compilers will be upset without these

template<int d> constexpr int LocalAssembling<d>::n_local_coefficients;
template<int d> constexpr int LocalAssembling<d>::n_total_coefficients;
template<int d> constexpr int LocalAssembling<d>::coeff_ref_coordinate_offset;
template<int d> constexpr int LocalAssembling<d>::coeff_ref_jacobian_offset;
template<int d> constexpr int LocalAssembling<d>::coeff_ref_inv_jacobian_offset;

std::ostream& operator<<(std::ostream& out, const LocalAssembling_type value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(LocalAssembling_type::p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(ConvDiff);
    PROCESS_VAL( Darcy );
    PROCESS_VAL(TCDStiffMassRhs);
    PROCESS_VAL(TCDDiffOnly);
    PROCESS_VAL(TCDMassOnly);
    PROCESS_VAL(TCDStiffOnly);
    PROCESS_VAL(TCDRhsOnly);
    PROCESS_VAL(NavierStokesAll);
    PROCESS_VAL(NavierStokesNL);
    PROCESS_VAL(NavierStokesLinear);
    PROCESS_VAL(TimeNavierStokesAll);
    PROCESS_VAL(TimeNavierStokesNL);
    PROCESS_VAL(TimeNavierStokesRhs);
    PROCESS_VAL(TimeNavierStokesMass);
    PROCESS_VAL(TimeNavierStokesExplNL);
    PROCESS_VAL(Custom);
    PROCESS_VAL(Poisson);
    default: s = "unknown LocalAssembling_type type"; break;
  }
#undef PROCESS_VAL
  return out << s;
}

template<int d>
typename Template_names<d>::MultiIndex_vector indices_up_to_order(int order);

template<>
Template_names<2>::MultiIndex_vector indices_up_to_order<2>(int order)
{
  switch(order)
  {
    case 0:
      return { MultiIndex2D::D00 };
      break;
    case 1:
      return { MultiIndex2D::D00, MultiIndex2D::D10, MultiIndex2D::D01 };
      break;
    case 2:
      return { MultiIndex2D::D00, MultiIndex2D::D10, MultiIndex2D::D01,
               MultiIndex2D::D20, MultiIndex2D::D11, MultiIndex2D::D02 };
      break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}

template<>
Template_names<3>::MultiIndex_vector indices_up_to_order<3>(int order)
{
  switch(order)
  {
    case 0:
      return { MultiIndex3D::D000 };
      break;
    case 1:
      return { MultiIndex3D::D000, MultiIndex3D::D100, MultiIndex3D::D010,
               MultiIndex3D::D001 };
      break;
    case 2:
      return { MultiIndex3D::D000,
               MultiIndex3D::D100, MultiIndex3D::D010, MultiIndex3D::D001,
               MultiIndex3D::D200, MultiIndex3D::D110, MultiIndex3D::D101,
               MultiIndex3D::D020, MultiIndex3D::D011, MultiIndex3D::D002 };
      break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}

template<int d>
typename Template_names<d>::MultiIndex_vector indices_of_order(int order);

template<>
Template_names<2>::MultiIndex_vector indices_of_order<2>(int order)
{
  switch(order)
  {
    case 0:
      return { MultiIndex2D::D00 };
      break;
    case 1:
      return { MultiIndex2D::D10, MultiIndex2D::D01 };
      break;
    case 2:
      return { MultiIndex2D::D20, MultiIndex2D::D11, MultiIndex2D::D02 };
      break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}

template<>
Template_names<3>::MultiIndex_vector indices_of_order<3>(int order)
{
  switch(order)
  {
    case 0:
      return { MultiIndex3D::D000 };
      break;
    case 1:
      return { MultiIndex3D::D100, MultiIndex3D::D010, MultiIndex3D::D001 };
      break;
    case 2:
      return { MultiIndex3D::D200, MultiIndex3D::D110, MultiIndex3D::D101,
               MultiIndex3D::D011, MultiIndex3D::D020, MultiIndex3D::D002 };
      break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}

template<int d>
typename Template_names<d>::MultiIndex_vector laplacian_indices();

template<>
Template_names<2>::MultiIndex_vector laplacian_indices<2>()
{
  return { MultiIndex2D::D20, MultiIndex2D::D02 };
}

template<>
Template_names<3>::MultiIndex_vector laplacian_indices<3>()
{
  return { MultiIndex3D::D200, MultiIndex3D::D020, MultiIndex3D::D002 };
}

template<int d>
LocalAssembling<d>::LocalAssembling(ParameterDatabase param_db,
                                    LocalAssembling_type t,
                                    std::vector<const FEFunction*> fefunctions,
                                    const CoeffFct& coeffs, int disctype)
 : db(default_local_assembling_database()), type(t),
   discretization_type(disctype), Coeffs(coeffs), fe_functions(fefunctions)
{
  Output::print<5>("Constructor of LocalAssembling3D: using type ", type);
  db.merge(param_db, true);
  Parameter disc_type{this->db["space_discretization_type"]};

  this->assembly_needs_reftrans = false;

  // the values below only matter if you need an existing finite element
  // function during your assembly. Change them in such a case
  this->N_Parameters = 0;
  this->N_ParamFct = 0;
  this->ParameterFct = {};
  this->N_FEValues = 0;
  this->N_PersistentDataPerPoint = 0;
  this->FEValue_FctIndex = {};
  this->FEValue_MultiIndex = {};
  this->BeginParameter = {};

  this->N_Spaces = 0; // is unused for built-in discretization types anyway

  // set all member variables according to type
  switch (this->type)
  {
    ///////////////////////////////////////////////////////////////////////////
    // stationary convection diffusion problems
    case LocalAssembling_type::ConvDiff:
    {
      this->N_Matrices = 1;
      this->RowSpace = { 0 };
      this->ColumnSpace = { 0 };
      this->N_Rhs = 1;
      this->RhsSpace = { 0 };
      this->N_Terms = d+1;
      this->Derivatives = indices_up_to_order<d>(1);
      this->Needs2ndDerivatives = new bool[1];
      this->Needs2ndDerivatives[0] = false;
      this->Manipulate = nullptr;
      this->local_assemblings_routines.push_back(BilinearAssembleGalerkin<d>);

      if ((disc_type.is("supg") || disc_type.is("gls")))
      {
        this->Derivatives = indices_up_to_order<d>(2);
        this->N_Terms = this->Derivatives.size();
        this->Needs2ndDerivatives[0] = true;
        if (disc_type.is("supg"))
        {
          this->local_assemblings_routines.push_back(
                    BilinearAssemble_SD<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
                    BilinearAssemble_GLS<d>);
        }
      }
      else if(!(disc_type.is("galerkin") || disc_type.is("dg")))
      {
        ErrThrow("currently the discretization type ", disc_type,
                 " is not supported by the class CD3D");
      }

      this->FESpaceNumber = std::vector<int>(this->Derivatives.size(), 0);
      break; // break for the type LocalAssembling3D_type::CD3D
    }

    case LocalAssembling_type::TCDStiffMassRhs:
    case LocalAssembling_type::TCDStiffRhs:
    case LocalAssembling_type::TCDMassOnly:
    case LocalAssembling_type::TCDStiffOnly:
    case LocalAssembling_type::TCDRhsOnly:
      this->set_parameters_for_tcd(type);
     break;

    ///////////////////////////////////////////////////////////////////////////
    case LocalAssembling_type::Darcy:
    {
      if(!disc_type.is("galerkin"))
      {
        ErrThrow("Darcy only supports a galerkin discretization currently, ",
                 disc_type);
      }
      // ( A B' )   ( 0 2 )
      // ( B C  )   ( 3 1 )
      this->N_Terms = 2*(d+1);
      //this->Derivatives = {D000, D000, D100, D010, D001, D100, D010, D001};
      auto fot = indices_up_to_order<d>(1);
      this->Derivatives = {fot[0]};
      this->Derivatives.insert(this->Derivatives.end(), fot.begin(), fot.end());
      this->Derivatives.insert(this->Derivatives.end(), fot.begin()+1,
                               fot.end());
      this->Needs2ndDerivatives = new bool[2];
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->N_Matrices = 4;
      this->RowSpace = { 0, 1, 0, 1 };
      this->ColumnSpace = { 0, 1, 1, 0 };
      this->N_Rhs = 2;
      this->RhsSpace = { 0, 1 };
      this->local_assemblings_routines.push_back(
              BilinearAssembleDarcyGalerkin<d>);
      if(d == 3)
        this->FESpaceNumber = { 0, 1, 0, 0, 0, 1, 1, 1 }; // 0: velocity, 1: pressure
      else
        this->FESpaceNumber = { 0, 1, 0, 0, 1, 1};
      this->Manipulate = nullptr;
      break;
    }

    ///////////////////////////////////////////////////////////////////////////
    // NSE3D: stationary Navier-Stokes problems
    case LocalAssembling_type :: NavierStokesAll:
    case LocalAssembling_type :: NavierStokesNL:
      if (!this->db["viscosity_mode"].is("newtonian"))
      {
        ErrThrow("Non-Newtonian viscosity is not currently implemented for "
          " stationary NSE.");
      }
      else if (disc_type.is("supg"))
      {
        this->set_parameters_for_nse_supg(type);
      }
      else if(disc_type.is("dg"))
      {
        this->set_parameters_for_nse_dg();
      }
      else
      {
        this->set_parameters_for_nse(type);
      }
      break;

    ////////////////////////////////////////////////////////////////////////////
    // TNSE3D: nonstationary Navier-Stokes problems
    case LocalAssembling_type::TimeNavierStokesAll:
    case LocalAssembling_type::TimeNavierStokesNL:
    case LocalAssembling_type::TimeNavierStokesRhs:
    case LocalAssembling_type::TimeNavierStokesMass:
    case LocalAssembling_type::NavierStokesLinear:
    case LocalAssembling_type::TimeNavierStokesExplNL:
      if (disc_type.is("residual_based_vms"))
      {
         this->set_parameters_for_tnse_residual_vms();
      }
      else if (disc_type.is("rbvms_time"))
      {
         this->set_parameters_for_tnse_rbvms_time();
      }
      else if (disc_type.is("vms_projection"))
      {
        if (!this->db["viscosity_mode"].is("newtonian"))
        {
          ErrThrow("Non-Newtonian viscosity is not currently implemented for "
            "TNSE with projection-based VMS.");
        }

        this->set_parameters_for_tnse_vms(type);
      }
      else if (disc_type.is("supg"))
      {
        if (!this->db["viscosity_mode"].is("newtonian"))
        {
          ErrThrow("Non-Newtonian viscosity is not currently implemented for "
            "TNSE with SUPG.");
        }

        this->set_parameters_for_tnse_supg(type);
      }
      else if (disc_type.is("dg"))
      {
        if (!this->db["viscosity_mode"].is("newtonian"))
        {
          ErrThrow("Non-Newtonian viscosity is not currently implemented for "
            "TNSE with DG.");
        }

        this->set_parameters_for_nse_dg();
      }
      else
      {
        this->set_parameters_for_tnse();
      }
      break;
    case LocalAssembling_type::Poisson:
      this->N_Matrices = 1;
      this->RowSpace = { 0 };
      this->ColumnSpace = { 0 };
      this->N_Rhs = 1;
      this->RhsSpace = { 0 };
      this->N_Terms = d+1;
      //this->Derivatives = { D000, D100, D010, D001 }; // or {D00, D10, D01}
      this->Derivatives = indices_up_to_order<d>(1);
      this->Needs2ndDerivatives = new bool[1];
      this->Needs2ndDerivatives[0] = false;
      this->Manipulate = nullptr;
      this->local_assemblings_routines.push_back(PoissonAssemble<d>);
      this->FESpaceNumber = std::vector<int>(this->Derivatives.size(), 0);
      break;

    default:
      ErrThrow("Unknown or unhandled LocalAssembling_type case. ", type);
  }

  AllOrigValues = new double** [N_Terms];
  OrigValues = new const double* [N_Terms];

  // some consistency checks
  if (Coeffs == nullptr)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }

  if (local_assemblings_routines.empty()
     || local_assemblings_routines[0] == nullptr)
  {
    ErrThrow("A local assembling routine was not set! ", this->type);
  }

  // find number of spaces
  int max = -1;

  for(int i = 0; i < N_Terms; i++)
  {
    int j = FESpaceNumber[i];

    if (j > max)
    {
      max = j;
    }
  }

  N_Spaces = max + 1;
}

//========================================================================
template<int d>
LocalAssembling<d>::LocalAssembling(
  int myN_Terms, MultiIndex_vector myDerivatives,
  std::vector<int> myFESpaceNumber, std::vector<int> myRowSpace,
  std::vector<int> myColumnSpace, std::vector<int> myRhsSpace,
  CoeffFct myCoeffs, std::vector<AssembleFctParam> myAssembleParam,
  ManipulateFct* myManipulate, int myN_Matrices, int myN_Rhs,
  int myN_ParamFct, std::vector<ParamFct*> myParameterFct,
  std::vector<int> myBeginParameter, int myN_Parameters,
  std::vector<const FEFunction*> myFEFunctions, int myN_FEValues,
  std::vector<int> myFEValue_FctIndex,
  MultiIndex_vector myFEValue_MultiIndex,
  int discretization_type_in)
 : db(default_local_assembling_database()), type{LocalAssembling_type::Custom},
   discretization_type{discretization_type_in}, N_Terms(myN_Terms),
   Derivatives(myDerivatives), FESpaceNumber(myFESpaceNumber),
   RowSpace(myRowSpace), ColumnSpace(myColumnSpace), RhsSpace(myRhsSpace),
   Coeffs(myCoeffs), local_assemblings_routines(myAssembleParam),
   Manipulate(myManipulate),
   N_Matrices(myN_Matrices), N_Rhs(myN_Rhs), N_ParamFct(myN_ParamFct),
   ParameterFct(myParameterFct), BeginParameter(myBeginParameter),
   N_Parameters(myN_Parameters), N_FEValues(myN_FEValues),
   fe_functions(myFEFunctions), FEValue_FctIndex(myFEValue_FctIndex),
   FEValue_MultiIndex(myFEValue_MultiIndex)
{
  this->N_PersistentDataPerPoint = 0;

  // Some data members get an extra treatment -
  // The auxiliary arrays (All)OrigValues are dynamically allocated with size
  // "N_Terms". "N_Spaces" is determined by finding the max in "FESpaceNumber"
  // (+1). "Needs2ndDerivative" is dynamically allocated to the size "N_Spaces"
  // and then filled according to the appearance of "D200", "D020", "D002",
  // "D110", "D101" or "D011" in "Derivatives".

  //Catch some things which might cause trouble.
  if ((int)myDerivatives.size() != N_Terms)
  {
    Output::print("Error: myDerivatives.size() != N_Terms.");
  }

  if ((int)myFESpaceNumber.size() != N_Terms)
  {
    Output::print("Error: myFESpaceNumber.size() != N_Terms.");
  }

  if ((int)myParameterFct.size() != N_ParamFct)
  {
    Output::print("Error: myParameterFct.size() != myN_ParamFct.");
  }

  if ((int)myBeginParameter.size() != N_ParamFct)
  {
    Output::print("Error: myBeginParameter.size() != myN_ParamFct.");
  }

  // Inform the world of what's going on.
  Output::print<5>("Constructor of LocalAssembling3D: using type ", type);

  assembly_needs_reftrans = false;

  // Dynamically allocate space for auxiliary arrays
  AllOrigValues = new double** [N_Terms];
  OrigValues = new const double* [N_Terms];

  // find number of spaces
  int max = -1;
  for (int i = 0; i < N_Terms; i++)
  {
    int j = FESpaceNumber[i];

    if (j > max)
    {
      max = j;
    }
  }

  N_Spaces = max + 1;

  //Fill the array Needs2ndDerivatives from the vector myNeeds2ndDerivatives
  Needs2ndDerivatives = new bool[N_Spaces];

  for (int i = 0; i < N_Spaces; i++)
  {
    Needs2ndDerivatives[i] = false;
  }

  for (int i = 0; i < N_Terms; i++)
  {
    auto alpha = Derivatives[i];
    auto sot = indices_up_to_order<d>(2);

    int j = FESpaceNumber[i];

    if (std::find_if(sot.begin() + d + 1, sot.end(),
      [alpha](decltype(alpha) mi)
      {
        return alpha == mi;
      })
      != sot.end())
    {
      Needs2ndDerivatives[j] = true;
    }
  }
}

template<int d>
LocalAssembling<d>::LocalAssembling(ParameterDatabase param_db)
: db(default_local_assembling_database()), type{LocalAssembling_type::Custom}
{
  db.merge(param_db, false);
}

//========================================================================
template<int d>
LocalAssembling<d>::~LocalAssembling()
{
  delete[] AllOrigValues;
  delete[] OrigValues;
  delete[] Needs2ndDerivatives;
}


//========================================================================
template<int d>
ParameterDatabase LocalAssembling<d>::default_local_assembling_database()
{
  ParameterDatabase db("default local assembling database");

  db.add("with_coriolis_force", false,
         "include the coriolis force for (Navier--)Stokes in 3D. This requires "
         "special local assemblings and the (pde-) coefficients of the example "
         "must include the coriolis force vector Omega.");

  db.add("laplace_type_deformation", false,
         "determine the way the laplacian is discretized.");

  db.add("nse_nonlinear_form", "convective",
         "Determine how the nonlinear term for Navier--Stokes is assembled. "
          "convective means ( (u.nabla) u, v), "
          "skew_symmetric means (1/2) [((u.nabla) u, v) - ((u.nabla) v, u)], "
          "rotational means ((nabla x u) x u, v), "
          "emac means (D(u)u + div(u)u, v)",
         {"convective", "skew_symmetric", "rotational", "divergence", "emac"});

  db.add("space_discretization_type", "galerkin",
         "The type of discretization. Note that not all types are possible for "
         "all problem classes.",
         {"galerkin", "supg", "upwind", "smagorinsky", "cip", "dg", "symm_gls",
          "nonsymm_gls", "pspg", "brezzi_pitkaeranta", "vms_projection",
          "vms_projection_expl", "local_projection", "local_projection_2_level",
          "residual_based_vms", "rbvms_time"});

  db.add("pspg_delta0", 0.1,
         "the stabilization parameter for pspg (Pressure Stabilization Petrov "
         "Galerkin) is delta0 * h^2 /nu, where h is a cell measure (e.g. "
         "diameter), nu is the inverse of the reynolds number, and delta0 is "
         " this parameter. This parameter is also used for Brezzi-Pitkaeranta",
         0., 10.);
  db.add("supg_delta0", 0.25,
         "the stabilization parameter for SUPG (Streamline-Upwind "
         "Petrov Galerkin) is delta0 * h^2, where h is a cell measure (e.g. "
         "diameter), delta0 is this parameter. ", 0., 10.);
  ///@todo add a parameter for the characteristic length L_0 (Brinkman case)
  db.add("graddiv_stab", 0.,
         "the stabilization parameter for Grad-Div is delta0 (nu + sigma L_0^2), "
   " where L is a characteristic length (taken equal to 1), nu is the "
   "inverse of the reynolds number,"
   " sigma is the inverse of permeability, and delta0 is this parameter.", 0., 10000.);

  db.add("gls_stab", 0.,
         "the stabilization parameter for GLS stabilization is "
   " delta0 h^2/(nu + sigma L_0^2), "
   " where L_0 is a characteristic length (taken equal to 1), nu is the "
   "inverse of the reynolds number,"
   " sigma is the inverse of permeability, and delta0 is this parameter.", 0., 10000.);

  db.add("L_0", 1.,
         "the parameter for relating Stokes with Darcy terms in Brinkman problem "
   " which is a characteristic length (by default equal to 1).", 0., 10000.);

  db.add("corner_stab", 0.,
         "the stabilization parameter needed for the Darcy limit of Brinkman in order"
         "to restrict normal jumps of the velocity across corners of the domain"
         "(by default equal to 0)", 0., 10000.);
  // the following 'local projection stabilization' (lps) parameters are not
  // used in this class, but are somewhat similar to other parameters here.
  db.add("lps_delta0", 0.1,
         "The stabilization parameter for local projection stabilization (lps) "
         "is delta0 * h^2 / nu, where h is a cell measure (e.g. diameter), nu "
         "is the inverse of the reynolds number, and delta0 is this parameter. "
         "Sometimes (in time-dependent problems) it is set to be "
         "delta0 * hK / nu.");
  db.add("lps_delta1", 0.1,
         "The stabilization parameter for local projection stabilization (lps) "
         "is delta0 * h^2 / nu, where h is a cell measure (e.g. diameter), nu "
         "is the inverse of the reynolds number, and delta0 is 'lps_delta0'. "
         "Sometimes (in time-dependent problems) it is set to be "
         "delta1 * hK / nu, where delta1 is this parameter.");

  db.add("lps_coeff_type", 0u, "Determine the way the local projection "
         "stabilization (lps) parameter is computed.", 0u, 5u);

  db.add("tensorial_diffusion", false, "This parameter enables the use of a tensorial "
           "diffusion coefficient instead or additional to a scalar (multiple of identity matrix). "
           "It is used, e.g., in the branch 'optimization' for describing longitudinal dispersion "
          "based on u*u^T. It can be scaled via coeff[d+3].");

  db.add("rbvms_param_mode", "g",
    "How to set the parameters \\tau_m, \\tau_c for the RBVMS discretization."
    "\ng: Using the reference map's Jacobian and the time step length."
    "\nh: \\delta_0 h_K^2 and \\delta_1.",
    { "g", "h", "codina" });

  db.add("rbvms_time_discretization", "backward_euler",
    "",
    { "backward_euler", "crank_nicolson" });

  db.add("rbvms_delta_0", 0.25,
    "Momentum base parameter for RBVMS in \"supg\" mode.",
    0.0, 1.e+10);

  db.add("rbvms_delta_1", 0.1,
    "Continuity (grad-div) parameter for RBVMS in \"supg\" mode.",
    0.0, 1.e+10);

  db.add("rbvms_tau_m_time", false,
    "",
    { true, false });

  db.add("rbvms_square_hk", true,
    "",
    { true, false });

  db.add("rbvms_c_inv", 1.0,
    "Inverse estimate constant for RBVMS in \"g\" mode. "
    "Practically relevant only with large viscosity.",
    0.0, 1.e+10);

  db.add("rbvms_tau_mul", 4.0,
    "Time step term multiplier for RBVMS in \"g\" mode. "
    "Should generally be 4.0 or zero.",
    0.0, 1.e+10);

  db.add("rbvms_explicit_time_derivative", false,
    "",
    { true, false });

  db.add("rbvms_momentum_pressure_coupling_B", false,
    "",
    { true, false });

  db.add("rbvms_momentum_pressure_coupling_C", false,
    "",
    { true, false });

  db.add("viscosity_mode", "newtonian",
    "Selects a viscosity model to use in Navier-Stokes simulations. "
    "Defaults to newtonian.",
    { "newtonian", "carreau_yasuda", "power_law", "cross", "herschel_bulkley" });

  db.add("viscosity_nu_0", 1.0,
    "In several non-Newtonian viscosity models, selects the "
    "limiting viscosity at zero shear rate."
    "\n- nu_0 is the limiting viscosity at zero shear rate in the "
    "Carreau-Yasuda and Cross fluid models."
    "\n- nu_0 is the minimum viscosity in the clamped power law model.",
    0.0, 1.e+10);

  db.add("viscosity_nu_infty", 1.0,
    "In several non-Newtonian viscosity models, selects the "
    "limiting viscosity at infinite shear rate."
    "\n- nu_infty is the limiting viscosity at infinite shear rate in the "
    "Carreau-Yasuda and Cross fluid models."
    "\n- nu_infty is the maximum viscosity in the clamped power law model.",
    0.0, 1.e+10);

  db.add("viscosity_tau_0", 0.0,
    "In the regularized Herschel-Bulkley model, selects the "
    "yield shear stress.",
    0.0, 1.e+10);

  db.add("viscosity_k", 1.0,
    "In several non-Newtonian viscosity models, selects the "
    "flow consistency index (or analogous parameter).",
    0.0, 1.e+10);

  db.add("viscosity_n", 1.0,
    "In several non-Newtonian viscosity models, selects the "
    "flow behavior index (or analogous parameter).",
    0.0, 1.e+10);

  db.add("viscosity_lambda", 1.0,
    "In the Carreau-Yasuda model, selects the relaxation time scale.",
    0.0, 1.e+10);

  db.add("viscosity_a", 2.0,
    "In the Carreau-Yasuda model, selects the Yasuda exponent. "
    "Defaults to a = 2 as in the Carreau model.",
    0.0, 1.e+10);

  db.add("viscosity_g_0", 1.0,
    "In the regularized Herschel-Bulkley model, selects the "
    "limiting shear rate. At shear rates below this threshold, viscosity is "
    "computed as for shear rate g_0.",
    0.0, 1.e+10);

  db.merge(ParameterDatabase::parmoon_default_database());

  return db;
}


//========================================================================
template<int d>
void LocalAssembling<d>::GetLocalForms(const TQuadFormula& qf,
                                       std::vector<const BaseFunctions*>& BaseFuncts,
                                       const TBaseCell *Cell, int cell_num,
                                       int N_Matrices, int N_Rhs,
                                       double ***LocMatrix, double **LocRhs,
                                       double factor)
{
  int N_Points = qf.GetN_QuadPoints();
  this->GetParameters(qf, cell_num);
  this->local_coefficients.resize(N_Points);

  if (N_PersistentDataPerPoint > 0)
  {
    if (persistent_data == nullptr)
    {
      ErrThrow("Tried to use a method that needs persistent pointwise data "
        "without providing storage!");
    }

    persistent_data->SetCurrentCell(cell_num, N_Points, N_PersistentDataPerPoint);
  }

  int i,j, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult;
  double *Coeff;
  const double *Param;
  const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

  // the following is only used in the call of Coeffs and Manipulate
  std::vector<const double *> parameters(N_Points);
  std::vector<double *> coefficients(N_Points);
  std::vector<int> N_BaseFuncts(BaseFuncts.size());
  std::transform(BaseFuncts.begin(), BaseFuncts.end(), N_BaseFuncts.begin(),
                 [](const BaseFunctions* bf){ return bf->GetDimension(); });

  for (int i = 0; i < N_Points; i++)
  {
    if (N_ParamFct == 0)
    {
      parameters[i] = nullptr;
    }
    else
    {
      parameters[i] = &this->parameter_functions_values[i * N_Parameters];
    }

    coefficients[i] = local_coefficients[i].data();
  }

  for (i = 0; i < N_Matrices; ++i)
  {
    CurrentMatrix = LocMatrix[i];
    N_Rows = BaseFuncts[RowSpace[i]]->GetDimension();

    N_Columns = BaseFuncts[ColumnSpace[i]]->GetDimension();

    for (j = 0; j < N_Rows; j++)
    {
      MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, sizeof(double)*N_Columns);
    }
  }

  for (i = 0; i < N_Rhs; ++i)
  {
    N_Rows = BaseFuncts[RhsSpace[i]]->GetDimension();
    memset(LocRhs[i], 0, sizeof(double)*N_Rows);
  }

  // *****************************************************
  // for 2Phase flow problems (Sashikumaar Ganesan)
  coefficients[0][0] = Cell->GetPhase_ID();
  coefficients[0][1] = Cell->GetRegionID();
  coefficients[0][2] = hK;
  // *****************************************************

  if (assembly_needs_reftrans)
  {
    FE_type element_type =
      TFESpace::get_element_for_shape(0, d)[Cell->GetShapeDesc()->GetType()];

    FiniteElement element(element_type);

    ReferenceTransformation_type ref_trans = element.GetRefTransID();
    BFRefElements ref_element = element.GetBaseFunct()->GetRefElement();

    if (TDatabase::ParamDB->USE_ISOPARAMETRIC
    && Cell->has_isoparametric_joint())
    {
      switch (ref_element)
      {
        case BFRefElements::BFUnitHexahedron:
          ref_trans = ReferenceTransformation_type::HexaIsoparametric;
          break;

        case BFRefElements::BFUnitTetrahedron:
          ref_trans = ReferenceTransformation_type::TetraIsoparametric;
          break;

        default:
          Output::root_warn("RefTrans", "Unknown isoparametric reference "
            "transformation for reference element ", ref_element);
          break;
      }
    }

#ifdef __3D__
    TRefTrans3D *F_K = FEDatabase::GetRefTrans3D(ref_trans);
#else
    TRefTrans2D *F_K = FEDatabase::GetRefTrans2D(ref_trans);
#endif

    F_K->SetCell(Cell);

    const double* xi = qf.get_xi();
    const double* eta = qf.get_eta();
    const double* zeta = qf.get_zeta();

    std::array<double, d * d> J;
    std::array<double, d * d> J_inv;

    for (int i = 0; i < N_Points; i++)
    {
      double* point_coeff = coefficients[i];

      auto p = qf.get_point(i);

      point_coeff[coeff_ref_coordinate_offset + 0] = xi[i];
      point_coeff[coeff_ref_coordinate_offset + 1] = eta[i];

      if (d == 3)
      {
        point_coeff[coeff_ref_coordinate_offset + 2] = zeta[i];
      }

      if (d == 2)
      {
#ifdef __2D__
        F_K->GetTransformationDerivatives(xi[i], eta[i], J.data());
        F_K->GetInverseTransformationDerivatives(xi[i], eta[i], J_inv.data());
#endif
      }
      else
      {
#ifdef __3D__
        F_K->GetTransformationDerivatives(xi[i], eta[i], zeta[i],
          J.data());

        F_K->GetInverseTransformationDerivatives(xi[i], eta[i], zeta[i],
          J_inv.data());
#endif
      }

      for (int j = 0; j < d * d; j++)
      {
        point_coeff[coeff_ref_jacobian_offset + j] = J[j];
        point_coeff[coeff_ref_inv_jacobian_offset + j] = J_inv[j];
      }
    }
  }

  if (Coeffs)
  {
    const double * X = qf.get_xi();
    const double * Y = qf.get_eta();
#ifdef __3D__
    const double * Z = qf.get_zeta();
    Coeffs(N_Points, X, Y, Z, parameters.data(), coefficients.data());
#else
    Coeffs(N_Points, X, Y, parameters.data(), coefficients.data());
#endif
  }

  if (Manipulate)
  {
    Manipulate(N_Points, coefficients.data(), parameters.data(), Cell);
  }

  for (i = 0; i < N_Terms; ++i)
  {
    AllOrigValues[i] = FEDatabase::GetOrigElementValues(
      *BaseFuncts[FESpaceNumber[i]], Derivatives[i]);
  }

  for (i = 0; i < N_Points; ++i)
  {
    if (N_PersistentDataPerPoint > 0)
    {
      persistent_data->SetCurrentPoint(i);
    }

    Mult = qf.get_weight(i) * factor;
    Coeff = local_coefficients[i].data();

    Param = parameters[i];

    for (j = 0; j < N_Terms; j++)
    {
      OrigValues[j] = AllOrigValues[j][i];
    }

    for (int k = 0; k < N_Rhs; ++k)
    {
      N_Rows = BaseFuncts[RhsSpace[k]]->GetDimension();
      for (int j = 0; j < N_Rows; j++)
      {
        if (LocRhs[k][j] != 0.0 && !std::isnormal(LocRhs[k][j]))
        {
          Output::root_warn<5>("LocalAssembling", "RHS ", k, ".", j,
            "point ", i, " prior: ", LocRhs[k][j]);
        }
      }
    }

    int m = 0;

    for (auto& lar: local_assemblings_routines)
    {
      lar(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts.data(), LocMatrix,
          LocRhs);

      for (int k = 0; k < N_Rhs; ++k)
      {
        N_Rows = BaseFuncts[RhsSpace[k]]->GetDimension();
        for (int j = 0; j < N_Rows; j++)
        {
          if (LocRhs[k][j] != 0.0 && !std::isnormal(LocRhs[k][j]))
          {
            Output::root_warn<5>("LocalAssembling", "RHS ", k, ".", j,
              "point ", i, " assembler ", m, ": ", LocRhs[k][j]);
          }
        }
      }
    }
  } // end loop over quadrature points
}

//========================================================================
template<int d>
void LocalAssembling<d>::GetParameters(const TQuadFormula& qf, int cellnum)
{
  if (N_ParamFct == 0)
  {
    return;
  }

  int n_points = qf.GetN_QuadPoints();

  if (N_PersistentDataPerPoint > 0)
  {
    if (persistent_data == nullptr)
    {
      ErrThrow("Tried to use a method that needs persistent pointwise data "
        "without providing storage!");
    }

    persistent_data->SetCurrentCell(cellnum, n_points, N_PersistentDataPerPoint);
  }

  this->parameter_functions_values.resize(n_points * this->N_Parameters, 0.0);

  const double *CurrValues;
  double *CurrOrigValues;
  const int *CurrIndex;
  std::vector<int> N_BaseFunct(N_FEValues);
  std::vector<const double *> Values(N_FEValues);
  std::vector<double **> orig_values(N_FEValues);
  std::vector<const int *> Index(N_FEValues);
  std::vector<int> base_vec_dims(N_FEValues);
  unsigned int param_length = d; // quad point has d entries


   // collect information
  for (int j = 0; j < N_FEValues; j++)
  {
    auto fefunction = fe_functions[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();

#ifdef __2D__
    auto fespace = fefunction->GetFESpace2D();
#else
    auto fespace = fefunction->GetFESpace3D();
#endif

    auto& fe = fespace->get_fe(cellnum);
    base_vec_dims[j] = fe.GetBaseFunct()->GetBaseVectDim();
    param_length += base_vec_dims[j];
    N_BaseFunct[j] = fe.GetN_DOF();
    orig_values[j] = FEDatabase::GetOrigElementValues(*fe.GetBaseFunct(),
                                                      FEValue_MultiIndex[j]);
    Index[j] = fespace->GetGlobalDOF(cellnum);
  }


  std::vector<double> Temp(param_length);

  // loop over all quadrature points
  for (int i = 0; i < n_points; i++)
  {
    if (N_PersistentDataPerPoint > 0)
    {
      persistent_data->SetCurrentPoint(i);
    }

    // first d parameters are the coordinates
    auto p = qf.get_point(i);
    Temp[0] = p.x;
    Temp[1] = p.y;

    if (d == 3)
    {
      Temp[2] = p.z;
    }

    // loop to calculate all FE values

    for (int k = d, j = 0; j < N_FEValues; j++)
    {
      double s = 0;
      int n = N_BaseFunct[j];

      CurrValues = Values[j];
      CurrOrigValues = orig_values[j][i];
      CurrIndex = Index[j];

      int baseVecDim = base_vec_dims[j];
      for (int m = 0; m < baseVecDim; ++m)
      {
        for (int l = 0; l < n; l++)
        {
          s += CurrValues[CurrIndex[l]] * CurrOrigValues[m * n + l];
        }

        Temp[k + m] = s;
      } // endfor m

      k += baseVecDim;
    }

    // loop to calculate all parameters
    for (int j = 0; j < N_ParamFct; j++)
    {
      //currparam = param + BeginParameter[j];
      double * currparam = &this->parameter_functions_values[
        i * N_Parameters + this->BeginParameter[j]];

      ParameterFct[j](Temp.data(), currparam);
    } // endfor j
  } // endfor i
}

//========================================================================
template<int d>
void LocalAssembling<d>::GetParameters(const TQuadFormula& qf, int cellnum, 
                                       double **Parameters)
{
  if (N_ParamFct == 0)
  {
    return;
  }

  int n_points = qf.GetN_QuadPoints();

  const double *CurrValues;
  double *CurrOrigValues;
  const int *CurrIndex;
  std::vector<int> N_BaseFunct(N_FEValues);
  std::vector<const double *> Values(N_FEValues);
  std::vector<double **> orig_values(N_FEValues);
  std::vector<const int *> Index(N_FEValues);
  std::vector<int> base_vec_dims(N_FEValues);
  unsigned int param_length = d; // quad point has d entries


   // collect information
  for (int j = 0; j < N_FEValues; j++)
  {
    auto fefunction = fe_functions[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();

#ifdef __2D__
    auto fespace = fefunction->GetFESpace2D();
#else
    auto fespace = fefunction->GetFESpace3D();
#endif

    auto& fe = fespace->get_fe(cellnum);
    base_vec_dims[j] = fe.GetBaseFunct()->GetBaseVectDim();
    param_length += base_vec_dims[j];
    N_BaseFunct[j] = fe.GetN_DOF();
    orig_values[j] = FEDatabase::GetOrigElementValues(*fe.GetBaseFunct(),
                                                      FEValue_MultiIndex[j]);
    Index[j] = fespace->GetGlobalDOF(cellnum);
  }


  std::vector<double> Temp(param_length);

  // loop over all quadrature points
  for (int i = 0; i < n_points; i++)
  {
    auto param = Parameters[i];

    // first d parameters are the coordinates
    auto p = qf.get_point(i);
    Temp[0] = p.x;
    Temp[1] = p.y;

    if (d == 3)
    {
      Temp[2] = p.z;
    }

    // loop to calculate all FE values

    for (int k = d, j = 0; j < N_FEValues; j++)
    {
      double s = 0;
      int n = N_BaseFunct[j];

      CurrValues = Values[j];
      CurrOrigValues = orig_values[j][i];
      CurrIndex = Index[j];

      int baseVecDim = base_vec_dims[j];
      for (int m = 0; m < baseVecDim; ++m)
      {
        for (int l = 0; l < n; l++)
        {
          s += CurrValues[CurrIndex[l]] * CurrOrigValues[m * n + l];
        }

        Temp[k + m] = s;
      } // endfor m

      k += baseVecDim;
    }

    // loop to calculate all parameters
    for (int j = 0; j < N_ParamFct; j++)
    {
      //currparam = param + BeginParameter[j];
      double * currparam = param + this->BeginParameter[j];

      ParameterFct[j](Temp.data(), currparam);
    } // endfor j
  } // endfor i
}
//========================================================================
template<int d>
void LocalAssembling<d>::add_parameter_function(int n_parameters,
                                                ParamFct* param_fct)
{
  N_ParamFct++;
  ParameterFct.push_back(param_fct);
  BeginParameter.push_back(N_Parameters);
  N_Parameters += n_parameters;
}

//========================================================================
template<int d>
void LocalAssembling<d>::add_fe_function(const FEFunction* f,
                                         MultiIndex_vector mi)
{
  fe_functions.push_back(f);
  auto last_fe_functions_index = fe_functions.size() - 1;
  auto n_mi = mi.size();
  N_FEValues += n_mi;

  bool second_derivatives = false;
  for (auto multi_index: mi)
  {
    FEValue_FctIndex.push_back(last_fe_functions_index);
    FEValue_MultiIndex.push_back(multi_index);
    if (int(multi_index) > d)
    {
      second_derivatives = true;
    }
  }

  bool *newSecondDer = new bool[N_Spaces + 1];

  std::copy(Needs2ndDerivatives,
            Needs2ndDerivatives + N_Spaces, newSecondDer);

  newSecondDer[N_Spaces] = second_derivatives;

  delete[] Needs2ndDerivatives;
  Needs2ndDerivatives = newSecondDer;
}

//========================================================================
template<int d>
void LocalAssembling<d>::add_derivatives(MultiIndex_vector derivatives,
                                         std::vector<int>  feSpaceNumber)
{
  auto n_derivatives = derivatives.size();
  auto n_feSpaceNumber = feSpaceNumber.size();

  if(n_derivatives!=n_feSpaceNumber)
  {
     ErrThrow("add_derivatives(): derivatives and feSpaceNumber do not have "
              "the same length (", n_derivatives , " and ", n_feSpaceNumber,
              " respectively)");
  }

  Derivatives.insert(Derivatives.end(), derivatives.begin(), derivatives.end());
  for(unsigned int i = 0; i < n_derivatives; ++i)
      FESpaceNumber.push_back(feSpaceNumber[i]);

  // adapt N_Terms, AllOrigValues and OrigValues
  N_Terms += n_derivatives;
  delete[] AllOrigValues;
  AllOrigValues = new double** [N_Terms];
  delete[] OrigValues;
  OrigValues = new const double* [N_Terms];
}

//========================================================================
template<int d>
void LocalAssembling<d>::replace_local_assembling(
  const std::vector<AssembleFctParam>& laf)
{
  local_assemblings_routines = laf;
}

//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_tcd(LocalAssembling_type type)
{
  this->N_Matrices = 2;
  this->RowSpace = { 0, 0 };
  this->ColumnSpace = { 0, 0 };
  this->N_Rhs = 1;
  this->RhsSpace = { 0 };
  this->N_Terms = d+1;
  //this->Derivatives = { D000, D100, D010, D001 }; // or {D00, D10, D01}
  this->Derivatives = indices_up_to_order<d>(1);
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->FESpaceNumber = std::vector<int>(d+1, 0);
  this->Manipulate = nullptr;

  Parameter tensorial_diffusion{this->db["tensorial_diffusion"]};
  Parameter disc_type{this->db["space_discretization_type"]};
  switch(type)
  {
    case LocalAssembling_type::TCDStiffMassRhs:
      // stiffness matrix and rhs
      this->local_assemblings_routines.push_back(TCDStiff<d>);
      // additional tensorial diffusion in stiffness matrix
      if(tensorial_diffusion.is(true))
       this->local_assemblings_routines.push_back(TCDStiff_TensorialDiffusionTerm<d>);
      // mass matrix
      this->local_assemblings_routines.push_back(TCDMass<d>);
      // rhs
      this->local_assemblings_routines.push_back(TCDRhs<d>);
      if(disc_type.is("supg"))
      {
        this->Derivatives = indices_up_to_order<d>(2);
        this->N_Terms = this->Derivatives.size();
        this->Needs2ndDerivatives[0] = true;
        this->Needs2ndDerivatives[1] = true;
        this->FESpaceNumber = std::vector<int>(N_Terms, 0);
        // stiffness matrix and rhs
        this->local_assemblings_routines.push_back(TCDStiffSUPG<d>);
        // mass matrix
        this->local_assemblings_routines.push_back(TCDMassSUPG<d>);
        // rhs
        this->local_assemblings_routines.push_back(TCDRhsSUPG<d>);

        if(tensorial_diffusion.is(true))
          ErrThrow("The use of a tensorial diffusion term is not yet "
                  "implemented in combination with SUPG.");
      }
      else if(!disc_type.is("galerkin"))
      {
        ErrThrow("currently the discretization type ", disc_type,
                 " is not supported by the class Time_CD2D");
      }
      break;
    case LocalAssembling_type::TCDStiffRhs:
      // stiff matrix, rhs
      // stiffness matrix and rhs
      this->local_assemblings_routines.push_back(TCDStiff<d>);

      // additional tensorial diffusion in stiffness matrix
      if(tensorial_diffusion.is(true))
       this->local_assemblings_routines.push_back(TCDStiff_TensorialDiffusionTerm<d>);

      // rhs
      this->local_assemblings_routines.push_back(TCDRhs<d>);

      if(disc_type.is("supg"))
      {
        // stiffness matrix
        this->local_assemblings_routines.push_back(TCDStiffSUPG<d>);
        // rhs
        this->local_assemblings_routines.push_back(TCDRhsSUPG<d>);
      }
      break;
    case LocalAssembling_type::TCDRhsOnly:
      N_Matrices = 0;
      this->local_assemblings_routines.push_back(TCDRhs<d>);
      if(disc_type.is("supg"))
      {
        this->local_assemblings_routines.push_back(TCDRhsSUPG<d>);
      }
      break;
    case LocalAssembling_type::TCDMassOnly:
      N_Rhs = 0;
      this->local_assemblings_routines.push_back(TCDMass<d>);
      break;
    case LocalAssembling_type::TCDStiffOnly:
      N_Rhs = 0;
      this->local_assemblings_routines.push_back(TCDStiff<d>);
      break;
    default:
      ErrThrow("unknown LocalAssembling_type ", this->type);
      break;
  }
}
//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_nse( LocalAssembling_type type)
{
  //bool with_coriolis = db["with_coriolis_force"];
  bool laplace_type_deformation = this->db["laplace_type_deformation"];
  std::string disc_type = this->db["space_discretization_type"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  bool galerkin = (disc_type == std::string("galerkin"));
  bool pspg = (disc_type == std::string("pspg"));
  bool symm_gls = (disc_type == std::string("symm_gls"));
  bool nonsymm_gls = (disc_type == std::string("nonsymm_gls"));
  bool brezzi_pitkaeranta = (disc_type == std::string("brezzi_pitkaeranta"));
  bool local_projection = (disc_type == std::string("local_projection"));
  int nstype = TDatabase::ParamDB->NSTYPE;
  int problem_type = db["problem_type"];

  if(laplace_type_deformation)
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("laplace_type_deformation is only supported for NSTYPE 3, 4, "
               "and 14");
    }
  }
  if(nonlin_form.is("rotational") || nonlin_form.is("emac"))
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("nse_nonlinear_form ", nonlin_form,
               " is only supported for NSTYPE 3, 4, and 14");
    }
  }
  if(!galerkin && !pspg && !symm_gls && !nonsymm_gls && !brezzi_pitkaeranta
     && !local_projection)
  {
    ErrThrow("unsupported space_discretization_type for NSE", d, "D: ",
             disc_type);
  }
  if((pspg || symm_gls || nonsymm_gls) && nstype != 14)
  {
    ErrThrow("for PSPG, symmetric GLS and non-symmetric GLS stabilization we "
             "need separate B and BT blocks as well as a C block, i.e., "
             "nstype 14");
  }
  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = d+1;
  this->RhsSpace = { 0, 0, 0, 1 };
  this->N_Terms = d+2;
  auto foi = indices_up_to_order<d>(1); // first_order_index
  this->Derivatives = indices_up_to_order<d>(0);
  this->Derivatives.insert(this->Derivatives.end(), foi.begin(), foi.end());
  // Derivatives = { D000, D000, D100, D010, D001 } or { D00, D00, D10, D01}
  this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure
  if(d == 3)
    this->FESpaceNumber.push_back(0);
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;
  this->N_Parameters = d;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = d;
  this->FEValue_FctIndex = std::vector<int>(d);
  std::iota(this->FEValue_FctIndex.begin(), this->FEValue_FctIndex.end(), 0);
  this->FEValue_MultiIndex = MultiIndex_vector(d, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  this->N_Matrices = (d+1)*(d+1); // 9, 16
  if(d==3)
  {
    this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
  }
  else
  {
    this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
  }
  if(laplace_type_deformation)
  {
    this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
    if(TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR && db["problem_type"].is(3))
    {
      this->Needs2ndDerivatives[0] = true;
      this->Needs2ndDerivatives[1] = true;
      this->local_assemblings_routines.push_back(OseenSingle<d>);
    }
  }
  else
  {
    this->Needs2ndDerivatives[0] = true;
    this->Needs2ndDerivatives[1] = true;
    if(TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR && db["problem_type"].is(3))
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGrad<d>);
      this->local_assemblings_routines.push_back(OseenSingle<d>);
    }
    else if(nstype == 1 || nstype == 2)
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGradSingle<d>);
    }
  }

  // add Darcy (resistance) term for Brinkman problem
  if (problem_type == 7)
  {
    if(nstype == 1 || nstype == 2)
    {
      this->local_assemblings_routines.push_back(NSResistanceMassMatrixSingle<d>);
    }
    else
    {
      this->local_assemblings_routines.push_back(NSResistanceMassMatrix<d>);
    }
  }

  // stabilization
  using namespace std::placeholders;
  double pspg_delta0 = db["pspg_delta0"];
  double gls_stab = db["gls_stab"];
  double characteristic_length = db["L_0"];

  if(pspg || symm_gls || nonsymm_gls) // need second derivatives
  {
    this->N_Terms = 2*d+2+d*(d+1)/2;
    auto soi = indices_up_to_order<d>(2);
    this->Derivatives.insert(this->Derivatives.end(), soi.begin()+1, soi.end());
    // Derivatives = { D000, D000, D100, D010, D001, D100, D010, D001,
    //                 D200, D110, D101, D020, D011, D002 }
    // or            { D00, D00, D10, D01, D10, D01, D20, D11, D02}
    for(int i = 0; i < d; ++i)
      this->FESpaceNumber.push_back(1);
    for(int i = 0; i < d*(d+1)/2; ++i)
      this->FESpaceNumber.push_back(0);
    // FESpaceNumber = {0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0}
    // or              {0, 1, 0, 0, 1, 1, 0, 0, 0}
    this->Needs2ndDerivatives[0] = true;
    if(pspg)
      this->local_assemblings_routines.push_back(
        std::bind(NSPSPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, pspg_delta0));
    else if(symm_gls)
      this->local_assemblings_routines.push_back(
        std::bind(NSsymmGLS<d>, _1, _2, _3, _4, _5, _6, _7, _8, gls_stab));
    else if(nonsymm_gls)
      this->local_assemblings_routines.push_back(
        std::bind(NSnonsymmGLS<d>, _1, _2, _3, _4, _5, _6, _7, _8, gls_stab, characteristic_length));
  }
  else if(brezzi_pitkaeranta)
  {
    this->N_Terms += d; // the pressure derivatives
    this->Derivatives.insert(this->Derivatives.end(), foi.begin()+1, foi.end());
    for(int i = 0; i < d; ++i)
      this->FESpaceNumber.push_back(1);
    this->local_assemblings_routines.push_back(
      std::bind(NS_BrezziPitkaeranta<d>, _1, _2, _3, _4, _5, _6, _7, _8,
                pspg_delta0));
  }

  // grad-div stabilization
  double graddiv_stab = db["graddiv_stab"];
  if(std::abs(graddiv_stab) > 1e-10)
  {
    this->local_assemblings_routines.push_back(
      std::bind(NSGradDiv<d>, _1, _2, _3, _4, _5, _6, _7, _8, graddiv_stab, characteristic_length));

    this->local_assemblings_routines.push_back(
        std::bind(NSGradDiv_RightHandSide<d>,
            _1, _2, _3, _4, _5, _6, _7, _8,
            graddiv_stab, characteristic_length));
  }

  switch(type)
  {
    case LocalAssembling_type::NavierStokesAll:
    {
    //this->local_assemblings_routines.push_back(NSDivergenceBlocks<d>);
    if(nonsymm_gls)
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>,
              _1, _2, _3, _4, _5, _6,
              _7, _8, -1));
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6,
              _7, _8, -1));
    }
    else
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6,
              _7, _8, 1));
    }


    if(nstype == 2 || nstype == 4 || nstype == 14)
    {
      this->local_assemblings_routines.push_back(NSGradientBlocks<d>);
    }

    if(pspg)
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSPSPG_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
              pspg_delta0));
    }
    else if(symm_gls)
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSsymmGLS_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
              gls_stab));
    }
    else if(nonsymm_gls)
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSnonsymmGLS_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7,
              _8, gls_stab, characteristic_length));
    }
    break;
  }
  case LocalAssembling_type::NavierStokesNL:
    {
      if(nonlin_form.is("convective"))
      {
        if(nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective<d>);
        }
      }
      else if(nonlin_form.is("skew_symmetric"))
      {
        if(nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
        }
      }
      else if(nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if(nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else if(nonlin_form.is("divergence"))
      {
        this->N_Parameters = d + d*d;
        this->N_ParamFct = 1;
        this->ParameterFct =  { NSParamsVelocityDerivatives<d> };
        this->N_FEValues = d + d*d;
        this->FEValue_FctIndex = std::vector<int>(d + d*d, 0);
        std::fill(this->FEValue_FctIndex.begin() + d+1,
                  this->FEValue_FctIndex.begin() + 2*d+2, 1);
        if(d == 3)
          std::fill(this->FEValue_FctIndex.begin()+2*d+2,
                    this->FEValue_FctIndex.end(), 2);
        auto ifo = indices_up_to_order<d>(1); // indices up to first order
        this->FEValue_MultiIndex = ifo;
        this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                        ifo.begin(), ifo.end());
        if(d == 3)
          this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                          ifo.begin(), ifo.end());
        this->BeginParameter = { 0 };

        this->local_assemblings_routines.push_back(NSNonlinearTerm_divergence<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }
      break;
    }
    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}
//========================================================================
template <int d>
void LocalAssembling<d>::set_parameters_for_nse_dg()
{
  // ( A B' )   ( 0 2 )
  // ( B C  )   ( 3 1 )
  this->N_Terms = 2*(d+1);
  //this->Derivatives = {D000, D000, D100, D010, D001, D100, D010, D001};
  auto fot = indices_up_to_order<d>(1);
  this->Derivatives = {fot[0]};
  this->Derivatives.insert(this->Derivatives.end(), fot.begin(), fot.end());
  this->Derivatives.insert(this->Derivatives.end(), fot.begin()+1,
                           fot.end());
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->N_Matrices = 4;
  this->RowSpace = { 0, 1, 0, 1 };
  this->ColumnSpace = { 0, 1, 1, 0 };
  this->N_Rhs = 2;
  this->RhsSpace = { 0, 1 };
  this->local_assemblings_routines.push_back(NSLaplaceGradGrad_dg<d>);


  if(type == LocalAssembling_type::NavierStokesNL)
  {
    this->local_assemblings_routines.push_back(NSNonlinearTerm_convective_dg<d>);
    this->N_Parameters = d;
        this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = 1;
  this->FEValue_FctIndex = std::vector<int>(1,0);
  this->FEValue_MultiIndex = MultiIndex_vector(1, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  }

  if(d == 3)
    this->FESpaceNumber = { 0, 1, 0, 0, 0, 1, 1, 1 }; // 0: velocity, 1: pressure
  else
    this->FESpaceNumber = { 0, 1, 0, 0, 1, 1};
  this->Manipulate = nullptr;
}

//========================================================================
template <int d>
void LocalAssembling<d>::set_parameters_for_nse_supg(LocalAssembling_type type)
{
  std::string disc_type = this->db["space_discretization_type"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  bool supg = (disc_type == std::string("supg"));

  int nstype = TDatabase::ParamDB->NSTYPE;

  if(nstype < 4)
  {
    ErrThrow("only NSTYPE 4 and 14 are supported");
  }
  if(!nonlin_form.is("convective"))
  {
    ErrThrow("Only convective form is supported yet");
  }

  if(!supg)
  {
    ErrThrow("unsupported space_discretization_type for NSE", d, "D: ",
             disc_type);
  }
  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = d+1;
  this->RhsSpace = { 0, 0, 0, 1 };

  this->Needs2ndDerivatives = new bool[d];
  // 2nd derivatives are not used in the local assembling routines at the moment
  this->N_Terms = 8;
  if(d==2)
  {
    for(int i=0; i<d;i++)
      this->Needs2ndDerivatives[i] = true;

    this->FESpaceNumber = { 0, 1, 0, 0, 1, 1, 0, 0 };
#ifdef __2D__         // 0   1    2    3    4    5    6    7
  this->Derivatives = { MultiIndex2D::D00, MultiIndex2D::D00, MultiIndex2D::D10,
                        MultiIndex2D::D01, MultiIndex2D::D10, MultiIndex2D::D01,
                        MultiIndex2D::D20, MultiIndex2D::D02};
#endif
  }
  else
  {
    for(int i=0; i<d;i++)
      this->Needs2ndDerivatives[i] = false;

    //                     u, p, ux,uy,uz, px,py,pz
    this->FESpaceNumber = { 0, 1, 0, 0, 0, 1, 1, 1};
#ifdef __3D__
  this->Derivatives = {MultiIndex3D::D000, // u
                       MultiIndex3D::D000, // p
                       MultiIndex3D::D100, // u_x
                       MultiIndex3D::D010, // u_y
                       MultiIndex3D::D001, // u_z
                       MultiIndex3D::D100, // p_x
                       MultiIndex3D::D010, // p_y
                       MultiIndex3D::D001};// p_z
#endif
  }

  this->Manipulate = nullptr;
  this->N_Parameters = d;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = d;
  this->FEValue_FctIndex = std::vector<int>(d);
  std::iota(this->FEValue_FctIndex.begin(), this->FEValue_FctIndex.end(), 0);
  this->FEValue_MultiIndex = MultiIndex_vector(d, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  this->N_Matrices = (d+1)*(d+1); // 9, 16
  if(d==3)
  {
    this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
  }
  else
  {
    this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
  }
  using namespace std::placeholders;
  double delta0 = db["supg_delta0"];
  double delta1 = db["graddiv_stab"];

  int problem_type = db["problem_type"];
  int fl_prb_type = TDatabase::ParamDB->FLOW_PROBLEM_TYPE;
  if(problem_type==3 || fl_prb_type==3)
  {
    this->local_assemblings_routines.push_back(
      std::bind(OseenSingleSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
    this->local_assemblings_routines.push_back(
      std::bind(OseenRhsSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
    this->local_assemblings_routines.push_back(
      std::bind(OseenGradBlockSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
  }
  else
  {
    // for both linear and nonlinear cases
    this->local_assemblings_routines.push_back(
      std::bind(NSLaplaceDeformationSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
    this->local_assemblings_routines.push_back(
      std::bind(NSRightHandSideSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
    this->local_assemblings_routines.push_back(
      std::bind(NSGradientBlocksSUPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
  }

  if(nstype == 14)
  {
    ErrThrow("nstype 14 is not yet implemented");
  }
  switch(type)
  {
    case LocalAssembling_type::NavierStokesAll:
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));
      break;
    }
    case LocalAssembling_type::NavierStokesNL:
    {
      break;
    }
    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}

//========================================================================

template<int d>
void LocalAssembling<d>::set_parameters_for_tnse()
{
  bool newtonian = this->db["viscosity_mode"].is("newtonian");
  bool laplace_type_deformation = this->db["laplace_type_deformation"];

  std::string disc_type = this->db["space_discretization_type"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);

  bool galerkin = (disc_type == std::string("galerkin"));
  bool smagorinsky = (disc_type==std::string("smagorinsky"));
  bool local_projection = (disc_type==std::string("local_projection"));
  bool jump_stab = (disc_type==std::string("jump_stab"));
  int nstype = TDatabase::ParamDB->NSTYPE;

  if (laplace_type_deformation)
  {
    if (nstype == 1 || nstype == 2)
    {
      ErrThrow("laplace_type_deformation is only supported for NSTYPE 3, 4, "
               "and 14");
    }
  }

  if (laplace_type_deformation != (TDatabase::ParamDB->LAPLACETYPE == 1))
  {
    Output::warn("inconsistent parameters in old and new database",
                 "The parameter 'laplace_type_deformation' is set to ",
                 laplace_type_deformation, " but LAPLACETYPE is ",
                 TDatabase::ParamDB->LAPLACETYPE, ". The latter is reset now!");

    TDatabase::ParamDB->LAPLACETYPE = (int)laplace_type_deformation;
  }

  if (nonlin_form.is("rotational") || nonlin_form.is("emac"))
  {
    if (nstype == 1 || nstype == 2)
    {
      ErrThrow("nse_nonlinear_form ", nonlin_form,
               " is only supported for NSTYPE 3, 4, and 14");
    }
  }

  if (!galerkin && !smagorinsky && !local_projection && !jump_stab)
  {
    ErrThrow("unsupported space_discretization_type for TNSE", d, "D: ",
             disc_type);
  }

  /*if (nstype == 14)
  {
    ErrThrow("NSTYPE", nstype, " is not supported yet");
  }*/

  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = d + 1;
  this->RhsSpace = { 0, 0, 0, 1 };
  this->N_Terms = d + 2;

  auto foi = indices_up_to_order<d>(1); // first_order_index
  this->Derivatives = indices_up_to_order<d>(0);
  this->Derivatives.insert(this->Derivatives.end(), foi.begin(), foi.end());
  // Derivatives = { D000, D000, D100, D010, D001 } or { D00, D00, D10, D01}

  this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure

  if (d == 3)
  {
    this->FESpaceNumber.push_back(0);
  }

  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;
  this->N_Parameters = d;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = d;
  this->FEValue_FctIndex = std::vector<int>(d);
  std::iota(this->FEValue_FctIndex.begin(), this->FEValue_FctIndex.end(), 0);
  this->FEValue_MultiIndex = MultiIndex_vector(d, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  this->N_Matrices = (d + 1) * (d + 1); // 9, 16

  if(db["space_discretization_type"].is("jump_stab"))
  {
    this->Needs2ndDerivatives[0] = true;
    this->Needs2ndDerivatives[1] = true;
    this->N_Parameters = 15;
    this->N_FEValues = 12;
    this->FEValue_FctIndex = 
       {
         0, 1, 
         0, 1, 
         0, 1, 
         0, 
         0, 0, 
         1, 
         1, 1
      };
    this->FEValue_MultiIndex = 
       { MultiIndex2D::D00, MultiIndex2D::D00, // u1, u2
         MultiIndex2D::D10, MultiIndex2D::D10, // u1x, u2x
         MultiIndex2D::D01, MultiIndex2D::D01, // u1y, u2y
         MultiIndex2D::D11,  // u1xy
         MultiIndex2D::D20, MultiIndex2D::D02, // u1xx, u1yy
         MultiIndex2D::D11, // u2xy
         MultiIndex2D::D20, MultiIndex2D::D02, // u2xx, u2yy
       };
    this->ParameterFct = {NSParamVelocityGradients};
  }
  if (d == 3)
  {
    this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
  }
  else
  {
    this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
  }

  using namespace std::placeholders;
  if ((type != LocalAssembling_type::TimeNavierStokesRhs
    && type != LocalAssembling_type::TimeNavierStokesMass)
    || type == LocalAssembling_type::NavierStokesLinear)
  {
    if (newtonian)
    {
      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
      }
      else
      {
        if (nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(NSLaplaceGradGradSingle<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(NSLaplaceGradGrad<d>);
        }
      }
    }
    else
    {
      ViscositySettings visc_settings(db);

      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceDeformationNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
      }
      else
      {
        if (nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceGradGradSingleNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
        }
        else
        {
          this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceGradGradNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
        }
      }
    }

    double characteristic_length = db["L_0"];

    // grad-div stabilization
    double graddiv_stab = db["graddiv_stab"];

    if (std::abs(graddiv_stab) > 1e-10)
    {
      if (nstype == 1 || nstype == 2)
      {
        ErrThrow("grad-div stabilization is only supported for NSTYPE 3, 4, "
                 "and 14");
      }

      this->local_assemblings_routines.push_back(
        std::bind(NSGradDiv<d>, _1, _2, _3, _4, _5, _6, _7, _8, graddiv_stab, characteristic_length));

      this->local_assemblings_routines.push_back(
        std::bind(NSGradDiv_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
                  graddiv_stab, characteristic_length));
    }

    // stabilization or turbulent models
    if (smagorinsky)
    {
      if (d == 2)
      {
        ErrThrow("Smagorinsky is not tested for 2D case");
      }

      this->N_Parameters = 15;
      this->N_ParamFct = 1;
      this->ParameterFct = { NSParamVelGradSmagorinsky<d> };
      this->N_FEValues = 12;
      this->FEValue_FctIndex = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 };
      // u1old, u2old, u3old, all derivatives of u1, u2, u3
#ifdef __3D__
      this->FEValue_MultiIndex = { MultiIndex3D::D000, MultiIndex3D::D000,
                                   MultiIndex3D::D000, MultiIndex3D::D100,
                                   MultiIndex3D::D100, MultiIndex3D::D100,
                                   MultiIndex3D::D010, MultiIndex3D::D010,
                                   MultiIndex3D::D010, MultiIndex3D::D001,
                                   MultiIndex3D::D001, MultiIndex3D::D001 };
#endif
      this->BeginParameter = { 0 };

      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(NSLaplaceDeformationSmagorinsky<d>);
      }
      else
      {
        if (nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(NSLaplaceGradGradSingleSmagorinsky<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(NSLaplaceGradGradSmagorinsky<d>);
        }
      }
    }
  }//endif LocalAssembling_type::TimeNavierStokesRhs

  switch (type)
  {
    case LocalAssembling_type::TimeNavierStokesAll:
    {
      this->local_assemblings_routines.push_back(
        std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));

      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>,
        _1, _2, _3, _4, _5, _6, _7, _8, -1));

      if (nstype == 2 || nstype == 4 || nstype == 14)
      {
        this->local_assemblings_routines.push_back(NSGradientBlocks<d>);
      }

      // the break statement is commented intentionally: The other option is to copy the
      // the complete lines of code for the "TimeNavierStokesNL". This is because
      // we need the linear and nonlinear matrices to assemble the system right hand side:
      // e.g., in the crank_nicolson time stepping scheme
      // THIS MEANS: we assemble all matrices at initial time
      // break;
      [[fallthrough]];
    }// NavierStokesAll

    case LocalAssembling_type::TimeNavierStokesNL:
    {
      if (nonlin_form.is("convective"))
      {
        if (nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective<d>);
        }
      }
      else if (nonlin_form.is("skew_symmetric"))
      {
        if (nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
        }
      }
      else if (nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if (nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else if (nonlin_form.is("divergence"))
      {
        this->N_Parameters = d + d * d;
        this->N_ParamFct = 1;
        this->ParameterFct =  { NSParamsVelocityDerivatives<d> };
        this->N_FEValues = d + d * d;
        this->FEValue_FctIndex = std::vector<int>(d + d * d, 0);

        std::fill(this->FEValue_FctIndex.begin() + d + 1,
                  this->FEValue_FctIndex.begin() + 2 * d + 2, 1);

        if (d == 3)
        {
          std::fill(this->FEValue_FctIndex.begin() + 2 * d + 2,
                    this->FEValue_FctIndex.end(), 2);
        }

        auto ifo = indices_up_to_order<d>(1); // indices up to first order
        this->FEValue_MultiIndex = ifo;
        this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                        ifo.begin(), ifo.end());
        if (d == 3)
        {
          this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                          ifo.begin(), ifo.end());
        }

        this->BeginParameter = { 0 };

        this->local_assemblings_routines.push_back(NSNonlinearTerm_divergence<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }
      break;
    }

    case LocalAssembling_type::TimeNavierStokesRhs:
      this->local_assemblings_routines.push_back(
        std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8, -1));
      break;

    case LocalAssembling_type::TimeNavierStokesMass:
      if (nstype == 1 || nstype == 2)
      {
        this->local_assemblings_routines.push_back(NSMassMatrixSingle<d>);
      }
      else
      {
        this->local_assemblings_routines.push_back(NSMassMatrix<d>);
      }
      break;

    case LocalAssembling_type::NavierStokesLinear:
      Output::print<5>("Nothing to do here!!!");
      break;

    case LocalAssembling_type::TimeNavierStokesExplNL:
    {
      // add the nonlinear term explicitly (using only already computed values
      // for the velocity) to the right-hand side
      // we need values and all derivatives of all velocity components
      this->N_Parameters = d + d * d;
      this->N_ParamFct = 1;
      this->ParameterFct =  { NSParamsVelocityDerivatives<d> };
      this->N_FEValues = d + d * d;
      this->FEValue_FctIndex = std::vector<int>(d + d * d, 0);
      std::fill(this->FEValue_FctIndex.begin() + d + 1,
                this->FEValue_FctIndex.begin() + 2 * d + 2, 1);

      if (d == 3)
      {
        std::fill(this->FEValue_FctIndex.begin() + 2 * d + 2,
                  this->FEValue_FctIndex.end(), 2);
      }

      auto ifo = indices_up_to_order<d>(1); // indices up to first order
      this->FEValue_MultiIndex = ifo;
      this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                      ifo.begin(), ifo.end());
      if (d == 3)
      {
        this->FEValue_MultiIndex.insert(this->FEValue_MultiIndex.end(),
                                        ifo.begin(), ifo.end());
      }

      this->BeginParameter = { 0 };
      this->local_assemblings_routines.push_back(NSRightHandSideExplicitNL<d>);
      break;
    }

    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}

//========================================================================

template<int d>
void LocalAssembling<d>::set_parameters_for_tnse_residual_vms()
{
  using namespace std::placeholders;

  if (!this->db["space_discretization_type"].is("residual_based_vms"))
  {
    ErrThrow("Unsupported space_discretization_type: ",
      this->db["space_discretization_type"]);
  }

  if (type == LocalAssembling_type::TimeNavierStokesExplNL)
  {
    ErrThrow("Fully explicit mode is not yet supported with residual_based_vms.");
  }

  bool newtonian = this->db["viscosity_mode"].is("newtonian");
  bool laplace_type_deformation = db["laplace_type_deformation"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);

  int nstype = TDatabase::ParamDB->NSTYPE;

  if (nstype != 14)
  {
    ErrThrow("residual_based_vms is only supported with NSTYPE = 14.");
  }

  if (nonlin_form.is("divergence"))
  {
    ErrThrow("nse_nonlinear_form ", nonlin_form,
     " is not supported with residual_based_vms.");
  }

  RBVMS_Settings param_settings(db);

  // we need cell-local information to compute the stabilization parameters
  this->assembly_needs_reftrans = true;

  this->N_Rhs = d + 1;
  if (d == 2)
  {
    this->RhsSpace = { 0, 0, 1 };
  }
  else
  {
    this->RhsSpace = { 0, 0, 0, 1 };
  }

  // u, p, Du, Dp, D²u
  // (each velocity component is from the same space, so only one copy)
  this->Derivatives = indices_of_order<d>(0); // 0: u

  // 1: p
  for (auto mi: indices_of_order<d>(0))
  {
    this->Derivatives.push_back(mi);
  }

  // 2: Du
  for (auto mi: indices_of_order<d>(1))
  {
    this->Derivatives.push_back(mi);
  }

  // d + 2: Dp
  for (auto mi: indices_of_order<d>(1))
  {
    this->Derivatives.push_back(mi);
  }

  // 2d + 2: D²u
  for (auto mi: indices_of_order<d>(2))
  {
    this->Derivatives.push_back(mi);
  }

  // 0: velocity, 1: pressure
  if (d == 2)
  {
    this->FESpaceNumber = { 0, 1, // u, p
      0, 0, 1, 1,                 // Du, Dp
      0, 0, 0, 0 };               // D²u
  }
  else
  {
    this->FESpaceNumber = { 0, 1,  // u, p
      0, 0, 0, 1, 1, 1,            // Du, Dp
      0, 0, 0, 0, 0, 0, 0, 0, 0 }; // D²u
  }
  this->N_Terms = this->Derivatives.size();

  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = true;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;

  // set the FEValues fields and parameter functions
  set_parameter_inputs_and_functions_for_tnse_residual_vms(false);

  if (type != LocalAssembling_type::TimeNavierStokesMass)
  {
    this->N_Matrices = (d + 1) * (d + 1); // 9, 16

    // spaces for matrix blocks:
    //
    // A  B^T
    // B  C
    //
    // where
    // - A has d \times d velocity-velocity subblocks Aij,
    // - B has d pressure-velocity column subblocks Bi,
    // - B^T has d velocity-pressure row subblocks Bj^T,
    // - C is a single pressure-pressure block.

    // dxd A blocks start at 0
    // 1 C block is at d * d
    // d B blocks start at d * d + 1
    // d B^T blocks start at d * d + 1 + d = (d + 1) * d + 1

    if (d == 3)
    {
      // V, V: A11, A12, A13, A21, A22, A23, A31, A32, A33,
      // Q, Q: C,
      // Q, V: B1, B2, B3,
      // V, Q: B1^T, B2^T, B3^T
      this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
    }
    else
    {
      // V, V: A11, A12, A21, A22,
      // Q, Q: C,
      // Q, V: B1, B2,
      // V, Q: B1^T, B2^T
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
    }
  }
  else
  {
    this->N_Matrices = d * d + d;

    // spaces for matrix blocks:
    //
    // M  0
    // MQ 0
    //
    // where
    // - M has d \times d velocity-velocity subblocks Mij,
    // - MQ has d pressure-velocity column subblocks Mqi.

    // dxd M blocks start at 0
    // d MQ blocks start at d * d

    if (d == 3)
    {
      // V, V: M11, M12, M13, M21, M22, M23, M31, M32, M33,
      // Q, V: MQ1, MQ2, MQ3
      this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,    0, 0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0 };
      // HACK: +4 irrelevant entries because some code elsewhere makes bad assumptions
    }
    else
    {
      // V, V: M11, M12, M21, M22,
      // Q, V: MQ1, MQ2
      this->RowSpace =    { 0, 0, 0, 0, 1, 1,    0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0,    0, 0, 0 };
      // HACK: +3 irrelevant entries because some code elsewhere makes bad assumptions
    }
  }

  if ((type != LocalAssembling_type::TimeNavierStokesRhs
    && type != LocalAssembling_type::TimeNavierStokesMass))
  {
    if (newtonian)
    {
      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
      }
      else
      {
        this->local_assemblings_routines.push_back(NSLaplaceGradGrad<d>);
      }
    }
    else
    {
      ViscositySettings visc_settings(db);

      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceDeformationNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
      }
      else
      {
        this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceGradGradNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
      }
    }

    // grad-div stabilization
    double graddiv_stab = db["graddiv_stab"];

    if (std::abs(graddiv_stab) > 1e-10)
    {
      ErrThrow("Grad-div stabilization is not supported in combination with "
        "residual-based VMS.");
    }
  }//endif LocalAssembling_type::TimeNavierStokesRhs

  switch (type)
  {
    case LocalAssembling_type::TimeNavierStokesAll:
    {
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>,
        _1, _2, _3, _4, _5, _6, _7, _8, -1));

      this->local_assemblings_routines.push_back(std::bind(NSVMSResiduals_RHS<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));

      // the break statement is commented intentionally: The other option is to copy the
      // the complete lines of code for the "TimeNavierStokesNL". This is because
      // we need the linear and nonlinear matrices to assemble the system right hand side:
      // e.g., in the crank_nicolson time stepping scheme
      // THIS MEANS: we assemble all matrices at initial time
      // break;
      [[fallthrough]];
    }// NavierStokesAll

    case LocalAssembling_type::TimeNavierStokesNL:
    {
      if (nonlin_form.is("convective"))
      {
        this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective<d>);
      }
      else if (nonlin_form.is("skew_symmetric"))
      {
        this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
      }
      else if (nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if (nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }

      this->local_assemblings_routines.push_back(
        std::bind(NSVMSResiduals<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          coeff_ref_inv_jacobian_offset, param_settings));

      break;
    }

    case LocalAssembling_type::TimeNavierStokesRhs:
      this->local_assemblings_routines.push_back(
        std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8, -1));

      this->local_assemblings_routines.push_back(std::bind(NSVMSResiduals_RHS<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));
      break;

    case LocalAssembling_type::TimeNavierStokesMass:
      this->local_assemblings_routines.push_back(std::bind(NSVMSResiduals_MassMatrix<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));
      break;

    case LocalAssembling_type::NavierStokesLinear:
      Output::print<5>("Nothing to do here!!!");
      break;

    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}

//========================================================================

template<int d>
void LocalAssembling<d>::set_parameters_for_tnse_rbvms_time()
{
  using namespace std::placeholders;

  if (!this->db["space_discretization_type"].is("rbvms_time"))
  {
    ErrThrow("Unsupported space_discretization_type: ",
      this->db["space_discretization_type"]);
  }

  if (type == LocalAssembling_type::TimeNavierStokesExplNL)
  {
    ErrThrow("Fully explicit mode is not yet supported with rbvms_time.");
  }

  bool newtonian = this->db["viscosity_mode"].is("newtonian");
  bool laplace_type_deformation = db["laplace_type_deformation"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);

  int nstype = TDatabase::ParamDB->NSTYPE;

  if (nstype != 14)
  {
    ErrThrow("rbvms_time is only supported with NSTYPE = 14.");
  }

  if (nonlin_form.is("divergence"))
  {
    ErrThrow("nse_nonlinear_form ", nonlin_form,
     " is not supported with rbvms_time.");
  }

  RBVMS_Settings param_settings(db);

  // we need cell-local information to compute the stabilization parameters
  this->assembly_needs_reftrans = true;

  // residual and subscale velocity per quadrature point
  this->N_PersistentDataPerPoint = 2 * d;

  this->N_Rhs = d + 1;
  if (d == 2)
  {
    this->RhsSpace = { 0, 0, 1 };
  }
  else
  {
    this->RhsSpace = { 0, 0, 0, 1 };
  }

  // u, p, Du, Dp, D²u
  // (each velocity component is from the same space, so only one copy)
  this->Derivatives = indices_of_order<d>(0); // 0: u

  // 1: p
  for (auto mi: indices_of_order<d>(0))
  {
    this->Derivatives.push_back(mi);
  }

  // 2: Du
  for (auto mi: indices_of_order<d>(1))
  {
    this->Derivatives.push_back(mi);
  }

  // d + 2: Dp
  for (auto mi: indices_of_order<d>(1))
  {
    this->Derivatives.push_back(mi);
  }

  // 2d + 2: D²u
  for (auto mi: indices_of_order<d>(2))
  {
    this->Derivatives.push_back(mi);
  }

  // 0: velocity, 1: pressure
  if (d == 2)
  {
    this->FESpaceNumber = { 0, 1, // u, p
      0, 0, 1, 1,                 // Du, Dp
      0, 0, 0, 0 };               // D²u
  }
  else
  {
    this->FESpaceNumber = { 0, 1,  // u, p
      0, 0, 0, 1, 1, 1,            // Du, Dp
      0, 0, 0, 0, 0, 0, 0, 0, 0 }; // D²u
  }
  this->N_Terms = this->Derivatives.size();

  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = true;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;

  // set the FEValues fields and parameter functions (same as for RB-VMS)
  set_parameter_inputs_and_functions_for_tnse_residual_vms(true);

  if (type != LocalAssembling_type::TimeNavierStokesMass)
  {
    this->N_Matrices = (d + 1) * (d + 1); // 9, 16

    // spaces for matrix blocks:
    //
    // A  B^T
    // B  C
    //
    // where
    // - A has d \times d velocity-velocity subblocks Aij,
    // - B has d pressure-velocity column subblocks Bi,
    // - B^T has d velocity-pressure row subblocks Bj^T,
    // - C is a single pressure-pressure block.

    // dxd A blocks start at 0
    // 1 C block is at d * d
    // d B blocks start at d * d + 1
    // d B^T blocks start at d * d + 1 + d = (d + 1) * d + 1

    if (d == 3)
    {
      // V, V: A11, A12, A13, A21, A22, A23, A31, A32, A33,
      // Q, Q: C,
      // Q, V: B1, B2, B3,
      // V, Q: B1^T, B2^T, B3^T
      this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
    }
    else
    {
      // V, V: A11, A12, A21, A22,
      // Q, Q: C,
      // Q, V: B1, B2,
      // V, Q: B1^T, B2^T
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
    }
  }
  else
  {
    this->N_Matrices = d * d + d;

    // spaces for matrix blocks:
    //
    // M  0
    // MQ 0
    //
    // where
    // - M has d \times d velocity-velocity subblocks Mij,
    // - MQ has d pressure-velocity column subblocks Mqi.

    // dxd M blocks start at 0
    // d MQ blocks start at d * d

    if (d == 3)
    {
      // V, V: M11, M12, M13, M21, M22, M23, M31, M32, M33,
      // Q, V: MQ1, MQ2, MQ3
      this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,    0, 0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0 };
      // HACK: +4 irrelevant entries because some code elsewhere makes bad assumptions
    }
    else
    {
      // V, V: M11, M12, M21, M22,
      // Q, V: MQ1, MQ2
      this->RowSpace =    { 0, 0, 0, 0, 1, 1,    0, 0, 0 };
      this->ColumnSpace = { 0, 0, 0, 0, 0, 0,    0, 0, 0 };
      // HACK: +3 irrelevant entries because some code elsewhere makes bad assumptions
    }
  }

  if ((type != LocalAssembling_type::TimeNavierStokesRhs
    && type != LocalAssembling_type::TimeNavierStokesMass))
  {
    if (newtonian)
    {
      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
      }
      else
      {
        this->local_assemblings_routines.push_back(NSLaplaceGradGrad<d>);
      }
    }
    else
    {
      ViscositySettings visc_settings(db);

      if (laplace_type_deformation)
      {
        this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceDeformationNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
      }
      else
      {
        this->local_assemblings_routines.push_back(std::bind(
          NSLaplaceGradGradNonNewtonian<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          visc_settings));
      }
    }

    // grad-div stabilization
    double graddiv_stab = db["graddiv_stab"];

    if (std::abs(graddiv_stab) > 1e-10)
    {
      ErrThrow("Grad-div stabilization is not supported in combination with "
        "residual-based VMS.");
    }
  }//endif LocalAssembling_type::TimeNavierStokesRhs

  switch (type)
  {
    case LocalAssembling_type::TimeNavierStokesAll:
    {
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>,
        _1, _2, _3, _4, _5, _6, _7, _8, -1));

      this->local_assemblings_routines.push_back(std::bind(TNSE_RBVMS_Time_SubgridTerms_RHS<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));

      // the break statement is commented intentionally: The other option is to copy the
      // the complete lines of code for the "TimeNavierStokesNL". This is because
      // we need the linear and nonlinear matrices to assemble the system right hand side:
      // e.g., in the crank_nicolson time stepping scheme
      // THIS MEANS: we assemble all matrices at initial time
      // break;
      [[fallthrough]];
    }// NavierStokesAll

    case LocalAssembling_type::TimeNavierStokesNL:
    {
      if (nonlin_form.is("convective"))
      {
        this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective<d>);
      }
      else if (nonlin_form.is("skew_symmetric"))
      {
        this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
      }
      else if (nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if (nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }

      this->local_assemblings_routines.push_back(
        std::bind(TNSE_RBVMS_Time_SubgridTerms_Matrices<d>, _1, _2, _3, _4, _5, _6, _7, _8,
          coeff_ref_inv_jacobian_offset, param_settings));

      break;
    }

    case LocalAssembling_type::TimeNavierStokesRhs:
      this->local_assemblings_routines.push_back(
        std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8, -1));

      this->local_assemblings_routines.push_back(std::bind(TNSE_RBVMS_Time_SubgridTerms_RHS<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));
      break;

    case LocalAssembling_type::TimeNavierStokesMass:
      this->local_assemblings_routines.push_back(std::bind(TNSE_RBVMS_Time_MassMatrix<d>,
        _1, _2, _3, _4, _5, _6, _7, _8,
        coeff_ref_inv_jacobian_offset, param_settings));
      break;

    case LocalAssembling_type::NavierStokesLinear:
      Output::print<5>("Nothing to do here!!!");
      break;

    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}

template <int d>
void RBVMS_Residuals_PlainAdvection(const double *in, double *out)
{
  NSParamsVMSResidualsWithoutLaplacianAndRHS<d>(in, out, false);
}

template <int d>
void RBVMS_Residuals_ExtendedAdvection(const double *in, double *out)
{
  NSParamsVMSResidualsWithoutLaplacianAndRHS<d>(in, out, true);
}

template<int d>
void LocalAssembling<d>::set_parameter_inputs_and_functions_for_tnse_residual_vms(bool extend_advection)
{
  // NSParamsVelocity:
  //   u (d components)
  // NSParamsVMSResidualsLaplacian:
  //   Lu (d components)
  // NSParamsVMSResidualsWithoutLaplacianAndRHS:
  //   res_m without laplacian and rhs (d components)
  //   res_c (1 component)
  // NSParamsOldVelocity:
  //   u_old (d components)
  // NSParamsVelocityGradient:
  //   Du (d^2 components)

  using namespace std::placeholders;

  this->N_Parameters = d + d + (d + 1) + d + d * d;
  this->N_ParamFct = 4;
  this->ParameterFct = { NSParamsVelocity<d>,
    NSParamsVMSResidualsLaplacian<d>,
    extend_advection ? RBVMS_Residuals_ExtendedAdvection<d> : RBVMS_Residuals_PlainAdvection<d>,
    NSParamsOldVelocity<d>,
    NSParamsVelocityGradient<d> };
  this->BeginParameter =
  {
    0,
    d,
    2 * d,
    3 * d + 1,
    4 * d + 1 };

  // FE values:
  //
  //  0,  1(,  2),
  // u1, u2(, u3),
  //
  //   d,  2d(,  3d),
  // Du1, Du2(, Du3),
  //
  // (d + 1) d, (d + 2) d(, (d + 3) d)
  //       Lu1,       Lu2(,       Lu3)
  //
  // (2d + 1) d,
  // Dp,
  //
  // (2d + 2)d, (2d + 2)d + 1(, (2d + 2)d + 2)
  //   u1_old,    u2_old(,         u3_old    )

  //                 u      Du        Lu    Dp  u_old
  this->N_FEValues = d + (d * d) + (d * d) + d + d;

  this->FEValue_FctIndex.clear();

  // u
  for (int i = 0; i < d; i++)
  {
    this->FEValue_FctIndex.push_back(i);
  }

  // Du
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      this->FEValue_FctIndex.push_back(i);
    }
  }

  // Lu
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      this->FEValue_FctIndex.push_back(i);
    }
  }

  // Dp
  for (int i = 0; i < d; i++)
  {
    this->FEValue_FctIndex.push_back(d);
  }

  // u_old
  for (int i = 0; i < d; i++)
  {
    this->FEValue_FctIndex.push_back(d + 1 + i);
  }

  this->FEValue_MultiIndex.clear();

  // u
  for (int i = 0; i < d; i++)
  {
    for (auto mi: indices_of_order<d>(0))
    {
      this->FEValue_MultiIndex.push_back(mi);
    }
  }

  // Du
  for (int i = 0; i < d; i++)
  {
    for (auto mi: indices_of_order<d>(1))
    {
      this->FEValue_MultiIndex.push_back(mi);
    }
  }

  // Lu
  for (int i = 0; i < d; i++)
  {
    for (auto mi: laplacian_indices<d>())
    {
      this->FEValue_MultiIndex.push_back(mi);
    }
  }

  // Dp
  for (auto mi: indices_of_order<d>(1))
  {
    this->FEValue_MultiIndex.push_back(mi);
  }

  // u_old
  for (int i = 0; i < d; i++)
  {
    for (auto mi: indices_of_order<d>(0))
    {
      this->FEValue_MultiIndex.push_back(mi);
    }
  }
}

//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_tnse_vms( LocalAssembling_type la_type)
{
  bool laplace_type_deformation = this->db["laplace_type_deformation"];
  std::string disc_type = this->db["space_discretization_type"];
  bool vms_projection = (disc_type == std::string("vms_projection"));
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  int nstype = TDatabase::ParamDB->NSTYPE;

  if(!laplace_type_deformation && nstype != 4)
  {
    ErrThrow("vms_projection is only supported for NSTYPE 4 !!!");
  }

  if(nonlin_form.is("rotational") || nonlin_form.is("emac"))
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("nse_nonlinear_form ", nonlin_form,
               " is only supported for NSTYPE 3, 4, and 14");
    }
  }
  if(!vms_projection)
  {
    ErrThrow("unsupported space_discretization_type for NSE", 3, "D: ",
             disc_type);
  }
  if(nstype == 14)
  {
    ErrThrow("NSTYPE", nstype, " is not supported yet");
  }
  // common for all NSTYPE, Discrete forms, etc

  this->N_Parameters = 16;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVariationalMSLargeScale<d> };
  this->N_FEValues = 13;
  this->FEValue_FctIndex = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 3};
#ifdef __3D__
  this->FEValue_MultiIndex = { MultiIndex3D::D000, // u1old,
                               MultiIndex3D::D000, // u2old
                               MultiIndex3D::D000, // u3old
                               MultiIndex3D::D100, // u1old_x
                               MultiIndex3D::D100, // u2old_x
                               MultiIndex3D::D100, // u3old_x
                               MultiIndex3D::D010, // u1old_y
                               MultiIndex3D::D010, // u2old_y
                               MultiIndex3D::D010, // u3old_y
                               MultiIndex3D::D001, // u1old_z
                               MultiIndex3D::D001, // u2old_z
                               MultiIndex3D::D001, // u3old_z
                               MultiIndex3D::D000  // pold
  };
  this->Derivatives = {MultiIndex3D::D000, MultiIndex3D::D000,
                       MultiIndex3D::D100, MultiIndex3D::D010,
                       MultiIndex3D::D001, MultiIndex3D::D000};
#endif
  this->N_Terms = 6;
  this->Needs2ndDerivatives = new bool[3];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Needs2ndDerivatives[2] = false;
  this->FESpaceNumber = { 0, 1, 0, 0, 0, 2 }; // 0: velocity, 1: pressure, 2: projection
  this->N_Matrices = 22;
  //                   ------------A------------, L, .......Bs....... , .GTilde,.,.G...
  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 2, 2 };
  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0};
  this->N_Rhs = 4;
  this->RhsSpace = { 0, 0, 0, 1};
  this->Manipulate = NULL;
  this->BeginParameter = { 0 };

  using namespace std::placeholders;
  if((la_type != LocalAssembling_type::TimeNavierStokesRhs
    && la_type != LocalAssembling_type::TimeNavierStokesMass)
    || la_type == LocalAssembling_type::NavierStokesLinear)
  {
    this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
    // additional part from the turbulence model
    this->local_assemblings_routines.push_back(NSLaplaceDeformationVariationalMS<d>);
  }
  switch(la_type)
  {
    case LocalAssembling_type::TimeNavierStokesAll:
    {
      this->local_assemblings_routines.push_back(
        std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));

      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6,
        _7, _8, -1));

      this->local_assemblings_routines.push_back(NSGradientBlocks<d>);

      this->local_assemblings_routines.push_back(NSVariationlMS_GMatrices<d>);
      // the break statement is commented intentionally: The other option is to copy the
      // the complete lines of code for the "TimeNavierStokesNL". This is because
      // we need the linear and nonlinear matrices to assemble the system right hand side:
      // e.g., in the crank_nicolson time stepping scheme
      // THIS MEANS: we assemble all matrices at initial time
      // break;
      [[fallthrough]];
    }// NavierStokesAll
    case LocalAssembling_type::TimeNavierStokesNL:
    {
      if(nonlin_form.is("convective"))
      {
        this->local_assemblings_routines.push_back(
               NSNonlinearTerm_convective<d>);
      }
      else if(nonlin_form.is("skew_symmetric"))
      {
        this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
      }
      else if(nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if(nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }
      // nonlinear matrices Gtilde
      this->local_assemblings_routines.push_back(NSVariationlMS_GTildeMatrices<d>);
      break;
    }
    case LocalAssembling_type::TimeNavierStokesRhs:
      this->local_assemblings_routines.push_back(
        std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8, -1));
      break;
    case LocalAssembling_type::TimeNavierStokesMass:
        this->local_assemblings_routines.push_back(NSMassMatrix<d>);
        // lumped mass matrix
        this->local_assemblings_routines.push_back(NSLumpMassMatrix<d>);
      break;
    case LocalAssembling_type::NavierStokesLinear:
      Output::print<5>("Nothing to do here!!!");
      break;
    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}

template<int d>
void LocalAssembling<d>::set_parameters_for_tnse_supg( LocalAssembling_type la_type)
{
  int nstype = TDatabase::ParamDB->NSTYPE;
  if(nstype < 4)
  {
    ErrThrow("SUPG method is only supported for NSTYPE 4, "
             "and 14");
  }
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    this->N_FEValues = d;
    if(d==2)
    {
      this->N_Parameters = 2;
      this->FEValue_FctIndex = { 0, 1 };
    }
    else
    {
      this->N_Parameters = 6; // why 6? I have to remember :(
      this->FEValue_FctIndex = { 0, 1, 2};
    }

    this->N_ParamFct = 1;
    this->ParameterFct =  { NSParamsVelocityDerivatives_SUPG_inf_sup<d> };

#ifdef __2D__
    this->FEValue_MultiIndex = { MultiIndex2D::D00, MultiIndex2D::D00 };
#else
    this->FEValue_MultiIndex = { MultiIndex3D::D000, MultiIndex3D::D000,
                                 MultiIndex3D::D000};

#endif
    this->BeginParameter = { 0 };
  }

  if(TDatabase::ParamDB->NSTYPE==14)
  {
    this->N_FEValues = d*2;
    if(d==2)
    {
      this->N_Parameters = 4;
      this->FEValue_FctIndex = { 0, 1, 2, 3 };
    }
    else
    {
      this->N_Parameters = 9; // why 9?
      this->FEValue_FctIndex = { 0, 1, 2,
                                 3, 4, 5 };
    }
    this->N_ParamFct = 1;
    this->ParameterFct =  { NSParamsVelocityDerivatives_SUPG_equal_order<d> };
#ifdef __2D__
    this->FEValue_MultiIndex = { MultiIndex2D::D00, MultiIndex2D::D00,
                                 MultiIndex2D::D00, MultiIndex2D::D00 };
#else
    this->FEValue_MultiIndex = { MultiIndex3D::D000, MultiIndex3D::D000,
                                 MultiIndex3D::D000, MultiIndex3D::D000,
                                 MultiIndex3D::D000, MultiIndex3D::D000 };
#endif
    this->BeginParameter = { 0 };
  }

  this->Needs2ndDerivatives = new bool[d];
  // 2nd derivatives are not used in the local assembling routines at the moment
  if(d==2)
  {
    for(int i=0; i<d;i++)
      this->Needs2ndDerivatives[i] = true;
  }
  else
  {
    for(int i=0; i<d;i++)
      this->Needs2ndDerivatives[i] = false;
  }

  this->N_Matrices = (d+1)*(d+1);
  this->N_Terms = 8; // Derivatives
  if(d==2)
  {
    this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
    //                      u, p, ux,uy,px,py,uxx,uyy
    this->FESpaceNumber = { 0, 1, 0, 0, 1, 1, 0, 0 };
#ifdef __2D__         // 0   1    2    3    4    5    6    7
  this->Derivatives = { MultiIndex2D::D00, MultiIndex2D::D00, MultiIndex2D::D10,
                        MultiIndex2D::D01, MultiIndex2D::D10, MultiIndex2D::D01,
                        MultiIndex2D::D20, MultiIndex2D::D02};
#endif
  }
  else
  {
    this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
    //                     u, p, ux,uy,uz, px,py,pz
    this->FESpaceNumber = { 0, 1, 0, 0, 0, 1, 1, 1};
#ifdef __3D__
  this->Derivatives = {MultiIndex3D::D000, // u
                       MultiIndex3D::D000, // p
                       MultiIndex3D::D100, // u_x
                       MultiIndex3D::D010, // u_y
                       MultiIndex3D::D001, // u_z
                       MultiIndex3D::D100, // p_x
                       MultiIndex3D::D010, // p_y
                       MultiIndex3D::D001};// p_z
#endif
  }
  this->N_Rhs = d+1;
  this->RhsSpace = { 0, 0, 0, 1 };
  this->Manipulate = nullptr;
  double delta0 = db["supg_delta0"];
  double delta1 = db["graddiv_stab"];
  // coefficient from BDF2 scheme scaled with tau;
  double tau = db["time_step_length"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);

  using namespace std::placeholders;
  switch(la_type)
  {
    case LocalAssembling_type::TimeNavierStokesAll:
      // velocity-velocity blocks
      if(nonlin_form.is("convective"))
      {
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
      }
      else if(nonlin_form.is("skew_symmetric"))
      {
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_skew_symmetric<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
      }
      else if(nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_rotational<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
      }
      else if(nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_emac<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }
      if(nstype == 4)
      {
        // pressure-velocity blocks
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_GradientBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
        this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));
        // right hand side
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_RightHandSide_InfSup<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
      }
      else // type 14
      {
        double factor = 3./(2.*tau);
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_EquOrder_Gradient_DivergenceBlocks<d>,
                    _1, _2, _3, _4, _5, _6, _7, _8, factor));

        factor = db["time_step_length"];
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_RightHandSide_EquOrder<d>, _1, _2, _3, _4, _5, _6, _7, _8, factor));
      }
      break;
    case LocalAssembling_type::TimeNavierStokesMass:
      // mass matrices are same for type 4 and 14 except the stabilization parameters
      this->local_assemblings_routines.push_back(
        std::bind(NS_SUPG_MassMatrix<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0));
      break;
    case LocalAssembling_type::TimeNavierStokesNL:
      if(nstype == 4)
      {
        // velocity-velocity blocks
        if(nonlin_form.is("convective"))
        {
          this->local_assemblings_routines.push_back(
            std::bind(NS_SUPG<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
        }
        else if(nonlin_form.is("skew_symmetric"))
        {
          this->local_assemblings_routines.push_back(
            std::bind(NS_SUPG_skew_symmetric<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
        }
        else if(nonlin_form.is("rotational"))
        {
          this->local_assemblings_routines.push_back(
            std::bind(NS_SUPG_rotational<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
        }
        else if(nonlin_form.is("emac"))
        {
          this->local_assemblings_routines.push_back(
            std::bind(NS_SUPG_emac<d>,_1, _2, _3, _4, _5, _6, _7, _8, delta0, delta1));
        }
        else
        {
          ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
        }
        // pressure-velocity blocks
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_GradientBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
        // right hand side
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_RightHandSide_InfSup<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
      }
      else
      {
        ErrThrow("NSTYPE 14 uses TimeNavierStokesAll to assemble, please check the main class ");
      }
      break;
    case LocalAssembling_type::TimeNavierStokesRhs:
      if(nstype == 4)
      {
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_RightHandSide_InfSup<d>, _1, _2, _3, _4, _5, _6, _7, _8, delta0));
      }
      else
      {
        double factor = db["time_step_length"];
        this->local_assemblings_routines.push_back(
          std::bind(NS_SUPG_RightHandSide_EquOrder<d>, _1, _2, _3, _4, _5, _6, _7, _8, factor));
      }
      break;
      default:
        ErrThrow("unknown LocalAssembling_type ", la_type);
      break;
  }
}

#ifdef __3D__
template class LocalAssembling<3>;
#endif
#ifdef __2D__
template class LocalAssembling<2>;
#endif
