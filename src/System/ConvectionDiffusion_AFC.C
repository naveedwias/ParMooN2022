#include "ConvectionDiffusion_AFC.h"
#include "Database.h"
#include "ConvDiff.h"
#include "AlgebraicFluxCorrection.h"
#include "LinAlg.h"
#include "anderson.h"
#include "MainUtilities.h"

#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "Example_CD2D.h"
 #include "AuxParam2D.h"
 #include "Upwind.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h" 
 #include "Example_CD3D.h"
 #include "Upwind3D.h"
 #include "FEFunction2D.h"
#endif

#ifdef _MPI
 #include "ParFECommunicator3D.h"
 #include "MumpsWrapper.h"
#endif

#include <cmath>

/* ************************************************************************* */
/**
 * @brief: While using AFC, system matrices is created as if
 *         all DOFs were active.
 */
const ParameterDatabase& test(const ParameterDatabase& db)
{
  return db;
}

/* ************************************************************************ */

template<int d>
ConvectionDiffusion_AFC<d>::ConvectionDiffusion_AFC(const TDomain& domain, 
						 const ParameterDatabase& param_db)
: ConvectionDiffusion<d>(domain, test(param_db), Example_CD(param_db)),
  afc_matrix_D(*this->ConvectionDiffusion<d>::systems.front().matrix.get_blocks_uniquely()[0]),
  afc_matrix_D_B(afc_matrix_D), 
  afc_alphas(afc_matrix_D.get_n_entries(), 1.0), diffusion_matrix(afc_matrix_D), 
  alphas_x_i(), old_solution(),
  initial_solution(), poisson_sol()
{
  // A default AFC database. Use to merge AFC database with Parameter DB
  auto afc_db(AlgebraicFluxCorrection::default_afc_database());
  afc_db.merge(param_db);
  this->db.merge(afc_db, true);
  
  set_AFC_parameters();
  //Currently, AFC doesn't support Multigrid
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("AFC doesn't support Multigrid!");
  }
  alphas_x_i = this->ConvectionDiffusion<d>::systems.front().solution;
  old_solution = this->ConvectionDiffusion<d>::systems.front().solution;
  initial_solution = BlockVector(this->ConvectionDiffusion<d>::systems.front().matrix, 
                                 false);
  
  poisson_sol = BlockVector(this->ConvectionDiffusion<d>::systems.front().matrix, 
                            false);
  //Set thes values to zero
  poisson_sol.reset();
  diffusion_matrix.reset();
  
  is_not_afc_fixed_point_rhs=1;
  rejected_steps=0;
  
  /*
   * Creates an array to store the interpolation of the exact solution. Useful
   * in finding AFC norm and d_h( ; , ) error.
   */
  
  auto& space = this->ConvectionDiffusion<d>::systems.front().fe_space;
  exact_interpolant.resize(space->get_n_dof());
  initial_fe_function = FEFunction(space, "c", initial_solution.get_entries());
  
  //creates a FEFUnction  
  FEFunction exact_interpolation(space, "interpolant",
                                 exact_interpolant.data());
  //Interpolates and store the value in exact_interpolant
  exact_interpolation.Interpolate(this->example.get_exact(0));
#ifdef __3D__
  TAuxParam3D aux;
  MultiIndex3D AllDerivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                     MultiIndex3D::D010, MultiIndex3D::D001 };
#else
  TAuxParam2D aux;
  MultiIndex2D AllDerivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                     MultiIndex2D::D01 };
#endif
  const FESpace*  fe_space_pointer = space.get();
  exact_interpolation.GetErrors(this->example.get_exact(0), d+1, AllDerivatives,
                                ConvectionDiffusion<d>::n_errors,
                                conv_diff_l2_h1_linf_error<d>,
                                this->example.get_coeffs(), &aux, 1,
                                &fe_space_pointer, interpolation_errors.data());
}


/** ************************************************************************ */
template<int d>
void ConvectionDiffusion_AFC<d>::set_AFC_parameters()
{
  if(!this->db["algebraic_flux_correction"].is("none"))
  {
    if(!this->db["algebraic_flux_correction"].is("afc"))
    {
      this->db["algebraic_flux_correction"].set("afc");
      Output::print("Only kind of algebraic flux correction"
        " for CD problems is AFC (afc).");
    }
    //make sure that galerkin discretization is used
    if (!this->db["space_discretization_type"].is("galerkin"))
    {                                            
      this->db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to", 
       "Galerkin because Algebraic Flux Correction is enabled.");
    }
  }
}
/** ************************************************************************ */
template<int d>
void ConvectionDiffusion_AFC<d>::assemble(const int iteration)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  //determine the local assembling type to be ConvectionDiffusion
  LocalAssembling_type laTypet = LocalAssembling_type::ConvDiff;
  AlgebraicFluxCorrection::Limiter limiter =
                       string_to_limiter(this->db["afc_limiter"]);
  
  if(iteration == 0)
    time_total = GetTime();
  
  auto & s = this->ConvectionDiffusion<d>::systems.front();
  std::vector<const FEFunction*> feFunctionPtr = {&s.fe_function};
  
  int afc_ini_supg = 0;
  if(!this->db["algebraic_flux_correction"].is("none") && iteration == 0 
    && this->db["afc_initial_iterate"].is("supg"))
  {
    afc_ini_supg = 1;   
    Output::print<2>("initial solution is SUPG solution");
    this->db["space_discretization_type"].set("supg");
  }
  
  LocalAssembling<d> laObject(this->db, laTypet, feFunctionPtr,
                                this->example.get_coeffs());
  
  //fetch stiffness matrix as block
  auto blocks = s.matrix.get_blocks_uniquely();
  SquareMatrixD* block = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
  
  /**
   * @brief: For Fixed Point RHS(FPR), the matrix and the RHS has to assembled 
   * only once. Becasue in FPR the rhs gets modified in 
   * do_algebraic_flux_correction, where the modification is only added and
   * hence we copy the old rhs.
   * 
   */
  if(is_not_afc_fixed_point_rhs == 1)  
  {
    bool do_afc = !this->db["algebraic_flux_correction"].is("none");
    ConvectionDiffusion<d>::call_assembling_routine(s, laObject, do_afc);
    // only important for MCL Steady state
    // @brief : We replace the right side by solution of the Poisson problem.
    //          This is done in AlgebrainFluxCorrection.C.
    // This is for homogeonous Neumann BC.
    // NOTE: If the Neumann BC is not zero, then only the RHS contribution is set to zero
    //      the Neumann condition is added.
    // If iteration == 0 and afc_initial_iterate is set to SUPG (or upwind), then one does not 
    // need to reset the RHS.
    if ((limiter == AlgebraicFluxCorrection::Limiter::MONOLITHIC_STEADY) && iteration != 0)
    {
      s.rhs.ResetActive();
    }
    if(iteration == 1)
      rhs_copy=s.rhs;    
  }
  //for FIXED_POINT_RHS and iteration>1
  else
  {
    s.rhs=rhs_copy;  
  }
  if(is_not_afc_fixed_point_rhs == 1)
  {
    if (!this->db["algebraic_flux_correction"].is("none")&&(iteration == 0)&&
          this->db["afc_initial_iterate"].is("upwind"))
    {
      Output::print<2>("upwind for convection-diffusion equation");
#ifdef __2D__
      UpwindForConvDiff(laObject.GetCoeffFct(), block, s.rhs.get_entries(),
                        s.fe_space.get(),nullptr, nullptr, false);
#else
      UpwindForConvDiff(block, &s.rhs.get_entries_vector()[0], 
                          s.fe_space.get(),laObject);  
#endif  
    }  
  }
  if (afc_ini_supg == 1)
    this->db["space_discretization_type"].set("galerkin");
  
  if (iteration == 0)
  {
   this->db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    /*// the following flags are for the BAIL proceedings
    // switch to fixed point rhs
    if (((int)db["afc_nonlinloop_switch_to_newton_scheme"]==10)||
        ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==20))
    {
      db["afc_fixed_point_matrix_weight"] = 0.0;
      db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    // switch to fixed point matrix 
    if (((int)db["afc_nonlinloop_switch_to_newton_scheme"]==11)||
       ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==21)) 
    {
      db["afc_fixed_point_matrix_weight"] = 1.0;
      db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }*/
  }

  if(!this->db["algebraic_flux_correction"].is("none"))
    do_algebraic_flux_correction(iteration, is_not_afc_fixed_point_rhs);
}

/* ************************************************************************ */
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(
                                          std::string afc_iteration_scheme);
/* ************************************************************************ */
template<int d>
bool ConvectionDiffusion_AFC<d>::solve(const int iteration)
{
  double t = GetTime();
  auto& s = this->systems.front();
  #ifdef _MPI
  const TParFECommunicator3D& comm = s.fe_space->get_communicator();
#endif
  BlockVector res = s.rhs;
  BlockVector current_rhs = s.rhs;
  BlockVector current_sol = s.solution;
  AlgebraicFluxCorrection::Iteration_Scheme it_scheme =
                       string_to_it_scheme(this->db["afc_iteration_scheme"]);
  
  
  /*
   * In case of not FPR, after the first iteration this has to be zero as then
   * we have to assemble the matrix again and again.
   */
                       
  is_not_afc_fixed_point_rhs=0;

  // To assemble for 1st iteration with FPR
  if(it_scheme != AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS
    || iteration == 1)
    is_not_afc_fixed_point_rhs=1;

  // Special treatment of first iteration
  if (iteration == 0)
  {
    // Compute residual vector
    s.matrix.BlockMatrix::apply_scaled_add(s.solution , res ,-1.0);

    // Compute the norm of the residual vector, residual_old 
    //  is the norm of residual r_k
#ifdef _MPI
    residual_old = res.norm({&comm});
#else
    residual_old = res.norm();
#endif
    old_solution = s.solution;
#ifdef _MPI
    if(this->ConvectionDiffusion<d>::solver.get_db()["solver_type"].is("direct"))
    {
      MumpsWrapper mumps_wrapper(s.matrix);
      mumps_wrapper.solve(s.rhs, s.solution);
    }
    else
      // same as sequential
    this->ConvectionDiffusion<d>::solver.solve(s.matrix, s.rhs, s.solution);
#else
    this->ConvectionDiffusion<d>::solver.solve(s.matrix, s.rhs, s.solution);
#endif
    
    initial_solution = s.solution;
    
    t = GetTime() - t;
    Output::print(" solving of a Convection Diffusion problem done in ",
                   t, " seconds");
    double h_min, h_max;
    
#ifdef __2D__
    const TFESpace2D & space = *this->systems.front().fe_space;
#else
    const TFESpace3D & space = *this->systems.front().fe_space;
#endif
    auto coll = space.GetCollection();
    coll->GetHminHmax(&h_min, &h_max);
    
    // For reqularized Newton method, the sigma is depending on h. Idea from
    // [BB17.CMAME]
    double sigma = (double)this->db["afc_newton_regu_sigma"];
    sigma *= h_max * h_max * h_max * h_max;
    this->db["afc_newton_regu_sigma"] = sigma;
    Output::print<3>("afc_newton_regu_sigma changed to ", 
                                    this->db["afc_newton_regu_sigma"]);   
    
    // The stopping tolerance is dependent on mesh refinment
    double dof = space.get_n_dof();
#ifdef _MPI
    dof = space.get_communicator().get_n_global_dof();
#endif
    double eps = (double)this->db["afc_nonlinloop_epsilon"];
    eps *= std::sqrt(dof);
    this->db["afc_nonlinloop_epsilon"] = eps;
    Output::print<2>("afc_nonlinloop_epsilon normalized to ", 
                     this->db["afc_nonlinloop_epsilon"]);
    Output::print<2>("nonlinear step ", iteration, " residual: ", residual_old);
    return(false);
  }

  if(!this->db["algebraic_flux_correction"].is("none"))
  {
    //Anderson acceleration
    if ((this->db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >=
      (int)this->db["afc_nonlinloop_anderson_acc_start"]))
    {
      int N_Unknowns = res.length();
      anderson_acceleration_damping(N_Unknowns,iteration, 
                                    solAnderson, deltaAnderson);
    }

    //Constant Damping
    if (this->db["afc_nonlinloop_damping_factor_constant"].is("yes"))
    {
      // compute proposal for the next solution
      FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
      AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, 
                                                       s.solution, this->db, one_block);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs;
      s.matrix.BlockMatrix::apply_scaled_add(s.solution , res ,-1.0);
      // compute the norm of the residual vector, residual_old is the 
      // norm of residual r_k
#ifdef _MPI
      residual = res.norm({&comm});
#else
      residual = res.norm();
#endif
      Output::print<4>("  residual for proposed new iterate ", residual);
    }
    
    //Dynamic Damping from [JK08.CMAME]
    else
      dynamic_damping(iteration);
  }

  old_solution.add_scaled(s.solution,-1.0);

  Output::print<2>("nonlinear step ", iteration, " residual: ", residual, 
                   " reduction ",residual/residual_old ," change of sol ", 
#ifdef _MPI
                    old_solution.norm({&comm})
#else
                    old_solution.norm()
#endif      
);
  
  // stopping criterion satisfied
  if ((residual < (double)this->db["afc_nonlinloop_epsilon"])&&(iteration >1))
  {
    time_total = GetTime() - time_total;
    Output::print<1>("NONLINEAR ITERATION: ite ", iteration, " res ", 
                     residual, " rejections ",
                     rejected_steps, " time ", time_total, " t/it ", 
                     time_total/(iteration+rejected_steps));
    return(true);
  }
  //maximal number of iterations
  if (iteration == (int)this->db["afc_nonlinloop_maxit"])
  {
    time_total = GetTime() - time_total;
    Output::print<1>("MAX_NONLINEAR ITERATION: ite ", iteration, " res ", 
                     residual, " rejections ", rejected_steps, " time ",
                     time_total, " t/it ", time_total/(iteration+rejected_steps));
    return(true);
  }
  // storage of old solution
  old_solution =s.solution;

  // Solve the linear system
  //Matrix is factorised only once and then stored in the system
#ifdef _MPI
    if(this->ConvectionDiffusion<d>::solver.get_db()["solver_type"].is("direct"))
    {
      MumpsWrapper mumps_wrapper(s.matrix);
      mumps_wrapper.solve(s.rhs, s.solution);
    }
    else
    {
      if(is_not_afc_fixed_point_rhs == 1)
       this->ConvectionDiffusion<d>::solver.update_matrix(s.matrix);
      this->ConvectionDiffusion<d>::solver.solve(s.rhs, s.solution);
    }
#else
  if(is_not_afc_fixed_point_rhs == 1)
   this->ConvectionDiffusion<d>::solver.update_matrix(s.matrix);
  this->ConvectionDiffusion<d>::solver.solve(s.rhs, s.solution);
#endif
  //Previous Implementation
  //this->solver.solve(s.matrix_, s.rhs_, s.solution_);
  
  t = GetTime() - t;
  Output::print<4>("  iteration done in ", t, " seconds");

  //Dynamic choice of Newton damping parameter
  if ((this->db["afc_iteration_scheme"].is("newton"))||
       (this->db["afc_iteration_scheme"].is("newton_regu")))
  {
    double threshold = (double)this->db["afc_change_method_threshold"];
    if (residual <= threshold)
    {
      double omega_derivative = 
                 (double)this->db["afc_fixed_point_derivative_weight_factor"];
      // first application
      if (omega_derivative < 1e-3)
        omega_derivative = 0.25;
      if (residual/residual_old > 0.99)
      {
        omega_derivative *=0.999;
        if (omega_derivative < 0.1)
          omega_derivative = 0.1;
      }
      if (residual/residual_old < 0.99)
      {
        omega_derivative *=1.001;
        if (omega_derivative >1.0)
          omega_derivative = 1.0;
      }
      this->db["afc_fixed_point_derivative_weight_factor"] = omega_derivative;
    }
    // go back to former method
    if (residual > 100*threshold)
    {
      this->db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    Output::print<2>("afc_fixed_point_derivative_weight_factor ",
		     (double)this->db["afc_fixed_point_derivative_weight_factor"]);
  }

  /* ************************************************************************ */
  // this scheme is for the BAIL proceedings
  /* ************************************************************************ */
  /*if ((db["afc_iteration_scheme"].is("newton"))&&
     ((int)db["afc_nonlinloop_switch_to_newton_scheme"]>=10))
  {
    double threshold = (double)db["afc_change_method_threshold"];
    // switch to formal Newton
    if (residual <= threshold)
    {
      db["afc_fixed_point_matrix_weight"] = 1.0;
      db["afc_fixed_point_derivative_weight_factor"] = 1.0;
    }
    else
    {
        if (residual > 10*threshold)
    {
        db["afc_fixed_point_matrix_weight"] = 0.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    }
    if (residual > 100*threshold)
    {
      if ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==20)
      {
        db["afc_fixed_point_matrix_weight"] = 0.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
      }
      if ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==21)
      {
        db["afc_fixed_point_matrix_weight"] = 1.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
      }
    }
    Output::print<2>("afc_fixed_point_matrix_weight ",  
                          (double)db["afc_fixed_point_matrix_weight"]);
  }*/
  
  //Storage of old residual
  residual_old=residual;
  return(false);
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::dynamic_damping(const int iteration)
{
  auto& s = this->systems.front();
#ifdef _MPI
  const TParFECommunicator3D& comm = s.fe_space->get_communicator();
#endif
  double omega = this->db["afc_nonlinloop_damping_factor"];
  int first_damp=1;
  BlockVector res = s.rhs;
  BlockVector current_rhs = s.rhs;
  BlockVector current_sol = s.solution;
  FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();

  while(1)
  {
    
    AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution, 
                                                       this->db, one_block);
    assemble(iteration);
    // calculation of residual r_k+1
    res = s.rhs;
    s.matrix.BlockMatrix::apply_scaled_add(s.solution , res ,-1.0);
    // compute the norm of the residual vector, 
    //    residual_old is the norm of residual r_k
#ifdef _MPI
    residual = res.norm({&comm});
#else
    residual = res.norm();
#endif
    Output::print<4>("  residual for proposed new iterate ", residual);
    // accept the first damping parameter
    if (iteration == 1)
      break;
    
    /*
     * if the norm of the residual vector decreases or if the damping parameter 
     * is already very small then accept the new iterate
     */
    if (residual < residual_old ||
               omega <= (double)this->db["afc_nonlinloop_damping_factor_min_tol"]*
      (double)this->db["afc_nonlinloop_damping_factor_min"])
    {
      /*
       * if the norm of the residual decreases without having decreased the
       * damping parameter in this step, then increase the damping parameter 
       * if possible
       */
      if(residual < residual_old && first_damp == 1)
      {
        this->db["afc_nonlinloop_damping_factor_max"]=
        std::min((double)this->db["afc_nonlinloop_damping_factor_max_global"],
          (double)this->db["afc_nonlinloop_damping_factor_max_increase"]
          *(double)this->db["afc_nonlinloop_damping_factor_max"]);
        omega=std::min((double)this->db["afc_nonlinloop_damping_factor_max"],
        (double)this->db["afc_nonlinloop_damping_factor_increase"]*omega);
      }
      Output::print<2>("  iterate accepted, damping factor ", omega, " ",
          this->db["afc_nonlinloop_damping_factor_max"]);
      this->db["afc_nonlinloop_damping_factor"] = omega;
      break;
    }
    else
    {
      if (this->db["afc_nonlinloop_anderson_acc"].is("no"))
      {
        rejected_steps++;
        // get starting situation back
        s.solution = current_sol;
        s.rhs = current_rhs;
      }
      // reduce damping factor
      omega = std::max((double)this->db["afc_nonlinloop_damping_factor_min"],
        omega*(double)this->db["afc_nonlinloop_damping_factor_decrease"]);
      // reduce maximal value for damping factor
      if(first_damp==1)
      {
        this->db["afc_nonlinloop_damping_factor_max"] =
        std::max((double)this->db["afc_nonlinloop_damping_factor_min"],
          (double)this->db["afc_nonlinloop_damping_factor_max_decrease"]
              *(double)this->db["afc_nonlinloop_damping_factor_max"]);
        first_damp = 0;
      }
      Output::print<2>("  iterate rejected, res old ", residual_old,
          " res new: ", residual, " damping factor :", omega);
      this->db["afc_nonlinloop_damping_factor"] = omega;
    }
    // if Anderson acceleration, then only reduction of the damping parameter
    // but no computation of different update
    if ((this->db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >=
      (int)this->db["afc_nonlinloop_anderson_acc_start"]))
      break;
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::anderson_acceleration_damping(
  int N_Unknowns, const int iteration,
  std::list<std::vector<double>> & solAnderson,
  std::list<std::vector<double>> & deltaAnderson)
{
  auto& s = this->systems.front();
  int k;

  std::vector <double> newSol(N_Unknowns);
  std::vector <double> newDelta(N_Unknowns);

  int lwork, ierr;
  lwork=-1;
  double temp=0.;

  int ndim, kdim;
  double *residuals, *fpast, *work;
  int kAnderson = (int)this->db["afc_nonlinloop_anderson_acc_vec"];
  BlockVector update, new_solution;
  update = s.solution;

  // update =\hat u^{\nu+1}-\tilde{u^{\nu}}
  update.add_scaled(old_solution,-1.0);

  //collect data for previous iterations
  for(int i = 0 ; i < N_Unknowns ; i++)
  {
    newDelta[i] = s.solution[i];
    newSol[i] = update[i];
  }
  //update list deltaAnderson
  //add new element at the end
  deltaAnderson.push_back(newDelta);
  //if list has sufficiently many entries, remove first entry
  if((int)deltaAnderson.size()>kAnderson)
    deltaAnderson.pop_front();

  //update list solAnderson
  solAnderson.push_back(newSol);
  if((int)solAnderson.size() > kAnderson)
    solAnderson.pop_front();

  //Use Anderson acceleration if there are sufficiently many entries
  if((iteration>kAnderson+(int)this->db["afc_nonlinloop_anderson_acc_start"]) &&
    ((int)solAnderson.size() == kAnderson))
  {
    ndim = N_Unknowns;
    kdim = kAnderson;
    //entries for solution
    residuals = new double[2*ndim*kdim];
    //entries for updates
    fpast = residuals+ndim*kdim;

    //copy arrays
    std::list <std::vector<double>>::iterator i_delta = deltaAnderson.begin();
    std::list <std::vector<double>>::iterator i_sol = solAnderson.begin();
    for(int i = 0 ; i < kdim ; i++)
    {
      //advance list iterators to next element
      for(int j = 0 ; j < ndim ; j++)
      {
        if(i < kdim)
        {
          k = i*ndim+j;
          fpast[k] = i_delta->at(j);
          residuals[k] = i_sol->at(j);
        }
      }
      std::advance(i_delta,1);
      std::advance(i_sol,1);
    }                                             //end of loop i
    //set initial length of work array=-1
    if(lwork == -1)
    {
      //the first call sets the length of work
      //anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
      //		    &temp, lwork, &ierr);
      anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
        &temp, lwork, &ierr, &alphas_x_i[0]);
    }
    lwork=(int)temp;
    work=new double[lwork];
    Output::print<2>("Call anderson acceleration, iteration=", 
		     iteration-kAnderson-(int)this->db["afc_nonlinloop_anderson_acc_start"]);
    anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
      work, lwork, &ierr, &alphas_x_i[0]);
    //s.solution = new_solution;
    delete[] work;
    delete[] residuals;
  }
}
/*******************************************************************/
/** Compute solution at x=0.5 and y=0.5 
 * Use to find thickness of layer 
 * Currently set for example from HMM1986. 
********************************************************************/
template<int d>
double ConvectionDiffusion_AFC<d>::ComputeCutLines()
{
#ifdef __2D__
  TFEFunction2D ufct = ConvectionDiffusion<d>::get_function();
  // coordinates of the domain
  double x_coord = 1 , y_coord = 1;
  // coordinates of the layer
  double x1 = 0, x2 = 0;
  // solution
  double u_sol[3];
  double h;
  //Linear Interpolation of the layer points
  //x1_h : x1 - h
  //x2_h : x2 - h
  double x1_h = 0.0 , x2_h = 0.0;
  //y1_h : Solution at x1-h
  //y2_h : Solution at x2-h
  //y1   : Solution at x1
  //y2   : Solution at x2
  double y1_h = 0.0 , y2_h = 0.0, y1 = 0.0, y2 = 0.0, old_sol = 0.0;
  //x_start : Starting point of the internal layer
  //x_end   : End point of the internal layer
  double x_start = 0.0, x_end = 0.0;
  // Bound points can be chosen arbitrary but choose a large
  // number to get more points
  int bound_points = 10001, x_init, x_final;
  
  h = 1.0/(bound_points-1);

  // smearing of the boundary layer at y = 0.25
  y_coord = 0.25;
  // Compute the solution between the strip (x_init, 0.25)
  // and (x_final, 0.25)
  x_init = (int) (0.175/h);
  x_final = (int) ((0.75)/h);
  for (int i = x_init ; i <= x_final ; i++)
  {
      x_coord = i*h;
      ufct.FindGradient(x_coord ,y_coord ,u_sol);
      if ((u_sol[0] >= 0.1) && (x1 == 0))
      {
        x1 = x_coord;
        //Store the solution for linear interpolation
        y1 = u_sol[0];
        y1_h = old_sol;
      } 
      if ((u_sol[0] >= 0.9) && (x2 == 0))
      {
        x2 = x_coord;
        //Store the solution for linear interpolation
        y2 = u_sol[0];
        y2_h = old_sol;
        break;
      }
      old_sol =  u_sol[0];
  }

  x1_h = x1 - h;
  x2_h = x2 - h;
  //Linear interpolation of point (x1_h, y1_h) and (x1,y1)
  // x = ((y-y1_h)*(x1-x1_h))/(y1-y1_h) + x1_h
  x_start = ((0.1 - y1_h)*(x1 - x1_h))/(y1 - y1_h) + x1_h;
  x_end   = ((0.9 - y2_h)*(x2 - x2_h))/(y2 - y2_h) + x2_h;
  //Output::print<1>("DEBUG: x1 = ", x1, "  x2 = ", x2);
  // Width of the layer
  double value = x_end - x_start;
  return value;  
#else
  Output::print<1>("Cut Lines not available in 3D");
  return 0.0;
#endif
}

/*******************************************************************/
/** Compute solution at (x,y=1.0) 
 * Currently set for example from Hemker problem.
 * This kind of computation is also done in the example file but
 * that code doesn't work with adaptive grids. This computation
 * is computationally inefficient but it works 
********************************************************************/
template<int d>
double ConvectionDiffusion_AFC<d>::ComputeCutLineY()
{
#ifdef __2D__
  TFEFunction2D ufct = ConvectionDiffusion<d>::get_function();
  // coordinates of the domain
  double x_coord = 1 , y_coord = 1;
  // solution
  double u_sol[3];
  double h;
  // Bound points can be chosen arbitrary but choose a large
  // number to get more points
  int bound_points = 2001, x_init, x_final;
  
  h = 1.0/(bound_points-1);

  // smearing of the boundary layer at y = 0.25
  y_coord = 1.0;
  // Compute the solution between the strip (x_init, 0.25)
  // and (x_final, 0.25)
  x_init = (int) (-2.0/h);
  x_final = (int) (8.0/h);
  for (int i = x_init ; i <= x_final ; i++)
  {
      x_coord = i*h;
      ufct.FindGradient(x_coord ,y_coord ,u_sol);
      Output::print("x coordinate : ", x_coord, " ,u_value : ", u_sol[0]);
  }
  return 0.0;  
#else
  Output::print<1>("Cut Lines not available in 3D");
  return 0.0;
#endif
}

/*******************************************************************/
/** Compute solution at (x=4.0,y) 
 * Currently set for example from Hemker. 
********************************************************************/
template<int d>
double ConvectionDiffusion_AFC<d>::ComputeCutLineX()
{
#ifdef __2D__
  TFEFunction2D ufct = ConvectionDiffusion<d>::get_function();
  // coordinates of the domain
  double x_coord = 1 , y_coord = 1;
  // solution
  double u_sol[3];
  double h;
  // Bound points can be chosen arbitrary but choose a large
  // number to get more points
  int bound_points = 10001, y_init, y_final;
  
  h = 1.0/(bound_points-1);

  // smearing of the boundary layer at y = 0.25
  x_coord = 4.0;
  y_init = (int) (-3.0/h);
  y_final = (int) (3.0/h);
  for (int i = y_init ; i <= y_final ; i++)
  {
      y_coord = i*h;
      ufct.FindGradient(x_coord ,y_coord ,u_sol);
      Output::print("y coordinate : ", y_coord, " ,u_value : ", u_sol[0]);
  }
  return 0.0;  
#else
  Output::print<1>("Cut Lines not available in 3D");
  return 0.0;
#endif
}

/*******************************************************************/
/** Use to find thickness of layer 
 * Currently set for example from Hemker. 
********************************************************************/
template<int d>
double ConvectionDiffusion_AFC<d>::SmearX()
{
#ifdef __2D__
  TFEFunction2D ufct = ConvectionDiffusion<d>::get_function();
  // coordinates of the domain
  double x_coord = 0.0 , y_coord = 0.0;
  // coordinates of the layer
  double y1 = 0, y2 = 0;
  // solution
  double u_sol[3];
  double h;
  //Linear Interpolation of the layer points
  //y1_h : y1 - h
  //y2_h : y2 - h
  double y1_h = 0.0 , y2_h = 0.0;
  //x1_h : Solution at y1-h
  //x2_h : Solution at y2-h
  //x1   : Solution at y1
  //x2   : Solution at y2
  double x1_h = 0.0 , x2_h = 0.0, x1 = 0.0, x2 = 0.0, old_sol = 0.0;
  //y_start : Starting point of the internal layer
  //y_end   : End point of the internal layer
  double y_start = 0.0, y_end = 0.0;
  // Bound points can be chosen arbitrary but choose a large
  // number to get more points
  int bound_points = 10001, y_init, y_final;
  
  h = 1.0/(bound_points-1);

  // smearing of the boundary layer at x = 4.0
  x_coord = 4.0;
  // Compute the solution between the strip (4.0, y_init)
  // and (4.0, y_final)
  y_init = (int) ((-3.0)/h);
  y_final = (int) ((3.0)/h);
  for (int i = y_init ; i <= y_final ; i++)
  {
      y_coord = i*h;
      ufct.FindGradient(x_coord ,y_coord ,u_sol);
      if ((u_sol[0] >= 0.1) && (y1 == 0))
      {
        y1 = y_coord;
        //Store the solution for linear interpolation
        x1 = u_sol[0];
        x1_h = old_sol;
      } 
      if ((u_sol[0] >= 0.9) && (y2 == 0))
      {
        y2 = y_coord;
        //Store the solution for linear interpolation
        x2 = u_sol[0];
        x2_h = old_sol;
        break;
      }
      old_sol =  u_sol[0];
  }

  y1_h = y1 - h;
  y2_h = y2 - h;
  //Linear interpolation of point (x1_h, y1_h) and (x1,y1)
  // x = ((y-y1_h)*(x1-x1_h))/(y1-y1_h) + x1_h
  y_start = ((0.1 - x1_h)*(y1 - y1_h))/(x1 - x1_h) + y1_h;
  y_end   = ((0.9 - x2_h)*(y2 - y2_h))/(x2 - x2_h) + y2_h;
  //Output::print<1>("DEBUG: x1 = ", x1, "  x2 = ", x2);
  // Width of the layer
  double value = y_end - y_start;
  return value;  
#else
  Output::print<1>("Cut Lines not available in 3D");
  return 0.0;
#endif
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::output(int i)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  
  ConvectionDiffusion<d>::output(i);
  
  double calculated_afc_error = afc_error();
  
#ifdef _MPI
  double thread_afc_error = calculated_afc_error;
  
  MPI_Reduce(&thread_afc_error,
             &calculated_afc_error, 1,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  
  if(this->db["output_compute_errors"])
  {
    if(my_rank == 0)
    {
        //AFC Errors
        if(this->db["algebraic_flux_correction"].is("afc"))
        {
        //NOTE: we assume sigma_0 to be 1. If the problem doesn't have reaction
        //      then change accordingly
        double epsilon = this->ConvectionDiffusion<d>::example.get_nu();
        Output::print<4>("sigma*L2     : ", 
                        setprecision(14), ConvectionDiffusion<d>::errors[0]);
        Output::print<4>("epsilon*H1-semi: ", 
                        setprecision(14), std::sqrt(epsilon)*
                        ConvectionDiffusion<d>::errors[1]);     
        Output::print<4>("Interpolation Energy Norm: ", setprecision(14), 
                        std::sqrt(interpolation_errors[0]*
                        interpolation_errors[0]
                        +epsilon*interpolation_errors[1]*
                        interpolation_errors[1]));    
        Output::print<4>("Interpolation SUPG_Norm: ", setprecision(14), 
                        interpolation_errors[2]);    
        double val = std::sqrt(ConvectionDiffusion<d>::errors[0]*
                        ConvectionDiffusion<d>::errors[0]
                        +epsilon*ConvectionDiffusion<d>::errors[1]*
                        ConvectionDiffusion<d>::errors[1]);
        Output::print<1>("Energy : ", setprecision(14), val);
        Output::print<1>("AFC    : ", setprecision(14), 
                        std::sqrt(val*val+calculated_afc_error));
//                         std::sqrt(val*val+afc_error()));
        //Output effectivity index only for a posteriori problems
        //if((int)TDatabase::ParamDB->INTERNAL_COERCIVITY > 0)
            Output::print<2>("eta-eff: ",  
                            TDatabase::ParamDB->INTERNAL_COERCIVITY/val);
        }
    }
  }
  //Thickness of internal layer
  if(this->db["compute_cut_lines"].is("yes"))
    Output::print<1>("smear  : ", SmearX());
}

/* ************************************************************************* */
template<int d>
double ConvectionDiffusion_AFC<d>::afc_error()
{
  auto & s = this->ConvectionDiffusion<d>::systems.front();
  FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
#ifdef _MPI
  const TParFECommunicator3D& comm = one_block.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  const int nDOFs = one_block.get_n_rows();
  const int* ColInd = one_block.get_vector_columns();
  const int* RowPtr = one_block.get_row_ptr();
  double* afc_matrix_D_entries = afc_matrix_D.GetEntries();
  double error = 0.0;
  int index;
  for(int i = 0 ; i < nDOFs ; i++)  
  {
#ifdef _MPI
    //skip the non-master rows in MPI case
    if(masters[i] != rank)
      continue;
#endif
    for(int j = RowPtr[i] ; j < RowPtr[i+1] ; j++)
    {
      index = ColInd[j];
      error += (1-afc_alphas[j])*afc_matrix_D_entries[j]*(s.solution[index]
                -exact_interpolant[index] - s.solution[i] 
                + exact_interpolant[i])*(s.solution[i]-exact_interpolant[i]);
    }
  }
  return error;
}

/* ************************************************************************* */
AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter);
/* ************************************************************************* */

template<int d>
void ConvectionDiffusion_AFC<d>::do_algebraic_flux_correction(
  const int iteration, const int is_not_afc_fixed_point_rhs)
{
  Output::print<4>("AFC: enter do_algebraic_flux_correction");
  auto & s = this->ConvectionDiffusion<d>::systems.front();
  bool compute_D_and_gamma = false;
  
  //Get the bilinear coeffecients for monolithic steady limiter
  //coeff[0] : Diffusion, coeff[1] : Convective field in x-direction, coeff[2] : Convective field in y-direction,
  //coeff[3] : Convective field in z-direction (if 3D) otherwise reaction coefficient.
#ifdef __2D__
  const CoeffFct2D bilinear_coeffs = this->ConvectionDiffusion<d>::example.get_coeffs();
#elif __3D__
  const CoeffFct3D bilinear_coeffs = this->ConvectionDiffusion<d>::example.get_coeffs();
#endif

  double x =0.0, y= 0.0, *coeffs;
  coeffs = new double[13];
#ifdef __2D__
  bilinear_coeffs(1, &x, &y, nullptr, &coeffs);
#elif __3D__
  double z = 0.0;
  bilinear_coeffs(1, &x, &y, &z, nullptr, &coeffs);
#endif
  
  //determine which kind of afc to use
  if(this->db["algebraic_flux_correction"].is("default") ||
      this->db["algebraic_flux_correction"].is("afc"))
  {
    // Determine which kind of limiter to use
    AlgebraicFluxCorrection::Limiter limiter = 
                                    string_to_limiter(this->db["afc_limiter"]);
    // Determine which kind of iteration scheme is active
    AlgebraicFluxCorrection::Iteration_Scheme it_scheme = 
                          string_to_it_scheme(this->db["afc_iteration_scheme"]);
    
    //get pointers/references to the relevant objects
    auto& feSpace = *s.fe_space;
    FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
#ifdef _MPI
    // NOTE: Taken from previous implementation in TCD FEM-FCT
    const TParFECommunicator3D& comm = s.fe_space->get_communicator();
    // solution is needed in consistency 2 for afc.
    comm.consistency_update(s.solution.block(0), 2);
#endif
    const std::vector<double>& solEntries = s.solution.get_entries_vector();
    std::vector<double>& rhsEntries = s.rhs.get_entries_vector();
    
    // fill a vector "neumannToDirichlet" with those rows that got
    // internally treated as Neumann although they are Dirichlet
    int firstDiriDof = feSpace.get_n_active();
    int nDiri = feSpace.get_n_dirichlet();
    
    std::vector<int> neumToDiri(nDiri, 0);
    std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);
    
    //AFC only applied after the 0th iteration or the initial iterate is zero.
    if (iteration > 0 || this->db["afc_initial_iterate"].is("afc_zero"))
    {
      /*
       * In fixed_point_rhs after the first iteration the matrix A is computed 
       * and stored, so is D and A+D(see matrix_copy). Hence, one doesn't need 
       * to compute D and gamma again (gamma doesn't depend on A or D) after the
       * first iteration. For other methods, A is getting assembled again and
       * hence D needs to be assembled again.
       */
      
      if(is_not_afc_fixed_point_rhs)
      {
        compute_D_and_gamma = true;  
        Output::print<4>("Reseting matrix D in AFC");
        afc_matrix_D.reset();
      }
      if ((limiter== AlgebraicFluxCorrection::Limiter::BJK17) 
                          && (afc_gamma.empty()))
      {
        Output::print<4>("AFC: vector gamma");
        afc_gamma.resize(feSpace.get_n_dof(),0.0);  
      }
      
      // apply AFC
      if(is_not_afc_fixed_point_rhs == 1)
      {
        AlgebraicFluxCorrection::steady_state_algorithm(
          one_block, diffusion_matrix,
          solEntries,rhsEntries,poisson_sol,
          neumToDiri, coeffs,
          afc_matrix_D, afc_matrix_D_B,
          afc_gamma, afc_alphas ,compute_D_and_gamma, this->db,
          limiter, it_scheme, is_not_afc_fixed_point_rhs);
        //performed only once in the whole iteration process
        if(iteration == 1 && 
           it_scheme == AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS)
        {
          //matrix_copy=A+D
          //In AlgebraicFluxCorrection.C we don't need to add D again
          matrix_copy = s.matrix;  
        }  
      }
      //case for fixed point rhs and iteration>1
      else
      {
        /*
         * the matrix that is used for the AFC scheme needs to have the correct 
         * Dirichlet entries and hence after the first iteration sending these 
         * values.
         */
        FEMatrix& one_block1 = *matrix_copy.get_blocks_uniquely().at(0).get();
        AlgebraicFluxCorrection::steady_state_algorithm(
          one_block1, diffusion_matrix,
          solEntries,rhsEntries,poisson_sol,
          neumToDiri, coeffs,
          afc_matrix_D, afc_matrix_D_B,
          afc_gamma, afc_alphas ,compute_D_and_gamma, this->db,
          limiter, it_scheme, is_not_afc_fixed_point_rhs);  
      }
      // Previous Implementation
      /*AlgebraicFluxCorrection::steady_state_algorithm(
       *         one_block,
       *         solEntries,rhsEntries,
       *         neumToDiri,
       *        s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma, db,
       *        limiter, it_scheme, is_not_afc_fixed_point_rhs);
       *        Output::print<1>("Matrix Norm: ",
       *                  s.matrix_.get_blocks_uniquely({{0,0}})[0]->GetNorm());
       */  
    }    
    //set the Dirichlet rows and Hanging rows correctly
    if (is_not_afc_fixed_point_rhs)
      AlgebraicFluxCorrection::correct_dirichlet_hanging_rows(one_block);    
    //...and in the right hand side, too, assume correct in solution vector
    s.rhs.copy_nonactive(s.solution);
    AlgebraicFluxCorrection::correct_dirichlet_hanging_rhs(one_block, s.rhs);
    
#ifdef _MPI
    // NOTE: Taken from previous implementation in TCD FEM-FCT
    // this update is necessary here, because copy_nonactive()
    // can disturb the consistency of a vector,
    // see Time_NSE3D::assemble_rhs() for more details
    comm.consistency_update(s.rhs.block(0), 2);
#endif
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme is unknown "
              "to class ConvectionDiffusion<",d,">.");  
  }
}


/* ************************************************************************* */

AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter)
{
  if (afc_limiter == std::string("kuzmin"))
    return AlgebraicFluxCorrection::Limiter::KUZMIN;
  else if (afc_limiter == std::string("BJK17"))
    return AlgebraicFluxCorrection::Limiter::BJK17;
  else if (afc_limiter == std::string("monolithic"))
    return AlgebraicFluxCorrection::Limiter::MONOLITHIC;
  else if (afc_limiter == std::string("monolithic_steady"))
    return AlgebraicFluxCorrection::Limiter::MONOLITHIC_STEADY;
  else if (afc_limiter == std::string("modified_kuzmin"))
    return AlgebraicFluxCorrection::Limiter::MUAS;
  else if (afc_limiter == std::string("MUAS"))
    return AlgebraicFluxCorrection::Limiter::MUAS;
  else if (afc_limiter == std::string("MUAS_Kno21"))
    return AlgebraicFluxCorrection::Limiter::MUAS_Kno21;
  else if (afc_limiter == std::string("MUAS_MAX"))
    return AlgebraicFluxCorrection::Limiter::MUAS_MAX;
  else if (afc_limiter == std::string("MUAS_MAX_ABS"))
    return AlgebraicFluxCorrection::Limiter::MUAS_MAX_ABS;
  else
  {
    ErrThrow("afc_limiter ", afc_limiter, " not implemented!!!");
  }
}

/* ************************************************************************* */
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(
                                              std::string afc_iteration_scheme)
{
  if (afc_iteration_scheme == std::string("fixed_point_rhs"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS;
  else if (afc_iteration_scheme == std::string("fixed_point_matrix"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_MATRIX;
  else if (afc_iteration_scheme == std::string("newton"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
  else if (afc_iteration_scheme == std::string("newton_no_damp"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
  else if (afc_iteration_scheme == std::string("newton_regu"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON_REGU;
  else
  {
    ErrThrow("afc_iteration_scheme",afc_iteration_scheme, "not implemented!!!");
  }
}

/* ************************************************************************* */
template <int d>
const std::vector< double >& ConvectionDiffusion_AFC<d>::get_afc_alphas() const
{
  return afc_alphas;
}

/* ************************************************************************* */
template <int d>
const FEMatrix& ConvectionDiffusion_AFC<d>::get_afc_D_entries() const
{
  AlgebraicFluxCorrection::Limiter limiter =
                       string_to_limiter(this->db["afc_limiter"]);
  //For modified Kuzmin the matrix is D-B 
  //Used specifically for a Posteriori error estimator
  if ((limiter == AlgebraicFluxCorrection::Limiter::MUAS)||
      (limiter == AlgebraicFluxCorrection::Limiter::MUAS_Kno21)||
      (limiter == AlgebraicFluxCorrection::Limiter::MUAS_MAX)||
      (limiter == AlgebraicFluxCorrection::Limiter::MUAS_MAX_ABS))
    return afc_matrix_D_B;
  else
    return afc_matrix_D;
}

/* ************************************************************************* */

template<int d> 
void ConvectionDiffusion_AFC<d>::assemble_poisson()
{
  LocalAssembling_type laType = LocalAssembling_type::Poisson;
  for(auto & s : this->ConvectionDiffusion<d>::systems)
  {
    std::vector<const FEFunction*> feFunctionPtr = {&s.fe_function};
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling<d> laObject(this->db, laType, feFunctionPtr,
                                this->ConvectionDiffusion<d>::example.get_coeffs());
    //This is pure Neumann problem and hence the matrix will be set to the correct one
    this->ConvectionDiffusion<d>::call_assembling_routine(s, laObject, false);
    //Store this for extracting the convective part
    diffusion_matrix = *s.matrix.get_blocks_uniquely().at(0).get();
   }
}

/* ************************************************************************* */
template<int d> 
void ConvectionDiffusion_AFC<d>::solve_poisson()
{
  auto s = this->ConvectionDiffusion<d>::systems.front();
  FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
  //Setting the last diagonal entry of row to very large value so as
  //to apply mean value codition
  int nDofs = one_block.get_n_rows();
  const int * ColInd = one_block.get_vector_columns();
  const int * RowPtr = one_block.get_row_ptr();
  double * Entries = one_block.GetEntries();
  //Last row
  int i = nDofs - 1;
  //Scaling value
  double scale_value = 1e30;
  int found_entry = 0;

  for(int j = RowPtr[i]; j<RowPtr[i+1]; j++)
  {
    if(i == ColInd[j])
    {
      found_entry = j;
      Entries[j] *= scale_value;
      break;
    }
  }
  this->solver.solve(s.matrix, s.rhs, s.solution);
  //Resetting to original value (this step is not requuired, but better
  //safe than sorry)
  Entries[found_entry] /= scale_value;
  poisson_sol = s.solution;
}


#ifdef __3D__
template class ConvectionDiffusion_AFC<3>;
#else
template class ConvectionDiffusion_AFC<2>;
#endif
