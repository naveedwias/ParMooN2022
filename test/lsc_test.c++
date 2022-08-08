/**
 * A test for our 2D implementation of the least squares commutator preconditioner.
 * The checks, whether the same number gmres iterations as reported by the
 * authors of the LSC in:
 *
 * Elman, Silvester & Wathen (2005): Finite Elements and Fast Iterative Solvers.
 * 1 ed. Oxford Science Publications. p.357, Table 8.2.
 *
 * is acheived for the lid driven cavity problem in a Q2/Q1 discretization,
 * with different Reynolds numbers and on different spatial discretization levels
 * (Picard iteration only).
 *
 * In a preparation step, the relative residual of the nonlinear problem
 * (relative to the norm of the initial right hand side) is reduced below 10^-5.
 * In the measuring step, the GMRES iteration must reduce the residual of the
 * linear problem below "10^-6 times the residual of the nonlinear problem".
 * The number of GMRES steps that are performed in the measuring step is the
 * target quantity of this test, because these are reported in the book.
 *
 * Clemens Bartsch, 2017/05/15
 */
#include <Domain.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include <LoopInfo.h>
#include "NavierStokes.h"
#include <Solver.h>
#include <Example_NSE2D.h>
#include "LocalAssembling.h"
#include "ParMooN.h"

#include <algorithm>

std::vector<int> reynolds_numbers = {10, 100};//,1000};

// Note: our ref level "3" corresponds to refinement level "4" in the book and
// so on. Level 6 is not included in the for the sake of computing time (but the
// target number of gmres iterations is reproduced correctly)
std::vector<int> n_ref_steps = {3, 4, 5};//, 6};

// Note: The reported gmres iteration numbers for RE=1000 are not hit exactly.
// The reason is probably the solution of the nonlinear system, i.e., the start
// solution of the "measuring step".
std::vector< std::vector<int> > target_n_gmres_steps =
    {
      {8, 11, 14},//, 18}, // Reynolds number   10
      {13, 16, 21},// 27}, // Reynolds number  100
      //{55, 62, 55, 45}   // Reynolds number 1000
    };

// Those are the numbers I found with the boundary corrected lsc. They are not
// yet checked against the 2nd edition of the book of Elman, but they show the
// expected behaviour: the number of iterations is independent of 'h'.
//std::vector< std::vector<int> > target_n_gmres_steps_bc_lsc =
//    {
//      {10, 10, 9, 10}, //Reynolds number 10
//      {14, 14, 14, 13} //Reynolds number 100
//    };


// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
    auto db = parmoon::parmoon_initialize();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(Solver<>::default_solver_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());

    db["problem_type"].set<size_t>(5);
    db["example"] = 1; //lid driven cavity example

    db.add("boundary_file", "DrivenCavitySquare", "");
    db.add("geo_file", "DrivenCavitySquareQuads", "");
    db.add("refinement_n_initial_steps", (size_t) 1, " ",(size_t) 1, (size_t)10);

    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->VELOCITY_SPACE = 2;
    TDatabase::ParamDB->PRESSURE_SPACE = 1;
    TDatabase::ParamDB->NSTYPE = 1;
    TDatabase::ParamDB->LAPLACETYPE = 0;

    for(size_t re_nr_index = 0; re_nr_index < reynolds_numbers.size(); ++re_nr_index)
    {
      //divided by 2, because: RE = 2/nu in paper, but RE=1/nu in ParMooN
      db["reynolds_number"] = 0.5*reynolds_numbers[re_nr_index];

      Output::print();
      Output::print("****** Checks with RE=", reynolds_numbers.at(re_nr_index), " ******" );

      for(size_t ref_st_index = 0; ref_st_index < n_ref_steps.size(); ++ref_st_index)
      {
        db["refinement_n_initial_steps"] = n_ref_steps[ref_st_index];
        Output::print("Refinement level ", n_ref_steps[ref_st_index]);

        // construct a domain object
        TDomain domain(db);

        // refine grid
        domain.refine_and_get_hierarchy_of_collections(db);

        Example_NSE2D example(db);

        //solver for nonlinear loop
        db["solver_type"]= "direct";
        db["nonlinloop_minit"] = 0;
        db["nonlinloop_maxit"] = 15;
        db["nonlinloop_slowfactor"] = 10000;
        db["nonlinloop_epsilon"] = 1.0e-5;
        db["nonlinloop_scale_epsilon_with_size"] = false;
        db["nonlinloop_residual_relative_to_rhs"] = true;

        // create an object of the Navier-Stokes class
        NavierStokes<2> ns(domain, db, example);

        ns.assemble_linear_terms();
        ns.stop_it(0);

        LoopInfo<Residuals> loop_info("nonlinear");
        loop_info.print_time_every_step = true;
        loop_info.verbosity_threshold = 1; // full verbosity
        loop_info.print(0, ns.get_residuals());

        //======================================================================
        // iterate with a direct solver until relative residual of 10^-5 is hit
        for(unsigned int k = 1;; k++)
        {
          ns.solve();
          ns.assemble_nonlinear_term();
          if(ns.stop_it(k))
          {
            loop_info.finish(k, ns.get_residuals());
            break;
          }
          else
            loop_info.print(k, ns.get_residuals());
        }

        {
          Output::print("Start measuring step");
          db["solver_type"] = "iterative";
          db["iterative_solver_type"] = "right_gmres";
          db["preconditioner"] = "least_squares_commutator";
          db["max_n_iterations"] = 100; //should not be hit
          db["gmres_restart"] = 100;    //should not be hit
          db["residual_tolerance"] = 1.0e-6 * ns.get_full_residual();
          db["residual_reduction"] = 1.0e-20;  //should not be hit

          NavierStokes<2> ns_2(domain, db, example);

          //insert solution from the nonlinear loop into the fresh object
          ns_2.get_solution() = ns.get_solution();

          ns_2.assemble_linear_terms();
          ns_2.assemble_nonlinear_term();

          ns.stop_it(0);
          ns_2.solve();
          ns_2.stop_it(1);

          int n_its = ns_2.get_it_solver_info().get_n_previous_iterations();

          // Now check the test quantity:
          if(n_its != target_n_gmres_steps[re_nr_index][ref_st_index])
            ErrThrow("Incorrect number of GMRES iterations. ", n_its, " ",
                     target_n_gmres_steps[re_nr_index][ref_st_index]);
        }


      }
    }
  parmoon::parmoon_finalize();
}
