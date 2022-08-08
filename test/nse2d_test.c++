/**
 * @brief A test program for the solving of NSE2D problems.
 *
 * This serves as a test for the solving of NSE2D problems. It is intended to
 * perform NSE2D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @author Naveed, Ulrich, Clemens
 *
 */
#include <cmath>
#include <Domain.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include "NavierStokes.h"
#include <Example_NSE2D.h>
#include <Multigrid.h>
#include <Chrono.h>
#include <algorithm>
#include "LocalAssembling.h"
#include "Utilities.h"
#include "ParMooN.h"

#include <ParMooN_repository_info.h>

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/data/mesh/";

const double accuracy = 1e-6;

void compare(const NavierStokes<2>& nse2d, std::array<double, int(5)> errors)
{
  auto computed_errors = nse2d.get_errors();

  // check the L2-error of the velcoity
  if( !utilities::are_equal(computed_errors[0], errors[0], accuracy) )
  {
    ErrThrow("L2 norm of velocity: ", std::setprecision(14), computed_errors[0],
        "  ", errors[0], "  ", std::abs(computed_errors[0] - errors[0]));
  }
  // check the L2-error of the divergence of the velocity
  if( !utilities::are_equal(computed_errors[1], errors[1], accuracy) )
  {
    ErrThrow("L2 norm of divergence: ", std::setprecision(14),
        computed_errors[1], "  ", errors[1], "  ", std::abs(computed_errors[1]
          - errors[1]));
  }
  // check the H1-error of the velcoity
  if( !utilities::are_equal(computed_errors[2], errors[2], accuracy) )
  {
    ErrThrow("H1 norm of velocity: ", std::setprecision(14), computed_errors[2],
        "  ", errors[2], "  ", std::abs(computed_errors[2] - errors[2]));
  }
  // check the L2-error of the pressure
  if( !utilities::are_equal(computed_errors[3], errors[3], accuracy))
  {
    ErrThrow("L2 norm of pressure: ", std::setprecision(14), computed_errors[3],
        "  ", errors[3], "  ", std::abs(computed_errors[3] - errors[3]));
  }
  // check the H1-error of the pressure
  if(!utilities::are_equal(computed_errors[4], errors[4], accuracy) )
  {
    ErrThrow("H1 norm of pressure: ", std::setprecision(14), computed_errors[4],
        "  ", errors[4], "  ", std::abs(computed_errors[4] - errors[4]));
  }
}

void compute(TDomain &domain, ParameterDatabase& db,
             std::array<double, int(5)> errors)
{
  NavierStokes<2> nse2d(domain, db);
  nse2d.assemble_linear_terms();
  // check stopping criterion
  nse2d.stop_it(0);
  for(unsigned int k=1;; k++)
  {
    Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
                     nse2d.get_residuals());
    nse2d.solve();
    // checking the first nonlinear iteration    
    nse2d.assemble_nonlinear_term();;
    if(nse2d.stop_it(k))
      break;
  }
  nse2d.output();
  // compare now the errors
  compare(nse2d, errors);
}

void check(TDomain &domain, ParameterDatabase db,
           int velocity_order, int nstype, int laplace_type,
           std::string nonlinear_form,
           std::array<double, int(5)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", nstype=", nstype, ", laplace_type=", laplace_type, 
                ", and nonlinear_form=", nonlinear_form);
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Multigrid::default_multigrid_database());
  db["problem_type"] = 5;
  db["solver_type"] = "direct";
  db["iterative_solver_type"] = "fgmres";
  db["residual_tolerance"] = 1.e-12;
  db["preconditioner"] = "least_squares_commutator";
  
  db["nonlinloop_maxit"] = 50;
  db["nonlinloop_epsilon"] = 1e-10;

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  db["nse_nonlinear_form"] = nonlinear_form;
  
  Chrono timer;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d direct solver,                    velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  
  db["solver_type"] = "iterative";
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(lsc preconditioner),       velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  

  db["preconditioner"] = "multigrid";
  db["multigrid_n_levels"] = db["refinement_n_initial_steps"].get<size_t>();
  //choose smoother on fine grid according to element
  std::vector<int> disc_p = {12,13,14,15,22,23,24};
  if(std::find(disc_p.begin(), disc_p.end(), velocity_order) != disc_p.end())
    db["multigrid_smoother"] = "cell_vanka_store";
  else
    db["multigrid_smoother"] = "patch_vanka_store";

  db["multigrid_type"] = "standard";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_n_pre_smooth"] = 0;
  db["multigrid_n_post_smooth"] = 1;
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_vanka_damp_factor"] = 1.0;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(multigrid preconditioner), velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
}

template <int laplace_type>
void check_one_element(TDomain& domain, ParameterDatabase db,
                       int velocity_order, std::string nonlinear_form,
                       std::array<double, int(5)> errors)
{
  if(laplace_type == 0 && nonlinear_form != "rotational"
     && nonlinear_form != "emac")
  {
    // NSTYPE = 1
    check(domain, db, velocity_order, 1, laplace_type, nonlinear_form, errors);
    // NSTYPE = 2
    check(domain, db, velocity_order, 2, laplace_type, nonlinear_form, errors);
  }
  // NSTYPE = 3
  check(domain, db, velocity_order, 3, laplace_type, nonlinear_form, errors);
  // NSTYPE is 4
  check(domain, db, velocity_order, 4, laplace_type, nonlinear_form, errors);
}

// =======================================================================
// routines including boundary errors
// =======================================================================


void compare_including_boundary(const NavierStokes<2>& nse2d, std::array<double, int(8)> errors)
{
  auto computed_errors = nse2d.get_errors();

  // check the L2-error of the velcoity
  if( std::abs(computed_errors[0]-errors[0]) > accuracy )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the L2-error of the divergence of the velocity
  if( std::abs(computed_errors[1] - errors[1]) > accuracy )
  {
    ErrThrow("L2 norm of divergence: ", std::setprecision(14), computed_errors[1], "  ", errors[1]);
  }
  // check the H1-error of the velcoity
  if( std::abs(computed_errors[2] - errors[2]) > accuracy )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[2], "  ", errors[2]);
  }
  // check the L2-error of the pressure
  if( std::abs(computed_errors[3] - errors[3]) > accuracy)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
  // check the H1-error of the pressure
  if(std::abs(computed_errors[4] - errors[4]) > accuracy )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[4], "  ", errors[4]);
  }
  // check the L2-error of the pressure
  if( std::abs(computed_errors[6] - errors[5]) > accuracy)
  {
    ErrThrow("L2 norm of velocity on the Nitsche boundary: ", computed_errors[6], "  ", errors[5]);
  }
  // check the H1-error of the pressure
  if(std::abs(computed_errors[7] - errors[6]) > accuracy )
  {
    ErrThrow("L2 norm of normal velocity on the Nitsche boundary: ", computed_errors[7], "  ", errors[6]);
  }
}

void compute_including_boundary(TDomain &domain, ParameterDatabase& db,
             std::array<double, int(8)> errors)
{
  NavierStokes<2> nse2d(domain, db);
  nse2d.assemble_linear_terms();
  // check stopping criterion
  nse2d.stop_it(0);

  nse2d.solve();

  nse2d.output();
  // compare now the errors
  compare_including_boundary(nse2d, errors);
}

void check_including_boundary(TDomain &domain, ParameterDatabase db,
           int velocity_order, int pressure_order, int nstype, int laplace_type,
           std::array<double, int(8)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", nstype=", nstype, ", laplace_type=", laplace_type);

  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Multigrid::default_multigrid_database());
  db["solver_type"] = "direct";
  //db["iterative_solver_type"] = "fgmres";
  //db["residual_tolerance"] = 1.e-12;
  //db["preconditioner"] = "least_squares_commutator";

  //db["nonlinloop_maxit"] = 50;
  //db["nonlinloop_epsilon"] = 1e-10;

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;

  Chrono timer;
  compute_including_boundary(domain, db, errors);
  timer.restart_and_print("nse2d-brinkman2d direct solver,  velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
}


template <int laplace_type>
void check_one_element_including_boundary(TDomain& domain, ParameterDatabase db,
                       int velocity_order, int pressure_order,
                       std::array<double, int(8)> errors)
{
  check_including_boundary(domain, db, velocity_order, pressure_order, 14, laplace_type, errors);
}




// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  parmoon::parmoon_initialize();
  bool testall = false;
  if (argv[1])
    testall = (std::string(argv[1]).compare("testall") == 0);

  Chrono chrono;

  /** Program 1
   *  This program tests direct solve with galerkin discretization
   * direct solver; Tests for Triangles
   * The VELOCITY_SPACE is fixed and tests for different NSTYPE's
   */
  {
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());

    db["problem_type"].set<size_t>(5);
    db["example"] = 2;

    db.add("refinement_n_initial_steps", (size_t) 2,"");

    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});

    // default construct a domain object
    TDomain domain(db);

    db["reynolds_number"] = 1;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE=5;
    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 0;

    // possibly parameters in the database
    check_parameters_consistency_NSE(db);
    // refine grid
    domain.refine_and_get_hierarchy_of_collections(db);
    std::array<double, int(5)> errors;

    //=========================================================================
    Output::print<1>("\nTesting the P2/P1 elements");
    errors = {{ 0.005456127622249, 0.085600792927298, 0.15659378093701,
            0.07108665290335, 1.3407010450111 }};
    // VELOCITY_SPACE = 2 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 2, "convective", errors);
    errors = {{ 0.005461517514555, 0.08567089916679, 0.156613079963,
            0.071844737124695, 1.3491718656292 }};
    check_one_element<0>(domain, db, 2, "skew_symmetric", errors);
    errors = {{ 0.006214685051306, 0.10151526663975, 0.16675196021083,
            0.11655892322349, 1.6480675404577 }};
    check_one_element<0>(domain, db, 2, "rotational", errors);
    errors = {{ 0.0057827846340661, 0.094178019255837, 0.16300544453388,
            0.1607781642646, 1.664867097309 }};
    check_one_element<0>(domain, db, 2, "emac", errors);
    errors = {{ 0.005461517514555, 0.08567089916679, 0.156613079963,
            0.071844737124695, 1.3491718656292 }};
    check_one_element<0>(domain, db, 2, "divergence", errors);

    //=========================================================================
    // This has noncontinuous pressure compared to the P2/P1 element
    Output::print<1>("\n Testing the lowest order nonconforming element");
    // VELOCITY_SPACE = -1 and the pressure space is chosen in the class NSE2D
    errors = {{0.068963470471598, 1.2485739511643e-15, 1.1714517336566,
      0.47722251294817, 2.2214414690792}};
    check_one_element<0>(domain, db, -1, "convective", errors);
    errors = {{0.070297368125753, 1.1266219631669e-15, 1.1746918577782,
      0.47707103735276, 2.2214414690792}};
    check_one_element<0>(domain, db, -1, "skew_symmetric", errors);
    errors = {{0.069043251632844, 8.2626436770319e-16, 1.1792880398777,
      0.49958277924691, 2.2214414690792}};
    check_one_element<0>(domain, db, -1, "rotational", errors);
    errors = {{0.074862605519878, 1.4994782176059e-15, 1.2074394025014,
      0.56239780077343, 2.2214414690792}};
    check_one_element<0>(domain, db, -1, "emac", errors);
    errors = {{0.068963470471598, 1.4086100394526e-15, 1.1714517336566,
      0.47722251294817, 2.2214414690792}};
    check_one_element<0>(domain, db, -1, "divergence", errors);


    if(testall)
    {
      //=======================================================================
      Output::print<1>("\nTesting the P3/P2 elements");
      errors = {{ 0.00028020829779641, 0.0061977013305506, 0.011391210186952,
              0.0056918081803184, 0.20769679500378 }};
      // VELOCITY_SPACE = 3 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 3, "convective", errors);
      errors = {{ 0.0002801683506044, 0.0061981261043577, 0.011391678805669,
              0.0056956176000077, 0.20780500907323 }};
      check_one_element<0>(domain, db, 3, "skew_symmetric", errors);
      errors = {{ 0.00034351517928429, 0.0096256533312287, 0.013780390342581,
              0.01555923134715, 0.39131977067853 }};
      check_one_element<0>(domain, db, 3, "rotational", errors);
      errors = {{ 0.00031706470111175, 0.0086454126016252, 0.012993003704192,
              0.014463947314389, 0.32790669662401 }};
      check_one_element<0>(domain, db, 3, "emac", errors);
      errors = {{ 0.00028020829779641, 0.0061977013305506, 0.011391210186952,
              0.0056956176000031, 0.20780500907308 }};
      check_one_element<0>(domain, db, 3, "divergence", errors);

      //=======================================================================
      Output::print<1>("\nTesting the P4/P3 elements");
      errors = {{ 1.1817023010728e-05, 0.00029575269112296, 0.0006435418450572,
              0.00050496270735108, 0.026998702772064 }};
      // VELOCITY_SPACE = 4 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 4, "convective", errors);

      //=======================================================================
      Output::print<1>("\nTesting the P5/P4 elements");
      errors = {{ 4.3466391168252e-07, 1.4197002168949e-05, 2.793323812439e-05,
              2.1824211773585e-05, 0.0016936362911126 }};
      // VELOCITY_SPACE = 5 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 5, "convective", errors);

      //=========================================================================
      Output::print<1>("\n Testing the MINI element");
      // VELOCITY_SPACE = 101 and the pressure space is chosen in the class NSE2D
      errors = {{0.08373116847134, 0.93872389099579, 1.4254617546805,
        0.92498130219029, 9.0076641032004}};
      check_one_element<0>(domain, db, 101, "convective", errors);
      errors = {{0.08374226497677, 0.93884772885449, 1.4257610818885,
        0.92469057864223, 9.0273708960303}};
      check_one_element<0>(domain, db, 101, "skew_symmetric", errors);
      errors = {{0.083374515807732, 0.94105838050009, 1.42644014532,
        0.87417152971762, 8.9316101204153}};
      check_one_element<0>(domain, db, 101, "rotational", errors);
      errors = {{0.084141011220491, 0.9378123002968, 1.4265321249254,
        0.98285538577688, 9.1582546158151}};
      check_one_element<0>(domain, db, 101, "emac", errors);
      errors = {{0.08374226497677, 0.93884772885449, 1.4257610818885,
        0.92469057864224, 9.0273708960303}};
      check_one_element<0>(domain, db, 101, "divergence", errors);

      //=========================================================================
      Output::print<1>("\nTesting the P2-bubble/P1-disc elements");
      errors = {{ 0.0071886299046824, 0.081842610256743, 0.21185558654462,
        0.36754876295023, 5.3058522557418 }};
      // VELOCITY_SPACE = 22 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 22, "convective", errors);
      errors = {{ 0.0071870075217203, 0.081812991627381, 0.21185522222301,
        0.36760565517932, 5.3064844653778 }};
      check_one_element<0>(domain, db, 22, "skew_symmetric", errors);
      errors = {{ 0.0071856730855546, 0.081567156856757, 0.2120767352661,
        0.37174134839162, 5.4006410196856 }};
      check_one_element<0>(domain, db, 22, "rotational", errors);
      errors = {{ 0.0071979487572762, 0.082123462115766, 0.21189447518786,
        0.37027289021106, 5.3911935934541 }};
      check_one_element<0>(domain, db, 22, "emac", errors);
      errors = {{ 0.0071870075217202, 0.081812991627381, 0.21185522222301,
        0.36760565517932, 5.3064844653777 }};
      check_one_element<0>(domain, db, 22, "divergence", errors);

      //=========================================================================
      Output::print<1>("\nTesting the P3-bubble/P2-disc elements");
      errors = {{ 0.00026037876329682, 0.005726948664285, 0.010856952083041,
              0.013370326976055, 0.39882948355395 }};
      // VELOCITY_SPACE = 23 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 23, "convective", errors);

      //=========================================================================
      Output::print<1>("\nTesting the P4-bubble/P3-disc elements");
      errors = {{ 1.2032722339771e-05, 0.00024376809624887, 0.00055164963287203,
              0.00063706731983293, 0.027783948983068 }};
      // VELOCITY_SPACE = 24 and the pressure space is chosen in the class NSE2D
      check_one_element<0>(domain, db, 24, "convective", errors);
    }
  }   // end program 1
  chrono.restart_and_print("test part 1");
  //=========================================================================
  /** Program 2
   *  This program tests direct solve with galerkin discretization
   * direct solver; Test for Quad's
   * The VELOCITY_SPACE is fixed and tests for different NSTYPE's
   */
  {
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db["problem_type"].set<size_t>(5);
    db["example"] = 2;

    db.add("refinement_n_initial_steps", (size_t) 2,"");

    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    // default construct a domain object
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);

    // parameters used for this test
    db["reynolds_number"] = 1;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE=5;
    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->LAPLACETYPE = 0;

    // possibly parameters in the database
    check_parameters_consistency_NSE(db);

    // refine grid
    domain.refine_and_get_hierarchy_of_collections(db);
    std::array<double, int(5)> errors;

    //=========================================================================
    Output::print<1>("\nTesting the Q2/Q1 elements");
    errors = {{ 0.004083204524442, 0.05101317284072, 0.10522635824261,
            0.017686667902813, 0.51182308944019 }};
    // VELOCITY_SPACE  = 2
    check_one_element<0>(domain, db, 2, "convective", errors);
    errors = {{ 0.0040835090784135, 0.051016841655666, 0.10523293839479,
            0.017611281280989, 0.51158180834432 }};
    check_one_element<0>(domain, db, 2, "skew_symmetric", errors);
    errors = {{ 0.0046884833250869, 0.070877525721747, 0.11742828649332,
            0.11040259260052, 1.0352007333739 }};
    check_one_element<0>(domain, db, 2, "rotational", errors);
    errors = {{ 0.0045805109059847, 0.066717957084682, 0.11493220231642,
            0.11124041590661, 0.94502906380428 }};
    check_one_element<0>(domain, db, 2, "emac", errors);
    errors = {{ 0.0040835090784135, 0.051016841655666, 0.10523293839479,
            0.017611281280989, 0.51158180834432 }};
    check_one_element<0>(domain, db, 2, "divergence", errors);

    if(testall)
    {
      //=========================================================================
      Output::print<1>("\nTesting the Q3/Q2 elements");
      errors = {{ 0.00019319433716041, 0.003563918945437, 0.0071078507849009,
              0.0018446328461379, 0.057123632497266 }};
      // VELOCITY_SPACE  = 3
      check_one_element<0>(domain, db, 3, "convective", errors);
      errors = {{ 0.00019316638866158, 0.0035638246364146, 0.0071078664938289,
              0.0018454762090097, 0.057115940978707 }};
      check_one_element<0>(domain, db, 3, "skew_symmetric", errors);
      errors = {{ 0.00026177332079517, 0.0079927161285491, 0.010243629124348,
              0.0094414547091167, 0.21505945593453 }};
      check_one_element<0>(domain, db, 3, "rotational", errors);
      errors = {{ 0.00024830802922935, 0.0072493956722643, 0.0096501647507086,
              0.0094088997257092, 0.18990292685865 }};
      check_one_element<0>(domain, db, 3, "emac", errors);

      //=========================================================================
      Output::print<1>("\nTesting the Q4/Q3 elements");
      errors = {{ 6.9434747253041e-06, 0.0001677182257252, 0.00035212311261646,
              9.4703269177756e-05, 0.0048368160352994 }};
      // VELOCITY_SPACE  = 4
      check_one_element<0>(domain, db, 4, "convective", errors);

      //=========================================================================
      Output::print<1>("\nTesting the Q5/Q4 elements");
      errors = {{ 2.3951237974726e-07, 6.8246180893733e-06, 1.4163394749406e-05,
              4.9000676557526e-06, 0.00030469183949993 }};
      // VELOCITY_SPACE  = 5
      check_one_element<0>(domain, db, 5, "convective", errors);
    }

    //=========================================================================
    Output::print<1>("\nTesting the Q2/P1-disc elements");
    errors = {{ 0.0040557267369371, 0.05094007770388, 0.10564123627325,
            0.030975452768144, 0.70842776239597 }};
    // VELOCITY_SPACE  = 22
    check_one_element<0>(domain, db, 22, "convective", errors);
    errors = {{ 0.0040565847234025, 0.050939919166369, 0.10564293309485,
            0.030954083022779, 0.70858252999859 }};
    check_one_element<0>(domain, db, 22, "skew_symmetric", errors);
    errors = {{ 0.0041643538522174, 0.055554927078881, 0.10980572753403,
            0.036275829042111, 0.75595849613121 }};
    check_one_element<0>(domain, db, 22, "rotational", errors);
    errors = {{ 0.0041393709334753, 0.053107233911226, 0.10760054343883,
            0.033157637926479, 0.72384438559736 }};
    check_one_element<0>(domain, db, 22, "emac", errors);
    errors = {{ 0.0040557267369371, 0.05094007770388, 0.10564293309485,
            0.030954083022779, 0.70858252999859 }};
    check_one_element<0>(domain, db, 22, "divergence", errors);
    if(testall)
    {
      //=========================================================================
      Output::print<1>("\nTesting the Q3/P2-disc elements");
      errors = {{ 0.00020642736694367, 0.0042373143298489, 0.0075794259144329,
              0.0039793392063018, 0.13655739379564 }};
      // VELOCITY_SPACE  = 23
      check_one_element<0>(domain, db, 23, "convective", errors);

      //=========================================================================
      Output::print<1>("\nTesting the Q4/P3-disc elements");
      errors = {{ 9.9544541389566e-06, 0.00032338687271106, 0.00045888380841742,
              0.00037625559674667, 0.01835259073302 }};
      // VELOCITY_SPACE  = 24
      check_one_element<0>(domain, db, 24, "convective", errors);

      //=========================================================================
      Output::print<1>("\nTesting the Q5/P4-disc elements");
      errors = {{ 5.3747286967048e-07, 2.4122655352531e-05, 2.7997005479668e-05,
              2.81141239001e-05, 0.0018176331985532 }};
      // VELOCITY_SPACE  = 25
      check_one_element<0>(domain, db, 25, "convective", errors);
    }
  }
  chrono.restart_and_print("test part 2");

  //=========================================================================
  /** Program 3
   *  This program tests direct solve with nonsymm_gls discretization and the Nitsche method
   * direct solver; Tests for Triangles
   */
  {
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(NavierStokes<2>::default_nse_database());
    //ParameterDatabase db = NavierStokes<2>::default_nse_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());

    std::array<double, int(8)> all_errors; // includes errors at the boundary (Nitsche)

    db["problem_type"].set<size_t>(7);

    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 0,"", (size_t) 0, (size_t) 20);

    db["output_compute_errors"] = true;

    TDatabase::ParamDB->NSTYPE = 14;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = true;

    // possibly parameters in the database
    check_parameters_consistency_NSE(db);

    int velocity_space, pressure_space;
    std::vector<size_t> neumann_id, nitsche_id;
    std::vector<double> neumann_value, nitsche_penalty;

    // ================================ EXAMPLE: Poiseuille_Flow =======================================

    db["boundary_file"].set("Default_UnitSquare");
    db["geo_file"].set("TwoTriangles");

    db["example"] = 6;
    db["refinement_n_initial_steps"].set<size_t>(3);

    velocity_space = 2;
    pressure_space = 1;

    // coefficients
    db["effective_viscosity"] = 0.00025;
    db["inverse_permeability"] = 0.004;

    // stabilization
    db["space_discretization_type"] = "galerkin";
    db["gls_stab"] = 0.;
    db["graddiv_stab"] = 0.;
    db["corner_stab"] = 0.;
    db["L_0"] = 0.;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 2;
    neumann_id = {1, 3};
    db["neumann_id"] = neumann_id;
    neumann_value = {-0.5, 0.5};
    db["neumann_value"] = neumann_value;

    db["n_nitsche_bd"] = 0;

    // default construct a domain object
    TDomain domain0(db);
    // refine grid
    domain0.refine_and_get_hierarchy_of_collections(db);

    Output::print<1>("\n Example:  Poiseuille_Flow: Testing P2/P1 (with gls_stab) for Brinkman2D via NSE2D");
    all_errors = {{ 0.074230400439798,  0.29287548273616, 4.0915592126458, 6.111062446254e-05,
            0.0013491125315994, 0., 0. }};

    check_one_element_including_boundary<0>(domain0, db, velocity_space, pressure_space, all_errors);

    // ================================ EXAMPLE: SinCos_DarcyFlow =======================================

    db["boundary_file"].set("Default_UnitSquare");
    db["geo_file"].set("TwoTriangles");

    db["example"] = 7;
    db["refinement_n_initial_steps"].set<size_t>(5);

    velocity_space = 1;
    pressure_space = 1;

    // coeficients
    db["effective_viscosity"] = 0.;
    db["inverse_permeability"] = 0.1;

    // stabilization
    db["space_discretization_type"] = "nonsymm_gls";
    db["gls_stab"] = 0.1;
    db["graddiv_stab"] = 0.1;
    db["corner_stab"] = 1.;
    db["L_0"] = 0.1;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 0;

    db["n_nitsche_bd"] = 4;
    nitsche_id = {0, 1, 2, 3};
    db["nitsche_id"] = nitsche_id;
    nitsche_penalty = {0., 0., 0., 0.};
    db["nitsche_penalty"] = nitsche_penalty;

    db["symmetric_nitsche_u"] = -1;
    db["symmetric_nitsche_p"] = -1;

    // default construct a domain object
    TDomain domain1(db);
    // refine grid
    domain1.refine_and_get_hierarchy_of_collections(db);

    Output::print<1>("\n Example: SinCos_DarcyFlow: Testing P1/P1 (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
    all_errors = {{ 0.014963869771925, 2.7518157106309, 3.9062600370977, 0.0001701830587647,
            0.044106028987814, 0.047886464620468, 0.039273676677846 }};

    check_one_element_including_boundary<0>(domain1, db, velocity_space, pressure_space, all_errors);

    // ================================ EXAMPLE: Discacciati_Flow =======================================

    db["boundary_file"].set(path_to_repo + "Rectangle.PRM", false);
    db["geo_file"].set(path_to_repo + "discacciati_11K.mesh", false);

    db["example"] = 8;
    db["refinement_n_initial_steps"].set<size_t>(0);

    velocity_space = 1;
    pressure_space = 1;

    // coefficients
    db["effective_viscosity"] = 0.01;
    db["inverse_permeability"] = 0.;

    // stabilization
    db["space_discretization_type"] = "nonsymm_gls";
    db["gls_stab"] = 0.1;
    db["graddiv_stab"] = 0.1;
    db["corner_stab"] = 0.1;
    db["L_0"] = 0.1;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 1;
    neumann_id = {1};
    db["neumann_id"] = neumann_id;
    neumann_value = {0.};
    db["neumann_value"] = neumann_value;

    db["n_nitsche_bd"] = 3;
    nitsche_id = {0, 2, 3};
    db["nitsche_id"] = nitsche_id;
    nitsche_penalty = {0., 0., 0.};
    db["nitsche_penalty"] = nitsche_penalty;

    db["symmetric_nitsche_u"] = -1;
    db["symmetric_nitsche_p"] = -1;

    // default construct a domain object
    TDomain domain2(db);
    // refine grid
    domain2.refine_and_get_hierarchy_of_collections(db);

    Output::print<1>("\n Example: Discacciati_Flow: Testing P1/P1 (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
    all_errors = {{ 0.92147504957481, 0.047176605078917, 1.3065668883004, 0.23682311417141,
            0.044785943200984, 0.0067226743785503, 0.0003305737320877 }};

    check_one_element_including_boundary<0>(domain2, db, velocity_space, pressure_space, all_errors);

    // ================================ EXAMPLE: Radial_Flow_with_Hole =======================================

    db["boundary_file"].set(path_to_repo + "ring.PRM", false);
    db["geo_file"].set(path_to_repo + "ring.mesh", false);

    db["example"] = 9;
    db["refinement_n_initial_steps"].set<size_t>(0);

    velocity_space = 1;
    pressure_space = 1;

    // coefficients
    db["effective_viscosity"] = 0.;
    db["inverse_permeability"] = 100.;

    // stabilization
    db["space_discretization_type"] = "nonsymm_gls";
    db["gls_stab"] = 0.1;
    db["graddiv_stab"] = 0.1;
    db["corner_stab"] = 1.;
    db["L_0"] = 0.1;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 1;
    neumann_id = {0};
    db["neumann_id"] = neumann_id;
    neumann_value = {0.};
    db["neumann_value"] = neumann_value;

    db["n_nitsche_bd"] = 0;

    // default construct a domain object
    TDomain domain3(db);
    // refine grid
    domain3.refine_and_get_hierarchy_of_collections(db);



    Output::print<1>("\n Example: Radial_Flow_with_Hole: Testing P1/P1 (with gls_stab and Dirichlet) for Brinkman2D via NSE2D");
    all_errors = {{ 0.00043540827744988, 8.6535486158172e-06, 0.0052597278868421, 0.13078275084169,
            0.043555740117776, 0, 0 }};

    check_one_element_including_boundary<0>(domain3, db, velocity_space, pressure_space, all_errors);


    db["n_nitsche_bd"] = 1;
    nitsche_id = {1};
    db["nitsche_id"] = nitsche_id;
    nitsche_penalty = {0.};
    db["nitsche_penalty"] = nitsche_penalty;

    db["symmetric_nitsche_u"] = -1;
    db["symmetric_nitsche_p"] = -1;

    Output::print<1>("\n Example: Radial_Flow_with_Hole: Testing P1/P1 (with gls_stab and Nitsche method) for Brinkman2D via NSE2D on a coarse grid");
    all_errors = {{  0.00043531021603966, 8.3872506583061e-06, 0.0052603289724414, 0.1307805526847,
            0.043526205098935, 4.3672111775544e-07, 4.3297741843194e-07 }};



    check_one_element_including_boundary<0>(domain3, db, velocity_space, pressure_space, all_errors);


    // ================================ EXAMPLE: circle_with_immersed_hole =======================================

    db["boundary_file"].set(path_to_repo + "disk.PRM", false);
    db["geo_file"].set(path_to_repo + "disk.mesh", false);

    db["example"] = 10;
    db["refinement_n_initial_steps"].set<size_t>(3);

    velocity_space = 1;
    pressure_space = 1;

    // coefficients
    db["effective_viscosity"] = 0.01;
    db["inverse_permeability"] = 100.;

    // stabilization
    db["space_discretization_type"] = "nonsymm_gls";
    db["gls_stab"] = 0.1;
    db["graddiv_stab"] = 0.1;
    db["corner_stab"] = 0.;
    db["L_0"] = 0.1;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 1;
    neumann_id = {0};
    db["neumann_id"] = neumann_id;
    neumann_value = {0.};
    db["neumann_value"] = neumann_value;

    db["n_nitsche_bd"] = 0;


    // default construct a domain object
    TDomain domain4(db);
    // refine grid
    domain4.refine_and_get_hierarchy_of_collections(db);

    Output::print<1>("\n Example: circle_with_immersed_hole: Testing P1/P1 (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
    all_errors = {{1.574236335318e-07, 4.4521747107522e-07, 3.1481651364273e-07,
                   6.4284369400791e-06, 1.574167346553e-05, 0, 0 }};

    check_one_element_including_boundary<0>(domain4, db, velocity_space, pressure_space, all_errors);

    // ================================ EXAMPLE: Exponential_Poiseuille_Flow =======================================

    db["boundary_file"].set("Default_UnitSquare");
    db["geo_file"].set("TwoTriangles");

    db["example"] = 12;
    db["refinement_n_initial_steps"].set<size_t>(5);

    velocity_space = 1;
    pressure_space = 1;

    // coefficients
    db["effective_viscosity"] = 0.1;
    db["inverse_permeability"] = 10.;

    // stabilization
    db["space_discretization_type"] = "nonsymm_gls";
    db["gls_stab"] = 0.1;
    db["graddiv_stab"] = 0.1;
    db["corner_stab"] = 1.;
    db["L_0"] = 0.1;
    db["pspg_delta0"] = 0.;

    // boundary conditions
    db["n_neumann_bd"] = 2;
    neumann_id = {1, 3};
    db["neumann_id"] = neumann_id;
    neumann_value = {-0.5, 0.5};
    db["neumann_value"] = neumann_value;

    db["n_nitsche_bd"] = 2;
    nitsche_id = {0, 2};
    db["nitsche_id"] = nitsche_id;
    nitsche_penalty = {0., 0.};
    db["nitsche_penalty"] = nitsche_penalty;

    db["symmetric_nitsche_u"] = -1;
    db["symmetric_nitsche_p"] = -1;

    // default construct a domain object
    TDomain domain6(db);
    // refine grid
    domain6.refine_and_get_hierarchy_of_collections(db);

    Output::print<1>("\n Example: Exponential_Poiseuille_Flow: Testing P1/P1 (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
    all_errors = {{  0.002262597973565, 0.0033740583845405,  0.038585873042087,
            0.00051923122268645, 0.034194588222972, 0.0059477846442968, 0.00011506854321844 }};

    check_one_element_including_boundary<0>(domain6, db, velocity_space, pressure_space, all_errors);


    // ================================ EXAMPLE: T-Shaped Cavity Flow =======================================

       db["boundary_file"].set(path_to_repo + "Tshape.PRM", false);
       db["geo_file"].set(path_to_repo + "Tshape_meshCodina_coarse.mesh", false);

       db["example"] = 16;
       db["refinement_n_initial_steps"].set<size_t>(1);

       velocity_space = 1;
       pressure_space = 1;

       // coefficients
       db["effective_viscosity"] = 1.;
       db["inverse_permeability"] = 10000.;

       // stabilization
       db["space_discretization_type"] = "nonsymm_gls";
       db["gls_stab"] = 0.1;
       db["graddiv_stab"] = 0.1;
       db["corner_stab"] = 0.1;
       db["L_0"] = 0.1;
       db["pspg_delta0"] = 0.;

       // boundary conditions
       db["n_neumann_bd"] = 1;
       neumann_id = {3};
       db["neumann_id"] = neumann_id;
       neumann_value = {0.};
       db["neumann_value"] = neumann_value;

       db["n_nitsche_bd"] = 6;
       nitsche_id = {0, 1, 2, 4, 6, 7};
       db["nitsche_id"] = nitsche_id;
       nitsche_penalty = {0., 0., 0., 0., 0., 0., 0.};
       db["nitsche_penalty"] = nitsche_penalty;

       db["symmetric_nitsche_u"] = -1;
       db["symmetric_nitsche_p"] = -1;

       // default construct a domain object
       TDomain domain7(db);
       // refine grid
       domain7.refine_and_get_hierarchy_of_collections(db);

       Output::print<1>("\n Example: Tshaped_Cavity_Flow: Testing P1/P1 (a) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
       all_errors = {{  3.6831198051424, 0.89212355828593,  2.2560497495104,
               349387.30630147, 36899.794706961, 2.6649800657205, 0.50524785606482 }};

       check_one_element_including_boundary<0>(domain7, db, velocity_space, pressure_space, all_errors);

       // coefficients
       db["effective_viscosity"] = 1.;
       db["inverse_permeability"] = 0.;

       Output::print<1>("\n Example: Tshaped_Cavity_Flow: Testing P1/P1 (b) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
       all_errors = {{  4.5070979045055, 0.088356859410387,  2.1559012409417,
               6.3688210021327, 1.4325118593401, 0.068212452525502, 0.039259574501929 }};

       check_one_element_including_boundary<0>(domain7, db, velocity_space, pressure_space, all_errors);

       // coefficients
       db["effective_viscosity"] = 0.;
       db["inverse_permeability"] = 1.;

       Output::print<1>("\n Example: Tshaped_Cavity_Flow: Testing P1/P1 (c) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
       all_errors = {{  3.6815404248534, 0.89632877902632,  2.2601798212013,
               34.913661701687, 3.6883042211571, 2.6737506794807, 0.50776168487161 }};

       check_one_element_including_boundary<0>(domain7, db, velocity_space, pressure_space, all_errors);


     /*  // ================================ EXAMPLE: Riverbed =======================================

             db["boundary_file"].set(path_to_repo + "doubleRiverbed.PRM", false);
             db["geo_file"].set(path_to_repo + "doubleRiverbed3.mesh", false);

             db["example"] = 15;
             db["refinement_n_initial_steps"].set<size_t>(0);

             velocity_space = 1;
             pressure_space = 1;

             // coefficients
             db["effective_viscosity"] = -2;
             db["inverse_permeability"] = -2;

             // stabilization
             db["space_discretization_type"] = "nonsymm_gls";
             db["gls_stab"] = 0.1;
             db["graddiv_stab"] = 0.1;
             db["corner_stab"] = 0.;
             db["pspg_delta0"] = 0.;

             // boundary conditions
             db["n_neumann_bd"] = 4;
             neumann_id = {0,2,3,5};
             db["neumann_id"] = neumann_id;
             neumann_value = {0.001, 0., 0., 0.001};
             db["neumann_value"] = neumann_value;

             db["n_nitsche_bd"] = 2;
             nitsche_id = {1, 4};
             db["nitsche_id"] = nitsche_id;
             nitsche_penalty = {0., 0.};
             db["nitsche_penalty"] = nitsche_penalty;

             db["symmetric_nitsche_u"] = -1;
             db["symmetric_nitsche_p"] = -1;

             db["L_0"] = 0.05;
             // default construct a domain object
             TDomain domain8(db);
             // refine grid
             domain8.refine_and_get_hierarchy_of_collections(db);

             Output::print<1>("\n Example: Riverbed with penalty-free Nitsche on top and bottom boundaries:"
                     " Testing P1/P1 (level 0) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
             all_errors = {{ 0.0098949002531999, 0.0078617337273354, 0.06541610153892,
                     0.0011773430743204,  0.0011553255168328, 0.0032349028425138, 0.00016157198010113 }};

             check_one_element_including_boundary<0>(domain8, db, velocity_space, pressure_space, all_errors);


             db["refinement_n_initial_steps"].set<size_t>(1);
             db["L_0"] = 0.088;

             domain8.refine_and_get_hierarchy_of_collections(db);

             Output::print<1>("\n Example: Riverbed: Testing P1/P1 (level 1) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
             all_errors = {{ 0.0087021237436807, 0.0038554325939445, 0.061231312739041,
                    0.0011788015851962, 0.0011926543200593, 0.00082239121634021, 2.7745707187569e-05}};

             check_one_element_including_boundary<0>(domain8, db, velocity_space, pressure_space, all_errors);

             db["L_0"] = 0.044;
             db["refinement_n_initial_steps"].set<size_t>(2);
             domain8.refine_and_get_hierarchy_of_collections(db);

             Output::print<1>("\n Example: Riverbed: Testing P1/P1 (level 2) (with gls_stab and Nitsche method) for Brinkman2D via NSE2D");
             all_errors = {{0.0083914212116807, 0.0021100202716435, 0.060053746416856, 0.0011797537118376,
                 0.0017617376138521, 0.00021049602899397, 5.9828796100087e-06 }};

             check_one_element_including_boundary<0>(domain8, db, velocity_space, pressure_space, all_errors);
*/
  } // end program 3
  chrono.stop_and_print("test part 3");
  chrono.print_total_time("entire test nse2d");
  parmoon::parmoon_finalize();
}
