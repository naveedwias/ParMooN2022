/**
 * @brief A test program for the solving of CD2D problems with adaptively
 * refined grids using hanging nodes.
 *
 * 
 * @date 2019/05/14
 * @author Abhinav Jha
 */
#include <cmath>

#include "Domain.h"
#include "Database.h"
#include "ConvectionDiffusion.h"
#include "NavierStokes.h"
#include "RefinementStrategy.h"
#include "LoopInfo.h"
#include "Chrono.h"
#include "LocalAssembling.h"
#include "Multigrid.h"
#include "ParMooN.h"

bool write_PS_and_VTU_files = false; // turn this on for debugging purposes
std::string PS_name(size_t curr_level, bool do_conforming_closure,
                    bool simplices, int ansatz_order);
double accuracy = 1e-11;

#ifdef __2D__
typedef ConvectionDiffusion<2> ConvDiff;
typedef NavierStokes<2> Stokes;
#else // 3D
typedef ConvectionDiffusion<3> ConvDiff;
typedef NavierStokes<3> Stokes;
#endif



void compare(const ConvDiff& cd, const Stokes& stokes,
             std::array<double, 9> errors)
{
  // check the cd errors
  if( std::abs(cd.get_L2_error() - errors[0]) > accuracy )
  {
    ErrThrow("Convection-Diffusion: L2 error not correct. ",
             cd.get_L2_error() - errors[0]);
  }
  if( std::abs(cd.get_H1_semi_error() - errors[1]) > accuracy )
  {
    ErrThrow("Convection-Diffusion: H1-semi error not correct. ",
             cd.get_H1_semi_error() - errors[1]);
  }
  if( std::abs(cd.get_SD_error() - errors[2]) > accuracy )
  {
    ErrThrow("Convection-Diffusion: sd error not correct.",
             cd.get_SD_error() - errors[2]);
  }
  if( std::abs(cd.get_L_inf_error() - errors[3]) > accuracy )
  {
    ErrThrow("Convection-Diffusion: L_inf error not correct.",
             cd.get_L_inf_error() - errors[3]);
  }

  auto stokes_errors = stokes.get_errors();
  if( std::abs(stokes_errors[0] - errors[4]) > accuracy )
  {
    ErrThrow("Stokes: L2 velocity error not correct. ",
             stokes_errors[0] - errors[4]);
  }
  if( std::abs(stokes_errors[1] - errors[5]) > accuracy )
  {
    ErrThrow("Stokes: L2 divergence error not correct. ",
             stokes_errors[1] - errors[5]);
  }
  if( std::abs(stokes_errors[2] - errors[6]) > accuracy )
  {
    ErrThrow("Stokes: H1 velocity error not correct. ",
             stokes_errors[2] - errors[6]);
  }
  if( std::abs(stokes_errors[3] - errors[7]) > accuracy )
  {
    ErrThrow("Stokes: L2 pressure error not correct. ",
             stokes_errors[3] - errors[7]);
  }
  if( std::abs(stokes_errors[4] - errors[8]) > accuracy )
  {
    ErrThrow("Stokes: H1 pressure error not correct. ",
             stokes_errors[4] - errors[8]);
  }
}


// indicate where to refine when using TDomain::RefineByIndicator.
// Cells which intersect with the zero set of this indicator function are 
// regularly refined.
void indicator(double x, double y, double* values)
{
  values[0] = 0.5*x*x + y - 0.5;
}

void indicator(double x, double y, double z, double* values)
{
  values[0] = 0.5*x*x + y + z - 0.5;
}


void test_indicator(bool do_conforming_closure, bool simplices, 
                    int ansatz_order, std::array<double, 9> errors)
{
  Output::print("\n", std::setfill('#'), std::setw(80), "#",
                "\n\nnext test:", std::boolalpha,
                "\n  conforming_closure: ", do_conforming_closure,
                "\n  simplices: ", simplices,
                "\n  ansatz_order: ", ansatz_order, "\n");
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
  ParameterDatabase db = TDomain::default_domain_parameters();
#ifdef __2D__
  db["geo_file"] = simplices ? "TwoTriangles" : "UnitSquare";
#else
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = simplices ? "Default_UnitCube_Tetra" : "Default_UnitCube_Hexa";
#endif
  db["conforming_closure"] = do_conforming_closure;
  
  TDomain domain(db);
  size_t n_adaptive_steps = 8;
  for(size_t curr_level = 0; curr_level < n_adaptive_steps; ++curr_level)
  {
    domain.RefineByIndicator(indicator);
    if(write_PS_and_VTU_files)
    {
      auto name = PS_name(curr_level, do_conforming_closure, simplices,
                          ansatz_order);
      domain.PS(name.c_str(), It_Finest, 0);
    }
  }
  // finally, solve two simple problems on this domain, only to check

  Chrono chrono;
  ConvDiff cd(domain, db);
  cd.assemble();
  cd.solve();
  cd.output();
  chrono.restart_and_print("solving dummy convection-diffusion problem");
  
  // Taylor-Hood elements
  if(ansatz_order == 5)
  {
    // hanging nodes for P6 are not supported, we solve a dummy problem instead
    ansatz_order = 1;
  }
  TDatabase::ParamDB->VELOCITY_SPACE = ansatz_order+1;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = 14;
  db.add("example", 2, "");
  Stokes stokes(domain, db);
  stokes.assemble_linear_terms();
  stokes.solve();
  stokes.output();
  chrono.restart_and_print("solving dummy Stokes problem");
  
  compare(cd, stokes, errors);
}


// =======================================================================
// main program
// =======================================================================
int main(int, char * argv[])
{
  parmoon::parmoon_initialize();
  bool testall = false;
  if (argv[1])
    testall = (std::string(argv[1]).compare("testall") == 0);
  
  Output::setVerbosity(2);
  
  std::array<double, 9> errors;
  bool do_conforming_closure;
  bool simplices;
  Chrono chrono_all;

  // ==========================================================================
  do_conforming_closure = true;
  simplices = true;
  errors = {{ 0.035770938439154,
              0.51608967114182, 0.51608967114182, 0.1410584617451,
              0.0033694651929196, 0.052739180565685, 0.10131900989499,
              0.061189866362196, 1.2606285806271}};
  test_indicator(do_conforming_closure, simplices, 1, errors);

  if(testall)
  {
    errors = {{ 0.0025402849496846,
                0.074968557695118, 0.074968557695118, 0.013210208083023,
                0.00015774407965653, 0.003025005307013, 0.006301255958945,
                0.004111734527521, 0.18407318866786 }};
    test_indicator(do_conforming_closure, simplices, 2, errors);

    errors = {{ 0.00017758187979405,
                0.0067654706203652, 0.0067654706203652, 0.0013618188684979,
                6.9047446340185e-06, 0.00018073621405227, 0.00039224779282074,
                0.00047269770076548, 0.027742712331764 }};
    test_indicator(do_conforming_closure, simplices, 3, errors);
    
    errors = {{ 1.4195556186299e-05,
                0.00066487351900192, 0.00066487351900192, 9.2887807569358e-05,
                2.1596605241172e-07, 6.2469912139526e-06, 1.4214603525187e-05,
                1.4601586281189e-05, 0.001362554550883 }};
    test_indicator(do_conforming_closure, simplices, 4, errors);
    
    errors = {{ 7.1275528937539e-07,
                3.9681826317622e-05, 3.9681826317622e-05, 6.0101510399313e-06,
                0.0033694651929196, 0.052739180565685, 0.10131900989499,
                0.061189866362196, 1.2606285806271 }};
    test_indicator(do_conforming_closure, simplices, 5, errors);
  }

  // ==========================================================================
  do_conforming_closure = false;
  simplices = true;
  errors = {{ 0.042380310691282,
              0.59129151149837, 0.59129151149837, 0.13495122911489,
              0.0040776830325563, 0.059620069966017, 0.11839522245947,
              0.062573587417992, 1.2777246574855 }};
  test_indicator(do_conforming_closure, simplices, 1, errors);

  if(testall)
  {
    errors = {{ 0.0027720696069861,
                0.081555432454459, 0.081555432454459, 0.012895816013502,
                0.00018256211506998, 0.0031675689911316, 0.0071989045250255,
                0.0041794269835082, 0.15759027621856 }};
    test_indicator(do_conforming_closure, simplices, 2, errors);

    errors = {{ 0.00020348123799333,
                0.0076917921372141, 0.0076917921372141, 0.0013602936366391,
                8.5639119491338e-06, 0.00019646294094567, 0.00046870777319375,
                0.00042894844203815, 0.024096330748339 }};
    test_indicator(do_conforming_closure, simplices, 3, errors);
    
    errors = {{ 1.5788182523811e-05,
                0.00073013922523117, 0.00073013922523117, 8.9588511335843e-05,
                2.6041352613623e-07, 6.5983274113732e-06, 1.6315134878125e-05,
                1.4121596969238e-05, 0.0011697764066892 }};
    test_indicator(do_conforming_closure, simplices, 4, errors);
    
    errors = {{ 8.1179355347452e-07, 
                4.4946135907158e-05, 4.4946135907158e-05, 5.8586848816766e-06,
                0.0040776830325563, 0.059620069966017, 0.11839522245947,
                0.062573587417992, 1.2777246574855 }};
    test_indicator(do_conforming_closure, simplices, 5, errors);
  }

  // ==========================================================================
  do_conforming_closure = true;
  simplices = false;
  errors = {{ 0.017283245814355,
              0.33311322626979, 0.33311322626979, 0.0762621835061,
              0.0018812306479113, 0.026942427615239, 0.055730550547945,
              0.019284748093558, 0.4513811252173 }};
  test_indicator(do_conforming_closure, simplices, 1, errors);

  if(testall)
  {
    errors = {{ 0.001279822868756,
                0.040891950610384, 0.040891950610384, 0.0073846377470416,
                0.00012203730364305, 0.0011495175851034, 0.0042198097032549,
                0.0020200030697402, 0.071262662273862 }};
    test_indicator(do_conforming_closure, simplices, 2, errors);

    errors = {{ 5.7208722496519e-05,
                0.0027018884574178, 0.0027018884574178, 0.00045784007898819,
                2.7533597086386e-06, 7.6842358375396e-05, 0.00015355517166441,
                9.9193006541922e-05, 0.0052795245628617 }};
    test_indicator(do_conforming_closure, simplices, 3, errors);

    errors = {{ 3.6366401311769e-06,
                0.00020633147020743, 0.00020633147020743, 2.670374549063e-05,
                1.4662623966544e-07, 1.9854294815438e-06, 8.0023460585573e-06,
                4.7196360998819e-06, 0.00034253670668658 }};
    test_indicator(do_conforming_closure, simplices, 4, errors);

    errors = {{ 1.278483456398e-07,
                8.9801283445683e-06, 8.9801283445683e-06, 1.1202919092679e-06,
                0.0018812306479113, 0.026942427615239, 0.055730550547945,
                0.019284748093558, 0.4513811252173 }};
    test_indicator(do_conforming_closure, simplices, 5, errors);
  }

  // ==========================================================================
  do_conforming_closure = false;
  simplices = false;
  errors = {{ 0.011874554092922,
              0.25451703865297, 0.25451703865297, 0.03750007808692,
              0.0031073214674175, 0.035484814538899, 0.082583770670472,
              0.013251140025838, 0.40341076538356 }};
  test_indicator(do_conforming_closure, simplices, 1, errors);

  if(testall)
  {
    errors = {{ 0.001153069846948,
                0.032083741477288, 0.032083741477288, 0.0032514526789689,
                0.00015446429247842, 0.0019587821219028, 0.0056359926662344,
                0.0012567553654574, 0.041291059431891 }};
    test_indicator(do_conforming_closure, simplices, 2, errors);

    errors = {{ 3.0346622444448e-05,
                0.0012592064796466, 0.0012592064796466, 0.0001375525439139,
                5.1334500295697e-06, 0.00010711536174364, 0.0002700320508768,
                8.3942466748833e-05, 0.0040347672657818 }};
    test_indicator(do_conforming_closure, simplices, 3, errors);

    errors = {{ 1.9880384093545e-06,
                0.00010071722035562, 0.00010071722035562, 7.4341387803489e-06,
                1.9387173962012e-07, 3.6907978999865e-06, 1.1300080311758e-05,
                4.2880353077882e-06, 0.00026178663645479 }};
    test_indicator(do_conforming_closure, simplices, 4, errors);

    errors = {{ 3.6993351781254e-08,
                2.3517749820065e-06, 2.3517749820065e-06, 2.0142433554415e-07,
                0.0031073214674175, 0.035484814538899, 0.082583770670472,
                0.013251140025838, 0.40341076538356 }};
    test_indicator(do_conforming_closure, simplices, 5, errors);
  }

  chrono_all.stop_and_print("full test 'adaptive_hanging' ");
  parmoon::parmoon_finalize();
}


std::string PS_name(size_t curr_level, bool do_conforming_closure,
                    bool simplices, int ansatz_order)
{
  std::string name{"domain_indicator_after"};
  name += do_conforming_closure ? "_" : "_no";
  name += "ConfClosure";
  name += simplices ? "_tria" : "_quad"; // writing PS makes little sense in 3D
  name += "_order" + std::to_string(ansatz_order);
  name += "_level" + std::to_string(curr_level);
  name += ".ps";
  return name;
}
