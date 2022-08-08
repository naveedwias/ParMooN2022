/**
 * @brief A test program for the checking the slope limiter for FE functions.
 *
 * This serves as a test program for TFEFunction2D::apply_limiter(). It tests
 * several setups: two functions with known limiters, a piecewise constant
 * function, two monotonic functions, a linear function and some reference
 * functions to check if the behaviour of the program has changed. Last but not
 * least, it is checked that if the limiter is apllied twice the result after
 * the first application did not change anymore.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then
 * this tests  shows you how many other program parts are affected by your
 * changes. If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 * See also: cd2d_test.c++
 *
 *
 * @author Derk Frerichs
 */

#include "ConvectionDiffusion.h"
#include "ParMooN_repository_info.h"
#include "Database.h"
#include "ParMooN.h"

// #include<math.h>

/* ########################################################################## */
/* Declare functions and constant expressions, see below for definition */
/* ########################################################################## */
void perform_test(bool testall, ParameterDatabase parmoon_db);

void check_errors(ParameterDatabase parmoon_db);
std::array<double, 5> get_reference_values(const int& domain, const
    std::string& limiter, const int& lvl, const int& degree);
void compare_errors(std::string error, double ref_err, double new_err);

std::vector<double> copy_fe_function(ConvectionDiffusion<2>* cd2d);
bool are_fe_functions_equal(std::vector<double> values_orig,
    double* values_new);

/* ########################################################################## */
/* Declare some often used variables */
/* ########################################################################## */
const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";
const int max_n_refinements = 4;


/* ########################################################################## */
/* Main program */
/* ########################################################################## */
int main(int, char* argv[])
{
  parmoon::parmoon_initialize();
  bool testall = false;
  if (argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }
  ParameterDatabase parmoon_db =
    ConvectionDiffusion<2>::default_cd_database(true);
  TDatabase::ParamDB->USE_ISOPARAMETRIC = 0;

  // BE CAREFUL:
  // The global parameter database can be really a pain. In this test I use
  // different examples. Some of these examples change the global parameter
  // INTERNAL_QUAD_RULE. For me this resulted in the fact that the ordering of
  // my experiments changed the norms of the solution since different quadrature
  // rules were used. This took me around 2 days together with Ulrich to figure
  // out where this error comes from. Without global databases this error may
  // have been circumvented.
  // Lessons learned: Be careful with experiments which set parameters in the
  // global database.
  //
  // I set now the quadrature rule to 97 to guarantee that in all
  // experiments this quadrature rule is used.
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;


  parmoon_db["space_discretization_type"].set("dg", true);
  parmoon_db["apply_limiter"].set(true, true);
  parmoon_db["M_lim"].set(1, true);
  parmoon_db["gamma_limiter"].set(1, true);
  perform_test(testall, parmoon_db);

  Output::print("\n\n-----------------------------");
  Output::print("         TEST PASSED\n");
  parmoon::parmoon_finalize();
}

/* ########################################################################## */
/* Declaration of functions */
/* ########################################################################## */
void perform_test(bool testall, ParameterDatabase parmoon_db)
{
  // This function calls all the single tests.
  parmoon_db["alpha_ref"] = 4;
  std::vector<std::string> geometries = {"UnitSquare_quads.mesh",
    "UnitSquare_tria.mesh", "UnitSquareIrregular.GEO", "UnitSquareWithSpherical"
      "InscribedRegion.mesh"};
  std::vector<std::string> boundaries(geometries.size(), "UnitSquare.PRM");
  if (!testall)
  {
    geometries.resize(2);
    boundaries.resize(geometries.size());
  }
  int max_approx_degree = (testall) ? 4 : 2;
  parmoon_db.add("max_degree", max_approx_degree, " ");

  int start_geom_i = 0;
  parmoon_db.add("geometry_nr", start_geom_i, " ", start_geom_i, (int)
      geometries.size());

  for (unsigned int geom_i = start_geom_i; geom_i < geometries.size(); ++geom_i)
  {
    parmoon_db["geometry_nr"].set((int) geom_i, true);
    parmoon_db["geo_file"].set(path_to_meshes+geometries[geom_i], false);
    parmoon_db["boundary_file"].set(path_to_meshes+boundaries[geom_i], false);
    Output::print("\nStart checking geometry ", geom_i, "\n");

    std::vector<std::string> limiters = {"LinTriaReco", "ConstTriaReco",
      "LinQuadReco", "ConstQuadReco", "LinQuadDeriv", "ConstQuadDeriv",
      "ConstJump", "ConstJumpMod"};

    for (auto method : limiters)
    {
      if ((method == "LinQuadDeriv" || method == "ConstQuadDeriv") &&
          geom_i != 0)
      {
        // These limiter do only make sense on quadliteral meshes
        continue;
      }
      Output::print("Start checking limiter ", method);
      parmoon_db["limiter_name"].set(method, true);
      check_errors(parmoon_db);

      Output::print('\n');
    }
  }
}

void check_errors(ParameterDatabase parmoon_db)
{
  // Get some values of the specific problem
  const auto domain_nr = parmoon_db["geometry_nr"].get<int>();
  const auto limiter_name = parmoon_db["limiter_name"].get<std::string>();

  // Set DG method for which slope limiters are designed and the problem
  parmoon_db["symmetry_DG"].set(-1, true);
  parmoon_db["face_sigma_DG"].set(.5e-8, true);
  parmoon_db["diffusion_coefficient"].set(1e-8, true);
  parmoon_db["example"].set(5,true);
  if (limiter_name == "LinQuadDeriv" || limiter_name == "ConstQuadDeriv")
  {
    // For this limiter another test is chosen since in the other problem this
    // limiter is not active
    parmoon_db["example"].set(6, true);
  }
  TDomain domain(parmoon_db);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  auto max_ref = std::min(max_n_refinements, 3);

  // Loop over refinements and polynomial degrees, compute norms and compare to
  // reference values
  domain.RegRefineAll();
  for (int ref_i = 1; ref_i < max_ref; ++ref_i)
  {
    for (int degree = 1; degree <=parmoon_db["max_degree"].get<int>(); ++degree)
    {
      TDatabase::ParamDB->ANSATZ_ORDER = -10-degree;

      ConvectionDiffusion<2> cd2d(domain, parmoon_db);
      cd2d.assemble();
      cd2d.solve();
      cd2d.output();

      // // With this code snippet a text file called "reference_norms.txt" is
      // // created where the norms of the problems are listed. The ordering of the
      // // norms is exactly as in the function get_reference_values(). I used this
      // // snippet together with some copy feature of my editor to easily fill the
      // // values in get_reference_values.
      // // Note, that this snippet APPENDS the values to the file
      // // reference_norms.txt. So you might want to delete it before using this
      // // code.
      // if (limiter_name != "LinQuadReco" && limiter_name != "ConstQuadReco")
      // {
      //   std::ofstream myfile;
      //   myfile.open ("reference_norms.txt", std::ios::app);
      //   myfile << std::setprecision(14) << cd2d.get_L2_error() << "," <<
      //     cd2d.get_H1_semi_error() << "," << cd2d.get_SD_error() << "," <<
      //     cd2d.get_DG_error() << "," << cd2d.get_L_inf_error() << "\n";
      //   myfile.close();
      // }

      // The limiter LinQuadReco and ConstQuadReco are nothing else than the limiter
      // LinTriaReco and ConstTriaReco. Therefore, a modified limiter_name is
      // created where the first mentioned limiter are mapped to the last mentioned
      // names.
      auto lim_name_mod = limiter_name;
      if (lim_name_mod == "LinQuadReco")
      {
        lim_name_mod = "LinTriaReco";
      }
      else if (lim_name_mod == "ConstQuadReco")
      {
        lim_name_mod = "ConstTriaReco";
      }
      auto errors = get_reference_values(domain_nr, lim_name_mod, ref_i,
          degree);
      compare_errors("L2", errors[0], cd2d.get_L2_error());
      compare_errors("H1-semi", errors[1], cd2d.get_H1_semi_error());
      compare_errors("SD", errors[2], cd2d.get_SD_error());
      compare_errors("DG", errors[3], cd2d.get_DG_error());
      compare_errors("L_inf", errors[4], cd2d.get_L_inf_error());
    }
    domain.RegRefineAll();
  }
  Output::print("Checked reference errors.");
}

void compare_errors(std::string error, double ref_err, double new_err)
{
  double eps = 1e-11;
  if( std::abs(ref_err - new_err) > eps )
  {
    ErrThrow(error, " error not correct. Computed: ",std::setprecision(14),
        new_err, ", reference: ", ref_err);
  }
}

std::array<double, 5> get_reference_values(const int& domain_nr, const
    std::string& limiter, const int& lvl, const int& degree)
{
  std::array<double, 5> error;
  switch (domain_nr)
  {
    case 0:
      if (limiter == "LinTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90067702177776,1.5752781618149,0.14575892751457,1.0005814404124,1.2030336019135};
                break;
              case 2:
                error = {0.92483088816734,2.4213326304871,0.26008367981119,1.0187596450685,1.3975169895754};
                break;
              case 3:
                error = {0.92074809662053,2.9978783537281,0.27878740515785,1.0129791405938,1.3811380697513};
                break;
              case 4:
                error = {0.9185028770492,3.4661930483607,0.29722333524095,1.0101973464975,1.3267939427059};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92202322334876,2.0715965245894,0.16630438175031,1.017409134129,1.3236646184256};
                break;
              case 2:
                error = {0.9172300966739,2.5465840481146,0.20265311019427,1.0227976758389,1.3975169027333};
                break;
              case 3:
                error = {0.92411920977256,3.0625156155674,0.2154920840282,1.035310111674,1.1960393409583};
                break;
              case 4:
                error = {0.92035280263915,3.5605608089947,0.21470118573541,1.0291522276406,1.192119482894};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else if (limiter == "ConstTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.89771094217658,1.1917443849423,0.10105414584775,1.0402655782962,1.061862099707};
                break;
              case 2:
                error = {0.91322564683252,9.3631498744457e-16,2.5042558772364e-16,1.0939834766004,1.0168371406548};
                break;
              case 3:
                error = {0.90611206689807,1.3488375704881e-15,3.8010388304055e-16,1.088491314908,1.0014290161463};
                break;
              case 4:
                error = {0.90664775854678,1.6937310814317,0.063937885160751,1.0665561938012,1.1776237561592};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92097550111603,1.6312090858758,0.077079903589412,1.0402439929557,1.1561684504222};
                break;
              case 2:
                error = {0.91396030498495,0.012342939398234,0.0010775566660273,1.0915985258115,1.0343647029061};
                break;
              case 3:
                error = {0.92316020471367,2.6896150640903,0.17756691986498,1.0680692019507,1.1960393409583};
                break;
              case 4:
                error = {0.92020833452276,3.4359538082253,0.21162417390045,1.06306556906,1.192119482894};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "LinQuadDeriv")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.27728683868785,2.3399843504679,1.3803141798212,0.43343787445925,0.72685370653841};
                break;
              case 2:
                error = {0.2741452426316,2.4080765543204,1.3527597987165,0.30689647908727,0.91190563619661};
                break;
              case 3:
                error = {0.2752707414325,2.5526895497731,1.321391024562,0.26522237402627,0.79279608758951};
                break;
              case 4:
                error = {0.27616439696742,2.946660667423,1.4008469590738,0.21670736828167,0.84679703394271};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.27358185901936,1.4074941805453,0.67146818059226,0.24535825538692,0.71957062705551};
                break;
              case 2:
                error = {0.27265924037296,1.7135969503664,0.67791058955989,0.21963012764213,0.7243251259149};
                break;
              case 3:
                error = {0.27335795235105,1.9699941905979,0.6507350767781,0.1886122562543,0.73205829672071};
                break;
              case 4:
                error = {0.27353184506941,2.1812764195907,0.65587377018579,0.15809782638112,0.75006696440514};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstQuadDeriv")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.27728683868785,2.3399843504679,1.3803141798212,0.43343787445925,0.72685370653841};
                break;
              case 2:
                error = {0.2741452426316,2.4080765543204,1.3527597987165,0.30689647908727,0.91190563619661};
                break;
              case 3:
                error = {0.2752707414325,2.5526895497731,1.321391024562,0.26522237402627,0.79279608758951};
                break;
              case 4:
                error = {0.27616439696742,2.946660667423,1.4008469590738,0.21670736828167,0.84679703394271};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.27358185901936,1.4074941805453,0.67146818059226,0.24535825538692,0.71957062705551};
                break;
              case 2:
                error = {0.27265924037296,1.7135969503664,0.67791058955989,0.21963012764213,0.7243251259149};
                break;
              case 3:
                error = {0.27335795235105,1.9699941905979,0.6507350767781,0.1886122562543,0.73205829672071};
                break;
              case 4:
                error = {0.27353184506941,2.1812764195907,0.65587377018579,0.15809782638112,0.75006696440514};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJump")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.9013678539929,1.6581226168701,0.17063831424449,0.99550261388981,1.2030336019135};
                break;
              case 2:
                error = {0.92505946123564,2.4820834880716,0.27023049317597,1.0147803250582,1.3975169895754};
                break;
              case 3:
                error = {0.92105861706005,3.1247296575907,0.28996793392293,1.0076498130241,1.3811380697513};
                break;
              case 4:
                error = {0.91886368955519,3.6623795731869,0.31236157545053,1.0038068280923,1.3267939427059};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92233844308042,2.1925719417929,0.19251075421185,1.0126817611495,1.3236646184256};
                break;
              case 2:
                error = {0.91944611475653,3.2996027737144,0.20659957071131,1.0051877081665,1.3975169027333};
                break;
              case 3:
                error = {0.92710340742257,4.0940844291701,0.21692684892085,1.0119806735009,1.1960393409583};
                break;
              case 4:
                error = {0.9237602307118,4.7924282717554,0.23272071970405,1.0073492662198,1.2039214516844};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJumpMod")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.89602643397125,0.90948360157965,0.080608169949882,1.0482689928919,1.1376828397671};
                break;
              case 2:
                error = {0.92179482013426,2.0060558574594,0.15419533500559,1.0404970980271,1.208762253717};
                break;
              case 3:
                error = {0.91655705735273,2.4472310071776,0.13919211452262,1.0375584929449,1.1981569644408};
                break;
              case 4:
                error = {0.91886368955519,3.6623795731869,0.31236157545053,1.0038068280923,1.3267939427059};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.91914068754568,0.17548497588279,0.025026106607363,1.0811847867813,1.1228086759861};
                break;
              case 2:
                error = {0.91472474460119,1.2439862091345,0.042414682358342,1.075913515949,1.1315408570026};
                break;
              case 3:
                error = {0.92074253441083,0.51273006879191,0.041654754907073,1.099014897672,1.1658153891537};
                break;
              case 4:
                error = {0.91738195097424,1.9056932788075,0.055000100666434,1.0805334099579,1.1583638223766};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else
      {
        ErrThrow("Unknown limiter in get_reference_values");
      }
      break;
    case 1:
      if (limiter == "LinTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90515690175852,1.7185369589499,0.17244552985564,0.99344799285992,1.3475455861269};
                break;
              case 2:
                error = {0.93157376841225,1.7219724099166,0.16388503028491,1.0453453404758,1.1528158225313};
                break;
              case 3:
                error = {0.90372506763915,2.3363363277981,0.23817328721676,1.0171792449404,1.1351873519538};
                break;
              case 4:
                error = {0.92325113709458,2.6440555315613,0.26129832664508,1.0398490196941,1.1763439505752};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90576461249383,1.9251899011445,0.063538910122629,1.000592604038,1.1689148235455};
                break;
              case 2:
                error = {0.92561554533425,2.6697521478882,0.15138723513786,1.0174006874391,1.2909606603809};
                break;
              case 3:
                error = {0.91977245161092,2.9901616977032,0.15113404688542,1.0147915147492,1.2874664111515};
                break;
              case 4:
                error = {0.91679649245065,3.4169023420634,0.19325804536748,1.0114793173686,1.248294620897};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else if (limiter == "ConstTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.89780011350307,1.0691213185771,0.0001069121352272,1.0399631182169,1.3399562622286};
                break;
              case 2:
                error = {0.92522559939361,7.2952500582364e-07,2.3565189148639e-07,1.0405497280775,1.0035190895287};
                break;
              case 3:
                error = {0.89345480208029,2.2432091916249e-06,6.295499844022e-07,1.0279059921978,1.0095170876596};
                break;
              case 4:
                error = {0.91462527348848,1.6918986132243e-14,5.5049969194921e-15,1.0351522560981,1.0003974110551};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.9049961320888,1.7526328022571,0.02915874124701,1.0047508936021,1.1689148235455};
                break;
              case 2:
                error = {0.9210263792404,1.3325138093195,0.040376591161716,1.0334181635672,1.1120969842799};
                break;
              case 3:
                error = {0.91617998374203,1.7858178885355,0.00017859075073945,1.0296475592642,1.1362664622467};
                break;
              case 4:
                error = {0.9127696453162,1.7532843094878,0.00017532851835421,1.02752074608,1.1265671529754};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJump")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90515690175852,1.7185369589499,0.17244552985564,0.99344799285992,1.3475455861269};
                break;
              case 2:
                error = {0.93531380845785,2.1573531760307,0.16370322862808,1.0219780682751,1.1528158225313};
                break;
              case 3:
                error = {0.90947588918337,2.7251058303026,0.23681192924149,0.99469762796524,1.1351873519538};
                break;
              case 4:
                error = {0.92762866219769,3.0579908855699,0.26116591761826,1.0126730145074,1.1763439505752};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90460866096828,1.6637432050438,0.073042536933907,1.0147416277667,1.1726528163721};
                break;
              case 2:
                error = {0.92601837085162,2.8095812829956,0.15326778775339,1.0122631880283,1.2909606603809};
                break;
              case 3:
                error = {0.92208878735113,3.5362198067076,0.15408903689696,1.0060677770563,1.2874664111515};
                break;
              case 4:
                error = {0.91906189779177,3.9458438115594,0.19926582632295,1.0027321455305,1.248294620897};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJumpMod")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90515690175852,1.7185369589499,0.17244552985564,0.99344799285992,1.3475455861269};
                break;
              case 2:
                error = {0.93531380845785,2.1573531760307,0.16370322862808,1.0219780682751,1.1528158225313};
                break;
              case 3:
                error = {0.90947588918337,2.7251058303026,0.23681192924149,0.99469762796524,1.1351873519538};
                break;
              case 4:
                error = {0.92762866219769,3.0579908855699,0.26116591761826,1.0126730145074,1.1763439505752};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90460866096828,1.6637432050438,0.073042536933907,1.0147416277667,1.1726528163721};
                break;
              case 2:
                error = {0.92331338814784,2.1493941833356,0.075590262005948,1.0299486079391,1.1120969842799};
                break;
              case 3:
                error = {0.92208878735113,3.5362198067076,0.15408903689696,1.0060677770563,1.2874664111515};
                break;
              case 4:
                error = {0.91906189779177,3.9458438115594,0.19926582632295,1.0027321455305,1.248294620897};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else
      {
        ErrThrow("Unknown limiter in get_reference_values");
      }
      break;
    case 2:
      if (limiter == "LinTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90451236024853,1.3694016918762,0.13633828158337,1.0161849309447,1.3727602665855};
                break;
              case 2:
                error = {0.93260611608361,1.5031679168513,0.15147911029967,1.0473174713678,1.1910556438901};
                break;
              case 3:
                error = {0.90565633084646,1.9053071393393,0.21866203334069,1.0159541990865,1.1746289776736};
                break;
              case 4:
                error = {0.92639067064333,2.6770740176127,0.2437175614081,1.0318004101081,1.3886154135979};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.9047921565775,1.7033397203601,0.10376301153674,1.0243822431461,1.174840310267};
                break;
              case 2:
                error = {0.92448211836621,2.1025844469941,0.19916369155049,1.062679402593,1.38504475677};
                break;
              case 3:
                error = {0.91979120662353,2.3981473307851,0.17058819858159,1.0563190609648,1.38769508584};
                break;
              case 4:
                error = {0.91680454001936,2.4648897967686,0.15641609695291,1.0641555360361,1.3703278540782};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else if (limiter == "ConstTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90111562135055,1.0662167897385,0.026019737442396,1.0599684025917,1.3594487831224};
                break;
              case 2:
                error = {0.9279952221145,0.56648487093553,0.088211946184057,1.0510576210398,1.0764856506601};
                break;
              case 3:
                error = {0.90183785880622,1.3532886877539,0.021880689988606,1.0446205981515,1.0733065235768};
                break;
              case 4:
                error = {0.92261391091808,1.5408216528389,0.00015408251583406,1.0475791032506,1.1494860644607};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90435130341067,1.5829559805147,0.088796885216417,1.0334033880979,1.174840310267};
                break;
              case 2:
                error = {0.92207450574026,1.3491111082006,0.011972041554096,1.0655798555001,1.0703833256206};
                break;
              case 3:
                error = {0.91787457297128,1.9596491234618,0.042903721836248,1.0484622779424,1.0963102500819};
                break;
              case 4:
                error = {0.91525722220984,2.4157164662627,0.069184659417499,1.0568647656052,1.122097937929};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJump")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.90747113260849,1.6067957874325,0.17852990440084,0.99395759870136,1.3727602665855};
                break;
              case 2:
                error = {0.93592409932973,2.2065173772896,0.17784928025753,1.0222868843463,1.1910556438901};
                break;
              case 3:
                error = {0.91114953879847,2.8079411219611,0.23105014576999,0.99497280486527,1.1746289776736};
                break;
              case 4:
                error = {0.92902594209001,3.3318446248894,0.25871249688354,1.0128659481669,1.2766601359111};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90536558152621,1.7473103413665,0.076362100493482,1.0191068335341,1.1883441052385};
                break;
              case 2:
                error = {0.9255402999411,2.6477898591936,0.10442438693998,1.0246216131138,1.3040298474693};
                break;
              case 3:
                error = {0.92114644150128,3.3779090328081,0.12140243431823,1.018344386213,1.1711772995234};
                break;
              case 4:
                error = {0.91987007837261,4.3533784491082,0.19378500932416,1.0027941658834,1.3703278540782};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJumpMod")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.89762518568754,3.1645485968137e-06,6.2937916740872e-07,1.0664960527393,1.0368524781562};
                break;
              case 2:
                error = {0.93140537559171,1.2242392060431,0.066810512878908,1.0597467139064,1.1910556438901};
                break;
              case 3:
                error = {0.90777412850526,2.3611871014598,0.083940394234566,1.0388252921609,1.1728751990648};
                break;
              case 4:
                error = {0.92367989387588,2.3917190636543,0.070751393352201,1.0475932230337,1.1494860644607};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.90259557248492,0.10510249467079,0.00035327890037192,1.0510638928228,1.0240899389815};
                break;
              case 2:
                error = {0.92159525049822,0.57136004205759,0.016784142996461,1.0699851998851,1.0566166845425};
                break;
              case 3:
                error = {0.91702336017954,1.9780393320298,0.025816456087747,1.0595082567767,1.0463431247567};
                break;
              case 4:
                error = {0.91357525030629,1.6543564937792,0.018839580320598,1.0570175257423,1.0471865907825};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else
      {
        ErrThrow("Unknown limiter in get_reference_values");
      }
      break;
    case 3:
      if (limiter == "LinTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.92368994787467,3.2179492984565,0.091820377404902,1.0253357669811,1.1740561186202};
                break;
              case 2:
                error = {0.92419215775554,3.6729598080818,0.096216481543991,1.0420804571558,1.4166709024905};
                break;
              case 3:
                error = {0.92452225726646,4.2171140428035,0.10492880246907,1.0431166225784,1.3478985062466};
                break;
              case 4:
                error = {0.92466378684472,4.7141600577474,0.11836771768207,1.0487363690339,1.3055566335578};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.9247207957591,4.1822620846832,0.039865897540068,1.0197993096424,1.2312169519479};
                break;
              case 2:
                error = {0.92550644532197,5.4753739296368,0.054632277327378,1.0287881043003,1.27172247735};
                break;
              case 3:
                error = {0.92561076592591,5.6515257299632,0.06006472915242,1.0393574583894,1.3504110819516};
                break;
              case 4:
                error = {0.92562680481279,5.6707601464527,0.079151274805569,1.0479493670665,1.5019243036894};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else if (limiter == "ConstTriaReco")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.92348290294181,2.4307482696941,0.032421586599118,1.0399738217399,1.1591171005769};
                break;
              case 2:
                error = {0.92396664788708,2.576274412208,0.013234089882576,1.0516736330367,1.1714963726555};
                break;
              case 3:
                error = {0.9243860596621,3.2416742460562,0.077195135599969,1.0550769959145,1.2377931063073};
                break;
              case 4:
                error = {0.92505484940331,4.989307049933,0.075612671652871,1.0485615929363,1.2097394159599};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92469374531352,3.8358528425452,0.025421976758009,1.0243237696149,1.2312169519479};
                break;
              case 2:
                error = {0.92542446238139,4.5479234588838,0.014713132632087,1.0363326488919,1.0895684689077};
                break;
              case 3:
                error = {0.92554196992572,4.6537500601324,0.011044767313145,1.0454767094288,1.0796286197357};
                break;
              case 4:
                error = {0.9255880873084,4.8313167564481,0.015872305377248,1.0589508829149,1.0946386857492};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJump")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.92357726055576,2.6384126383792,0.029555923981586,1.0392997346195,1.1664422507503};
                break;
              case 2:
                error = {0.92432625188211,3.8921000028337,0.027578301672559,1.0330975850478,1.2961252910136};
                break;
              case 3:
                error = {0.92496929183209,4.9543135420868,0.038211332111361,1.0243225606377,1.2470480642781};
                break;
              case 4:
                error = {0.92511708426513,5.7382819499796,0.029488820321712,1.0339659664559,1.1369494330467};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92459139904562,2.9715843839577,0.019327003386573,1.030325070689,1.1085399240202};
                break;
              case 2:
                error = {0.92546746100693,5.0656947867799,0.023924669385826,1.0287374891669,1.0895684689077};
                break;
              case 3:
                error = {0.92571912786252,6.3587145599471,0.019183360002138,1.0256042157884,1.1087343804843};
                break;
              case 4:
                error = {0.92592020912667,7.5717837254072,0.022446871722961,1.0257593835859,1.0955005168542};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }

      }
      else if ( limiter == "ConstJumpMod")
      {
        switch (lvl)
        {
          case 1:
            switch (degree)
            {
              case 1:
                error = {0.92319570282797,0.14340397300844,0.0024366759397559,1.0486876858536,1.0569346445811};
                break;
              case 2:
                error = {0.92375847474813,0.81464591412187,0.0029998989612368,1.0599453078948,1.0668677577835};
                break;
              case 3:
                error = {0.92411461853778,2.0738552918849,0.0070602840248946,1.0574580369331,1.0606967159616};
                break;
              case 4:
                error = {0.9244392970578,3.3641079456262,0.0099992604214746,1.0583302913574,1.0948387864489};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          case 2:
            switch (degree)
            {
              case 1:
                error = {0.92451006581576,0.11439706862164,0.0012696409354076,1.0388557570207,1.0539593314599};
                break;
              case 2:
                error = {0.9252109659,0.14211721510113,0.0017347913601593,1.0521858076272,1.0737312968661};
                break;
              case 3:
                error = {0.92534896950571,1.7451731803474,0.002530953669586,1.0576853818057,1.0686977101649};
                break;
              case 4:
                error = {0.92543446186421,2.7975017610963,0.0039720297918225,1.0601401857922,1.0821683320394};
                break;
              default:
                ErrThrow("Unknown polynomial degree in get_reference_values");
            }
            break;
          default:
            ErrThrow("Unknown refinement lvl in get_reference_values");
        }
      }
      else
      {
        ErrThrow("Unknown limiter in get_reference_values");
      }
      break;
    default:
      ErrThrow("Unknown domain_nr for reference values.");
  }
  return error;
}
