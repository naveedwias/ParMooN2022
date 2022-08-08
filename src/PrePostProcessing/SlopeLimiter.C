#include <algorithm>
#include <numeric>
#include <string>
#include <cmath>
#include "Database.h"
#include "FEDatabase.h"
#include "BaseCell.h"
#include "MainUtilities.h"
#include "MooNMD_Io.h"
#include "Utilities.h"
#include "SlopeLimiter.h"
#include "QuadratureFormulaDatabase.h"

void TSlopeLimiter::check_limiter_name(const std::string& name_of_limiter)
{
  Output::print<5>("Verify that the given limiter is implemented at all.");
  // Check if it is an allowed limiter
  auto allowed_limiters = possible_limiters;
  for (auto limiter : allowed_limiters)
  {
    if (name_of_limiter == limiter)
    {
      // The given limiter is really implemented
      break;
    }
    else if (limiter == allowed_limiters.back())
    {
      // The limiter does not exist since we checked the last element in the
      // array
      ErrThrow("Unknown limiter! You have to specify a valid limiter, see also",
          " documentation -> TSlopeLimiter");
    }
  }
}

TSlopeLimiter::TSlopeLimiter(const ParameterDatabase& param_db)
{
  Output::print<5>("Construct a limiter using given database.");
  auto db = default_slope_limiter_database();
  db.merge(param_db, false, false);

  apply_limiter = db["apply_limiter"];

  limiter_name = db["limiter_name"].get<std::string>();
  check_limiter_name(limiter_name);

  set_m_lim(db["M_lim"]);
  set_gamma(db["gamma_limiter"]);
  set_alpha_ref(db["alpha_ref"]);
  set_char_length(db["characteristic_length"]);
  set_char_solution_scale(db["characteristic_solution_scale"]);
  set_C0_CJM(db["C0_CJM"]);
  set_C1_CJN(db["C1_CJN"]);
}


void TSlopeLimiter::set_m_lim(const double& value)
{
  Output::print<5>("Set M_lim.");
  if (value < 0)
  {
    ErrThrow("M_lim has to be greater or equal to 0.");
  }
  else
  {
    m_lim = value;
  }
}


void TSlopeLimiter::set_gamma(const double& value)
{
  Output::print<5>("Set gamma.");
  if (value < 0)
  {
    ErrThrow("gamma has to be greater or equal to 0.");
  }
  else
  {
    gamma = value;
  }
}


void TSlopeLimiter::set_alpha_ref(const double& value)
{
  Output::print<5>("Set alpha_ref.");
  if (value < 1 || value > 5)
  {
    ErrThrow("alpha_ref has to be in the interval [1, 5].");
  }
  else
  {
    alpha_ref = value;
  }
}

void TSlopeLimiter::set_char_length(const double& value)
{
  Output::print<5>("Set characteristic_length.");
  if (value <= 0)
  {
    ErrThrow("characteristic_length has to be greater than 0.");
  }
  else
  {
    characteristic_length = value;
  }
}

void TSlopeLimiter::set_char_solution_scale(const double& value)
{
  Output::print<5>("Set characteristic_solution_scale.");
  if (value == 0)
  {
    ErrThrow("characteristic_solution_scale has to be non 0.");
  }
  else
  {
    characteristic_solution_scale = value;
  }
}

void TSlopeLimiter::set_C0_CJM(const double& value)
{
  Output::print<5>("Set C0_CJM.");
  if (value < 0)
  {
    ErrThrow("C0_CJM has to be greater than 0.");
  }
  else
  {
    C0_CJM = value;
  }
}

void TSlopeLimiter::set_C1_CJN(const double& value)
{
  Output::print<5>("Set C1_CJN.");
  C1_CJN = value;
}


ParameterDatabase TSlopeLimiter::default_slope_limiter_database()
{
  Output::print<5>("\nCreating default slope limiter database:");
  ParameterDatabase db("Slope Limiter database");;

  Output::print<5>("Add apply_limiter");
  std::string description = "This boolean parameter determines whether a slope "
    "limiter shall be applied to the solution after the problem is solved. "
    "The limiter is chosen with database[\"limiter_name\"]";
  db.add<bool>("apply_limiter", false, description);

  Output::print<5>("Add limiter_name");
  description = "This string specifies the limiter that is applied if "
    "database[\"apply_limiter\"] is true. For a description of the limiter see "
    "TSlopeLimiter.";
  std::set<std::string> range(std::begin(possible_limiters),
      std::end(possible_limiters));
  db.add<std::string>("limiter_name", "Galerkin", description, range);

  Output::print<5>("Add M_lim");
  description = "This value is the maximal allowed value for the absolute value"
    " of the xi and eta coefficients of the FE function. It is only needed for "
    "the limiter \"LinQuadDeriv\" and \"ConstQuadDeriv\", see TSlopeLimiter.";
  db.add("M_lim", 1., description, 0., 1.e100);

  Output::print<5>("Add gamma_limiter");
  description = "This value is the factor of decreasing the difference of the "
    "constant coefficients in the minmod case of the limiters "
    "\"LinQuadDeriv\" and \"ConstQuadDeriv\", see TSlopeLimiter.";
  db.add("gamma_limiter", 1., description, 0., 1.e100);

  Output::print<5>("Add alpha_ref");
  description = "When \"ConstJumpMod\" is used as the limiter, the cell gets "
    "marked if alpha_E < alpha_ref, where alpha_E is the alpha computed for "
    "each edge of the cell, see TSlopeLimiter.";
  db.add("alpha_ref", 4., description, 1., 5.);

  Output::print<5>("Add characteristic_length");
  description = "Characteristic length of the domain used in ConstJumpMod. "
    "This should be larger than the longest edge in the triangulation. See "
    "also TSlopeLimiter.";
  db.add("characteristic_length", 1., description, 0., 1.e100);

  Output::print<5>("Add characteristic_solution_scale");
  description = "Characteristic measure of the solution used in ConstJumpMod, "
    "see also TSlopeLimiter.";
  db.add("characteristic_solution_scale", 1., description, -1e100, 1.e100);

  Output::print<5>("Add C0_CJM");
  description = "Reference constant in ConstJumpMod, see also TSlopeLimiter.";
  db.add("C0_CJM", 1., description, 0., 1.e100);

  Output::print<5>("Add C1_CJN");
  description = "Constant in ConstJumpNorm, see also TSlopeLimiter.";
  db.add("C1_CJN", 0., description, -1.e100, 1.e100);

  return db;
}

std::array<std::string, 12> TSlopeLimiter::possible_limiters = {"Galerkin",
  "LinTriaReco", "ConstTriaReco", "LinQuadReco", "ConstQuadReco",
  "LinQuadDeriv", "ConstQuadDeriv", "ConstJump", "ConstJumpMod",
  "ConstJumpL1Norm", "ConstJumpL2Norm", "ConstJumpLinftyNorm"};


std::vector<bool> TSlopeLimiter::get_is_cell_to_limit(const FEFunction*
    fe_function)
{
  // The cells have to be computed in case they haven't been computed yet or
  // with a different fe_function. To this extend also the features have to be
  // computed.
  compute_features(fe_function);
  compute_cells_to_limit(fe_function);
  return is_cell_to_limit;
}

std::vector<std::vector<double>> TSlopeLimiter::get_features(const
    FEFunction* fe_function)
{
  // The features have to be computed in case they haven't been computed yet or
  // with a different fe_function
  compute_features(fe_function);
  return features;
}

void TSlopeLimiter::info()
{
  Output::print<2>("\nSlope limiter of type:    ", limiter_name);
  Output::print<3>("M_lim:                    ", m_lim);
  Output::print<3>("gamma:                    ", gamma);
  Output::print<3>("alpha_ref:                ", alpha_ref);
  Output::print<3>("characteristic_length:    ", characteristic_length);
  Output::print<3>("characteristic_solution_scale:   ", characteristic_solution_scale);
  if (features.empty())
  {
    Output::print<3>("Size of features:         ", "features not computed yet");
  }
  else
  {
    Output::print<3>("Size of features:         ", features.size(), " x ",
        features[0].size());
  }
  if (is_cell_to_limit.empty())
  {
    Output::print<3>("Size of is_cell_to_limit: ", "is_cell_to_limit not ",
        "computed yet");
  }
  else
  {
    Output::print<3>("Size of is_cell_to_limit: ", is_cell_to_limit.size());
  }
}

void TSlopeLimiter::limit_function(FEFunction* const fe_function)
{
  // Limiter consist basically of three steps:
  // 1. For each cell, compute features of the solution and decide based on that
  //    features whether to mark the cell or not
  // 2. Change / limit solution in the marked cells
  //
  // Step 1 can be further divided into two partial steps, namely
  // 1a) Compute the features for each cell
  // 1b) Decide whether to mark the cell or not.
  // The different limiters differ in the details of the steps, but the
  // structure is the same. Therefore, I try to follow this structure and
  // distinguish in between the specific methods.
  //
  // For information about the limiters, see
  // https://doi.org/10.1016/j.cam.2021.113487
  // To understand the code I would suggest that you read and hopefully
  // understand especially Section 3 of this paper.
  // The limiters are:
  // 1. Galerkin (no limiter at all, see also p. 12 in the paper)
  // 2. LinTriaReco (see Section 3.1)
  // 3. ConstTriaReco (see Section 3.1)
  // 4. LinQuadReco (see Section 3.1)
  // 5. ConstQuadReco (see Section 3.1)
  // 6. LinQuadDeriv (see Section 3.2)
  // 7. ConstQuadDeriv (see Section 3.2)
  // 8. ConstJump (see Section 3.3)
  // 9. ConstJumpMod (see Section 3.3)
  // 10. ConstJumpL1Norm (see documentation)
  // 11. ConstJumpL2Norm (see documentation)
  // 12. ConstJumpLinftyNorm (see documentation)

  if (limiter_name == "Galerkin" || !apply_limiter)
  {
    // Galerkin means do nothing.
    return;
  }
#ifdef __3D__
    ErrThrow("The limiter is intended to work in 2D. It's not clear if it works"
        " in 3D.");
#endif

    auto t = GetTime();
    compute_features(fe_function);
    compute_cells_to_limit(fe_function);
    limit_fe_function(fe_function);
    t = GetTime() - t;
    Output::print<3>("time for applying slope limiter: ", t, " seconds");
}

void TSlopeLimiter::compute_features(const FEFunction* const fe_function)
{
#ifdef __3D__
    ErrThrow("The limiter is intended to work in 2D. It's not clear if it works"
        " in 3D.");
#endif
  if (fe_function->GetFESpace() != nullptr &&
      !fe_function->GetFESpace()->is_discontinuous())
  {
    ErrThrow("The limiter are designed to work with discontinuous spaces.");
  }

  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    // Except for LinQuadDeriv and ConstQuadDeriv all methods rely on straight
    // lines, e.g. to get quadrature points, to evaluate in the midpoint of the
    // line or to mirror along a straight line. Therefore, isoparametric
    // elements are forbidden. This can be solved but I don't have the time to
    // do so.
    ErrThrow("This method is only implemented for non-isoparametric ",
        "transformation.");
  }

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

  if (lim_name_mod == "Galerkin")
  {
    // Galerkin means do nothing.
    return;
  }

  // Check if limiter and grid correspond.
  // Since we either have triangular or quadrilateral grids, it is enough to
  // check the first cell.
  if ((lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv") &&
      fe_function->GetFESpace()->GetCollection()->GetCell(0)->GetN_Vertices() !=
      4)
  {
    ErrThrow(lim_name_mod, " does only work on quadrilateral grids.");
  }

  // Starting step 1a)
  // We want to compute cell wise features and store them for later use.
  // Depending on the method we've got a different number of features that are
  // needed. The number is just hard coded here and might make sense if you know
  // the features of the specific limiter. In total we've got
  // (n_cells x n_features) features since we've got features for every cell.
  const auto fe_space = fe_function->GetFESpace();
  const auto coll = fe_space->GetCollection();
  const auto n_cells = coll->GetN_Cells();
  int n_features = 0;
  if (lim_name_mod == "LinTriaReco" || lim_name_mod == "ConstTriaReco")
  {
    // Features are: the integral mean and the evaluation of the solution in the
    // edge midpoints respectively the integral mean along the edges, i.e.
    // 1 + number of edges features per cell. We assume that in the mesh there
    // is only a specific type of cells, i.e. either triangles or quads, and
    // therefore the number of edges per cell is constant.
    auto n_edges = coll->GetCell(0)->GetN_Joints();
    n_features = 1 + n_edges;
    // For the limiter ConstTriaReco we also need the integral mean for virtual
    // cells appearing for boundary edges. Therefore, at most n_edges more
    // feature. The ordering is that after the other features for each edge
    // a virtual cell average can be stored even though it might be unused.
    n_features += (lim_name_mod == "ConstTriaReco") ? n_edges : 0;
  }
  else if (lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv")
  {
    // Features are: first three coefficients of the following decomposition,
    // i.e. 3 features.
    // For a sufficiently large degree, the solution can be locally on each cell
    // K expanded as
    // u_h(x,y) = a0_K + a1_K psi(x) + a2_K xi(y) + higher_order_terms
    // where psi(x,y) equals half of the first component of the inverse of the
    // transformation from the reference element to K, and xi(x,y) equals half
    // of the second component of the inverse of the transformation. The basis
    // functions of the higher order can be chosen almost arbitrarily. BUT they
    // have to satisfy
    // N_0(basis_function) := 1/|K| * \int_K basis_function = 0,
    // N_1(basis_function) := c1 * \int_K^ basis_function(F(x^,y^)) * x^ = 0,
    // N_2(basis_function) := c2 *\int_K^ basis_function(F(x^,y^)) * y^ = 0,
    // where K^ is the reference element, x^ (y^) is the first (second)
    // coordinate on the reference element, and F is the transformation of the
    // reference element to the physical cell. The scaling factors c1 and c2 are
    // determined by N_1(x^)=N_2(y^) = 1 on the reference element, i.e. F = id.
    // The coefficients a0_K, a1_K and a2_K of the decomposition mentioned above
    // are also determined by these nodal functions, i.e.
    // a0_K := NF_0(u_h),
    // a1_K := NF_1(u_h),
    // a2_K := NF_2(u_h).
    // In other words, a1_K and a2_K collect all the x^- and y^-part of the
    // solution and are independent of the basis of the higher order terms. The
    // coefficient a0_K is the cell average of u_h These three coefficients are
    // also the feature in the method.
    n_features = 3;
  }
  else if (lim_name_mod == "ConstJump")
  {
    // Features are: the integral of the jump along the edges, i.e. n_edges
    // features for each cell. As before, the number of edges of a cell is
    // constant in the mesh.
    auto n_edges = coll->GetCell(0)->GetN_Joints();
    n_features = n_edges;
  }
  else if (lim_name_mod == "ConstJumpMod")
  {
    // Features are: the alpha value for each edge. As before, the number of
    // edges of a cell is constant in the mesh.
    auto n_edges = coll->GetCell(0)->GetN_Joints();
    n_features = n_edges;
  }
  else if (lim_name_mod == "ConstJumpL1Norm"
      || lim_name_mod == "ConstJumpL2Norm"
      || lim_name_mod == "ConstJumpLinftyNorm")
  {
    // Features are: the (L2 / L1 / Linfty)-norm of the jump of uh along each
    // edge. Therefore, there are n_edges value in each cell. As before, the
    // number of edges of a cell is constant in the mesh.
    auto n_edges = coll->GetCell(0)->GetN_Joints();
    n_features = n_edges;
  }
  else
  {
    ErrThrow("Unknown limiter! You have to specify a limiter according to ",
        "database[\"lim_name_mod\"], see Documentation.");
  }
  features.resize(n_cells, std::vector<double>(n_features, 0));
  // The following limiter need some more storage, namely to compute the mean
  // and the standard deviation of all norms.
  if (lim_name_mod == "ConstJumpL1Norm" || lim_name_mod == "ConstJumpL2Norm"
      || lim_name_mod == "ConstJumpLinftyNorm")
  {
    features.resize(n_cells + 1, std::vector<double>(n_features, 0));
    // features[n_cells][0] counts the number of edges, [1] contains the mean of
    // integrals, [2] the standard deviation
    features[n_cells].resize(3);
  }

  // To compute the features in LinQuadDeriv and ConstQuadDeriv we need to
  // integrate over the original and reference cell. In ConstTriaReco we need to
  // integrate over reflected cells of boundary cells. To this extend, we
  // prepare a quadrature rule on the reference cell.  The other methods don't
  // need this, so in these cases you can ignore the following lines.
  int max_degree = 0; // maximal needed degree
  if (lim_name_mod == "ConstTriaReco" ||
      lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv")
  {
    // We want to find out the order of the used finite element to have an exact
    // quadrature rule
    auto used_elements = fe_space->GetUsedElements();
    auto n_used_ele = fe_space->GetN_UsedElements();
    for (int ele_i = 0; ele_i < n_used_ele; ++ele_i)
    {
      auto elem_id = used_elements[ele_i];
      FiniteElement fe_temp(elem_id);
      auto degree_temp = fe_temp.GetBaseFunct()->GetPolynomialDegree();
      max_degree = (degree_temp > max_degree) ? degree_temp : max_degree;
    }
  }
  // Prepare quadrature rule: For LinQuadDeriv and ConstQuadDeriv we need to
  // integrate the basis functions times a linear function. Hence, the
  // quadrature rule has to be exact for max_degree + 1. For ConstTriaReco we
  // only need to integrate the FE basis functions themselves.
  auto degree = (lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv") ?
    max_degree + 1 : max_degree;
  const auto joint_type = (coll->GetCell(0)->GetN_Joints() == 4) ?
    BFRefElements::BFUnitSquare : BFRefElements::BFUnitTriangle;
  const auto& qf_ref = *QuadratureFormulaDatabase::qf_from_degree(degree,
                                                                  joint_type);
  const auto xi = qf_ref.get_xi();
  const auto eta = qf_ref.get_eta();
  const auto n_quad_points = qf_ref.GetN_QuadPoints();

  // Loop over cells to fill the features vector.
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    auto cell = coll->GetCell(cell_i);
    // Fill the features vector accordingly to the method
    if ( lim_name_mod == "ConstTriaReco" || lim_name_mod == "LinTriaReco" )
    {
      // Features are the integral mean in the cell and the value of u_h at the
      // edge midpoint respectively the integral mean along the edge. The
      // ordering in features[cell_i] is as listed in the previous sentence,
      // i.e. first the integral mean, then the respective value for the edges.

      // Compute the integral mean
#ifdef __2D__
      features[cell_i][0] = fe_function->compute_cell_average(cell_i);
#endif

      if ( lim_name_mod == "ConstTriaReco" ||
          (lim_name_mod == "LinTriaReco" && !cell->IsBoundaryCell()) )
      {
        // Boundary cells in the method LinTriaReco are not treated at all.
        // Therefore, we don't need features of those cells except the cell
        // average since this is used by neighboring interior cells.

        // Compute the edge values
        auto n_edges = cell->GetN_Edges();
        for (int joint_j = 0; joint_j < n_edges; ++joint_j)
        {
          auto feature = &features[cell_i][1 + joint_j];
          if (lim_name_mod == "LinTriaReco")
          {
            // compute the value at the midpoint of the edge
            auto midp = cell->ComputeMidOfJoint(joint_j);
#ifdef __2D__
            fe_function->FindValueLocal(cell, cell_i, midp.x, midp.y, feature);
#endif
          }
          else
          {
            // compute the integral mean along the edge
            double length;
            double val_integral;
#ifdef __2D__
            fe_function->compute_edge_integral_and_length(cell_i, joint_j,
                val_integral, length);
#endif
            *feature = val_integral / length;
            if (cell->IsBoundaryCell())
            {
              // Here we need a virtual cell and the integral mean of u_h over
              // the virtual cell

              // Compute quadrature points on original cell, i.e. prepare
              // quadrature rule on physical cell
              auto ref_trans = fe_space->get_fe(cell_i).GetRefTransID();
              FEDatabase::SetCellForRefTrans(cell, ref_trans);
              TQuadFormula qf_orig(qf_ref);
              FEDatabase::GetOrigFromRef(ref_trans, qf_ref, qf_orig);

              // Compute the integral and the cell_measure
              double integral = 0;
              double cell_measure = 0;
              for (int pt_i = 0; pt_i < n_quad_points; ++pt_i)
              {
                auto pt = qf_orig.get_point(pt_i);
                double x_reflected;
                double y_reflected;
                cell->ComputeReflectedPoint(joint_j, pt.x, pt.y, x_reflected,
                    y_reflected);
                double val_uh = 0;
#ifdef __2D__
                fe_function->FindValueLocal(cell, cell_i, x_reflected,
                    y_reflected, &val_uh);
#endif
                auto weight_orig = qf_orig.get_weight(pt_i);

                cell_measure += weight_orig;
                integral += weight_orig * val_uh;
              }
              features[cell_i][1 + n_edges + joint_j] = integral / cell_measure;
            }
          }
        } // endfor joint_j
      }
    }
    else if (lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv")
    {
      // Features are the coefficients of the decomposition mentioned above. The
      // coefficients can be computes using the functionals given above.
      //
      // Compute quadrature points on original cell, i.e. prepare quadrature
      // rule on physical cell
      auto ref_trans = fe_space->get_fe(cell_i).GetRefTransID();
      if (ref_trans != ReferenceTransformation_type::QuadAffin)
      {
        // For non-affine transformation the basis functions and the nodal
        // functionals have to be chosen differently. It is not clear yet how to
        // do this. A possible way is to use only the affine part of the
        // transformation. Another way would be to use exactly the same
        // functionals but with non-affine transformation. However, both methods
        // are neither implemented nor tested yet.
        ErrThrow("This method is only implemented for affine transformation.");
      }
      FEDatabase::SetCellForRefTrans(cell, ref_trans);
      TQuadFormula qf_orig(qf_ref);
      FEDatabase::GetOrigFromRef(ref_trans, qf_ref, qf_orig);

      // Evaluate the functionals
      double cell_measure = 0;
      for (int pt_i = 0; pt_i < n_quad_points; ++pt_i)
      {
        double val_uh = 0;
#ifdef __2D__
        auto point = qf_orig.get_point(pt_i);
        fe_function->FindValueLocal(cell, cell_i, point.x, point.y, &val_uh);
#endif
        auto weight_ref = qf_ref.get_weight(pt_i);
        auto weight_orig = qf_orig.get_weight(pt_i);

        cell_measure += weight_orig;
        features[cell_i][0] += weight_orig * val_uh; // N_0(u_h)
        features[cell_i][1] += weight_ref * val_uh * xi[pt_i]; // N_1(u_h)
        features[cell_i][2] += weight_ref * val_uh * eta[pt_i]; // N_2(u_h)
      }
      // correct scaling of the functionals for Kronecker delta property of the
      // basis functions. The numbers come from N_1(x^)=N2(y^) = 1 on the
      // reference cell.
      features[cell_i][0] *= 1 / cell_measure;
      features[cell_i][1] *= 2.0 / 3;
      features[cell_i][2] *= 2.0 / 3;
    }
    else if (lim_name_mod == "ConstJump" || lim_name_mod == "ConstJumpMod"
        || lim_name_mod == "ConstJumpL1Norm"
        || lim_name_mod == "ConstJumpL2Norm"
        || lim_name_mod == "ConstJumpLinftyNorm")
    {
      // Compute the jump integrals and possibly edge alphas
      auto n_edges = cell->GetN_Joints();

      // Construct a 1D quadrature rule to integrate along edges
      auto degree = fe_space->getFEDegree(cell);
      for (int edge_j = 0; edge_j < n_edges; ++edge_j)
      {
        // Check for neighbours to distinguish boundary and inner edges
        auto joint = cell->GetJoint(edge_j);
        auto neigh = joint->GetNeighbour(cell);
        int cell_nr_neigh = 0;
        // These limiter are only defined for interior edges, since it's not
        // clear how the jump at the boundary is defined and how it can be used
        // as limiting criterion -> Exclude boundary edges from computations
        if (!neigh)
        {
          // To exclude ConstJumpMod the feature has to be at least alpha_ref,
          // see also compute_cells_to_limit
          if (lim_name_mod == "ConstJumpMod")
          {
            features[cell_i][edge_j] = alpha_ref;
          }
          continue;
        }
        cell_nr_neigh = neigh->GetCellIndex(); // cell number of neighbor
        if (cell_i > cell_nr_neigh)
        {
          // We consider the edges only once since the integral of the jump
          // along an edge is independent of the cell where this particular
          // limiter belongs to.
          continue;
        }

        // needed degree is the maximum of the degrees of the fe
        // functions in the two adjacent cells
        auto degree_neigh = fe_space->getFEDegree(neigh);
        degree = (degree < degree_neigh) ? degree_neigh : degree;

        // Prepare quadrature rule for integration along edges
        // We either integrate FEFunction or FEFunction * FEFunction, hence the
        // accuracy has to be degree or 2 * degree
        degree = (lim_name_mod == "ConstJumpL1Norm") ? degree : 2 * degree;
        auto quad_form_1D = QuadratureFormulaDatabase::qf_from_degree(
            2 * degree, BFRefElements::BFUnitLine);
        auto n_quad_pts = quad_form_1D->GetN_QuadPoints();
        auto fe = fe_space->get_fe(cell_i);
        auto ref_element = fe.GetBaseFunct()->GetRefElement();
        std::vector<std::vector<double>> quad_points_ref(2,
            std::vector<double>(n_quad_pts));
        std::vector<std::vector<double>> quad_points_orig(2,
            std::vector<double>(n_quad_pts));

        for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
        {
          auto zeta_1D = quad_form_1D->get_point(quad_pt_i).x;
          auto xi_eta_2D = transform(ref_element, edge_j, zeta_1D);
          quad_points_ref[0][quad_pt_i] = xi_eta_2D.x;
          quad_points_ref[1][quad_pt_i] = xi_eta_2D.y;
        } // endfor quad_pt_i
        FEDatabase::GetRefTrans2D(fe.GetRefTransID())->SetCell(cell);
        // compute points on original cell
        FEDatabase::GetOrigFromRef(fe.GetRefTransID(), n_quad_pts,
            quad_points_ref[0].data(), quad_points_ref[1].data(),
            quad_points_orig[0].data(), quad_points_orig[1].data());

        // compute the integral of the square of the jump
        double integral = 0;
        // For neighboring elements the jump is given by the difference of the
        // function values at the quadrature points. Since we compute the square
        // of the jump the "direction" of the jump does not matter.  For
        // boundary elements the jump is here defined to be 0 such that these
        // cells don't contribute to the marking criterion.
        for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
        {
          double jump = 0;
          double val_cell = 0;
          double val_neigh = 0;
#ifdef __2D__
          fe_function->FindValueLocal(cell, cell_i,
              quad_points_orig[0][quad_pt_i], quad_points_orig[1][quad_pt_i],
              &val_cell);
          fe_function->FindValueLocal(neigh, cell_nr_neigh,
              quad_points_orig[0][quad_pt_i], quad_points_orig[1][quad_pt_i],
              &val_neigh);
#endif
          jump = val_cell - val_neigh;
          auto weight_orig = quad_form_1D->get_weight(quad_pt_i);
          if (lim_name_mod == "ConstJumpL1Norm")
          {
            integral += weight_orig * std::abs(jump); // L1-norm
          }
          else if (lim_name_mod == "ConstJumpLinftyNorm")
          {
            auto abs = std::abs(jump);
            integral = (abs > integral) ? abs : integral; // Linfty-norm
          }
          else // ConstJump, ConstJumpMod or ConstJumpL2Norm
          {
            integral += weight_orig * jump * jump; // L2-norm (squared)
          }
        }

        // For ConstJumpLinftyNorm include the values at the vertices
        if (lim_name_mod == "ConstJumpLinftyNorm")
        {
          auto v1 = cell->GetVertex(edge_j);
          double val_cell = 0;
          double val_neigh = 0;
#ifdef __2D__
          fe_function->FindValueLocal(cell, cell_i, v1->GetX(), v1->GetY(),
              &val_cell);
          fe_function->FindValueLocal(neigh, cell_nr_neigh, v1->GetX(),
              v1->GetY(), &val_neigh);
#endif
          auto absjump = std::abs(val_cell - val_neigh);
          integral = (absjump > integral) ? absjump : integral; // Linfty-norm
          auto v2 = cell->GetVertex((edge_j + 1 ) % n_edges);
#ifdef __2D__
          fe_function->FindValueLocal(cell, cell_i, v2->GetX(), v2->GetY(),
              &val_cell);
          fe_function->FindValueLocal(neigh, cell_nr_neigh, v2->GetX(),
              v2->GetY(), &val_neigh);
#endif
          absjump = std::abs(val_cell - val_neigh);
          integral = (absjump > integral) ? absjump : integral; // Linfty-norm
        }

        // The above calculated weight has to be transformed to be original
        // cell which was not done yet. The transformation is given by half
        // of the length of the joint. Therefore, adjust the integral.
        double hE = 0;
        if (lim_name_mod != "ConstJumpLinftyNorm")
        {
          hE = cell->ComputeDiameterOfJoint(edge_j);
          integral *= hE / 2;
        }
        if (lim_name_mod == "ConstJumpL1Norm" || lim_name_mod == "ConstJumpL2Norm")
        {
          // Divide by length of edge according to formula
          integral /= hE;
          // For L2 norm compute the square root
          integral = (lim_name_mod == "ConstJumpL2Norm") ? std::sqrt(integral)
            : integral;
        }

        // save computed integral in this cell
        features[cell_i][edge_j] = integral;

        if (lim_name_mod == "ConstJumpMod")
        {
          // For this limiter the edge alpha has to be computed and can be
          // stored in the features vector.
          auto scaled_length = hE / characteristic_length;
          if (scaled_length >= 1)
          {
            ErrThrow("Choose characteristic_length larger than longest edge in",
                " the triangulation.  hE: ", hE, "  characteristic_length: ",
                characteristic_length);
          }
          double C = C0_CJM * characteristic_solution_scale
            * characteristic_solution_scale * characteristic_length;
          auto alpha = std::log(1. / C * integral) / std::log(scaled_length);
          features[cell_i][edge_j] = alpha;
        }
        // Save feature also in the neighboring cell
        // calculate local edge number in neighbor cell
        auto edge_nr_neigh = joint->get_joint_nr_in_cell(neigh);
        features[cell_nr_neigh][edge_nr_neigh] = features[cell_i][edge_j];

        if (lim_name_mod == "ConstJumpL1Norm"
            || lim_name_mod == "ConstJumpL2Norm"
            || lim_name_mod == "ConstJumpLinftyNorm")
        {
          features[n_cells][0] += 1;  // count number of interior edges
          features[n_cells][1] += features[cell_i][edge_j]; // sum up integrals
        }
      } // endfor edge_j
    }
  } //endfor cell_i
  // Compute integral mean for some limiter
  if (lim_name_mod == "ConstJumpL1Norm"
      || lim_name_mod == "ConstJumpL2Norm"
      || lim_name_mod == "ConstJumpLinftyNorm")
  {
    // compute mean of integrals
    features[n_cells][1] /= features[n_cells][0];

    // compute standard deviation and store it in features[n_cells][2]
    for (int cell_i = 0; cell_i < n_cells; ++cell_i)
    {
      auto cell = coll->GetCell(cell_i);
      // Compute the jump integrals and possibly edge alphas
      auto n_edges = cell->GetN_Joints();

      for (int edge_j = 0; edge_j < n_edges; ++edge_j)
      {
        // Check for neighbours to distinguish boundary and inner edges
        auto joint = cell->GetJoint(edge_j);
        auto neigh = joint->GetNeighbour(cell);
        // These limiter are only defined for interior edges, since it's not
        // clear how the jump at the boundary is defined and how it can be used
        // as limiting criterion -> Exclude boundary edges from computations
        if (!neigh)
        {
          continue;
        }
        int cell_nr_neigh = neigh->GetCellIndex(); // cell number of neighbor
        if (cell_i > cell_nr_neigh)
        {
          // Consider the edges only once
          continue;
        }

        auto summand = features[cell_i][edge_j] - features[n_cells][1];
        summand *= summand;
        features[n_cells][2] += summand;
      } // endfor edge_j
      features[n_cells][2] /= features[n_cells][0];
      features[n_cells][2] = std::sqrt(features[n_cells][2]);
    }
  }
}

void TSlopeLimiter::compute_cells_to_limit(const FEFunction* const fe_function)
{
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

  if (lim_name_mod == "Galerkin")
  {
    // Galerkin means do nothing.
    return;
  }

  // Starting step 1b)
  // After having computed the features needed, we can mark cells according
  // to the respective limiter.
  const auto fe_space = fe_function->GetFESpace();
  const auto coll = fe_space->GetCollection();
  const auto n_cells = coll->GetN_Cells();
  is_cell_to_limit.resize(n_cells, false);



  // Loop over cells to fill the array
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    auto cell = coll->GetCell(cell_i);
    if ( lim_name_mod == "ConstTriaReco" ||
        (lim_name_mod == "LinTriaReco" && !cell->IsBoundaryCell()) )
    {
      auto n_edges = cell->GetN_Joints();
      for (int joint_j = 0; joint_j < n_edges; ++joint_j)
      {
        auto neigh = cell->GetJoint(joint_j)->GetNeighbour(cell);
        if (neigh)
        {
          auto cell_nr_neigh = neigh->GetClipBoard(); // neighbor
          if (cell_nr_neigh < 0)
          { // mesh without conforming closure is not implemented here yet
            ErrThrow("Meshes with hanging nodes are not allowed in this ",
                "method.");
          }

          // check whether the feature lies between the value at the centroid of
          // the cell and the centroid of the neighbor
          double tol = 1e-11; // tolerance for rounding errors
          auto low = std::min(features[cell_i][0],
              features[cell_nr_neigh][0]);
          auto high = std::max(features[cell_i][0],
              features[cell_nr_neigh][0]);
          if ( features[cell_i][1+joint_j] < low-tol ||
              features[cell_i][1+joint_j] > high+tol)
          { // limiting has to be done
            is_cell_to_limit[cell_i] = true;
            break;
          } // endif limiting has to be done
        }
        else
        {
          // boundary edge for ConstTriaReco. Note, that for LinTriaReco this case
          // cannot happen due to the if a few lines earlier

          // check whether the feature lies between the value at the centroid of
          // the cell and the centroid of the virtual neighbor
          double tol = 1e-11; // tolerance for rounding errors
          auto low = std::min(features[cell_i][0],
              features[cell_i][1 + n_edges + joint_j]);
          auto high = std::max(features[cell_i][0],
              features[cell_i][1 + n_edges + joint_j]);
          if ( features[cell_i][1+joint_j] < low-tol ||
              features[cell_i][1+joint_j] > high+tol)
          { // limiting has to be done
            is_cell_to_limit[cell_i] = true;
            break;
          } // endif limiting has to be done
        }
      }
    }
    else if (lim_name_mod == "LinQuadDeriv" || lim_name_mod == "ConstQuadDeriv")
    {
      // It is checked whether the linear coefficients a1_K and a2_K are "too
      // large" in the sense of the paper mentioned in the beginning of this
      // method.

      double a1_K_bar;
      if (std::abs(features[cell_i][1]) <= m_lim)
      {
        a1_K_bar = features[cell_i][1];
      }
      else
      {
        // Now minmod has to be called with three neighbors. In the
        // case of a cell having (possibly several) boundary edges the cell has
        // less neighbors than edges. To be precisely, the numbers of neighbors
        // equals the numbers of edges minus the number of boundary edges. In
        // the case of boundary cells, minmod is called with an dummy entry that
        // does not effect the evaluation of minmod. As dummy entry the
        // coefficient a1_K is chosen since this value is part of the evaluation
        // of minmod and therefore does not effect the result of minmod.
        std::vector<double> minmod_vec_x(3, features[cell_i][1]);

        auto n_edges = cell->GetN_Edges();
        auto neigh_right = cell->GetJoint(1 % n_edges)->GetNeighbour(cell);
        if (neigh_right)
        {
          auto neigh_nr_right = neigh_right->GetCellIndex();
          minmod_vec_x[1] = gamma *
            (features[neigh_nr_right][0] - features[cell_i][0]);
        }

        auto neigh_left = cell->GetJoint(3 % n_edges)->GetNeighbour(cell);
        if (neigh_left)
        {
          auto neigh_nr_left = neigh_left->GetCellIndex();
          minmod_vec_x[2] = gamma *
            (features[cell_i][0] - features[neigh_nr_left][0]);
        }
        a1_K_bar = utilities::minmod(minmod_vec_x);
      }

      double a2_K_bar;
      if (std::abs(features[cell_i][2]) <= m_lim)
      {
        a2_K_bar = features[cell_i][2];
      }
      else
      {
        // Call minmod as before. Here the dummy value is a2_K.
        std::vector<double> minmod_vec_y(3, features[cell_i][2]);

        auto n_edges = cell->GetN_Edges();
        auto neigh_top = cell->GetJoint(2 % n_edges)->GetNeighbour(cell);
        if (neigh_top)
        {
          auto neigh_nr_top = neigh_top->GetCellIndex();
          minmod_vec_y[1] = gamma *
            (features[neigh_nr_top][0] - features[cell_i][0]);
        }

        auto neigh_bottom = cell->GetJoint(0 % n_edges)->GetNeighbour(cell);
        if (neigh_bottom)
        {
          auto neigh_nr_bot = neigh_bottom->GetCellIndex();
          minmod_vec_y[2] = gamma *
            (features[cell_i][0] - features[neigh_nr_bot][0]);
        }
        a2_K_bar = utilities::minmod(minmod_vec_y);
      }
      auto tol = 1e-12;
      if (std::abs(a1_K_bar - features[cell_i][1]) > tol ||
          std::abs(a2_K_bar - features[cell_i][2]) > tol )
        // Actually a_1_K_bar == expansion_coeffs[cell_i][1] is the criterion
        // that has to be checked. But to avoid the effect of round-off errors
        // we check if a_1_K_bar lies within tol around
        // expansion_coeffs[cell_i][1]. Analogously for a2_K_bar
      {
        is_cell_to_limit[cell_i] = true;
        if (lim_name_mod == "LinQuadDeriv")
        {
          features[cell_i][1] = a1_K_bar;
          features[cell_i][2] = a2_K_bar;
        }
      }

    }
    else if (lim_name_mod == "ConstJump")
    {
      double indicator = 0;
      auto n_edges = cell->GetN_Joints();
      for (int edge_j = 0; edge_j < n_edges; ++edge_j)
      {
        indicator += features[cell_i][edge_j];
      }
      auto measure = cell->GetMeasure();
      auto hK = cell->GetDiameter();
      indicator /= (std::pow(measure, 0.75) * hK);
      if (indicator > 1)
      {
        is_cell_to_limit[cell_i] = true;
      }
    }
    else if ( lim_name_mod == "ConstJumpMod" )
    {
      double alpha_min = features[cell_i][0];
      for (auto alpha : features[cell_i])
      {
        alpha_min = (alpha < alpha_min) ? alpha : alpha_min;
      }
      if (alpha_min < alpha_ref)
      {
        is_cell_to_limit[cell_i] = true;
      }
    }
    else if (lim_name_mod == "ConstJumpL1Norm"
        || lim_name_mod == "ConstJumpL2Norm"
        || lim_name_mod == "ConstJumpLinftyNorm")
    {
      // mark cells where at least one edge has a larger norm than the mean of
      // the norm
      auto mean = features[n_cells][1];
      auto standard_deviation = features[n_cells][2];
      auto ref_val = mean + C1_CJN * standard_deviation;
      auto n_edges = cell->GetN_Joints();
      // To exclude continuous solutions the mean should be larger than a given
      // threshold. Here 1e-13 is chosen but the choice is arbitrary
      if (ref_val > 1e-13)
      {
        if (cell->GetCellIndex() == 0)
        {
          Output::info("TSlopeLimiter", "mean:  ", mean, "\t\t",
              "standard deviation:  ", standard_deviation, " (",
              standard_deviation/mean*100, " %)\t\t", "ref_val:  ", ref_val);
        }
        for (int edge_j = 0; edge_j < n_edges; ++edge_j)
        {
          if (features[cell_i][edge_j] > ref_val)
          {
            is_cell_to_limit[cell_i] = true;
            break;
          }
        }
      }
    } // end lim_name_mod
  }
}

void TSlopeLimiter::limit_fe_function(FEFunction* const fe_function)
{
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

  if (lim_name_mod == "Galerkin")
  {
    // Galerkin means do nothing.
    return;
  }

  // Starting step 2: Limit the solution in the marked cells
  // This step can be decomposed into computing a new function, and replacing
  // the solution by an interpolation of the previous computed function.
  // This new function strongly depends on the used method. The interpolation
  // and the replacement is then for all methods equal. Both steps are conducted
  // locally on each cell.

  // Get access to the FE functions values to replace entries later
  auto values = fe_function->GetValues();

  // All methods have in common that the new function is at most linear and
  // can therefore in 2D has 3 degrees of freedom, i.e. a constant value and
  // a linear values in x resp. y direction. These coefficients are stored
  // in the following array. Since it's done only locally, storage for only
  // three coefficients is needed and can be reused on each cell
  std::array<double, 3> coeffs = {0, 0, 0};

  const auto fe_space = fe_function->GetFESpace();
  const auto coll = fe_space->GetCollection();
  const auto n_cells = coll->GetN_Cells();
  is_cell_to_limit.resize(n_cells);
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    if (!is_cell_to_limit[cell_i])
    {
      // Nothing to be done
      continue;
    }

    auto cell = coll->GetCell(cell_i);
    // Fill the coeffs array according to the specific limiter
    if (lim_name_mod == "LinTriaReco" && !cell->IsBoundaryCell())
    {
      // Still boundary cells are not touched by this algorithm

      // Please consult the paper mentioned in the beginning of the method for
      // details about the limiting procedure. I try to explain why and how
      // things are happening but without explaining the method itself.

      // We need the number of edges since for each interior edge a limiter is
      // constructed
      auto n_edges = cell->GetN_Joints();

      // For this method we need the centroids of the cell and its neighbors.
      // These centroids consist of 2 coordinates.
      std::vector<parmoon::Point> centroids_neigh(n_edges,
          parmoon::Point((unsigned int) 2));
      auto centroid_cell = cell->getCentroid();
      // We further also need the cell numbers of the neighbors
      std::vector<int> cell_nr_neighs(n_edges);
      for (int joint_j = 0; joint_j < n_edges; ++joint_j)
      {
        auto neigh = cell->GetJoint(joint_j)->GetNeighbour(cell);
        centroids_neigh[joint_j] = neigh->getCentroid();
        cell_nr_neighs[joint_j] = neigh->GetCellIndex();
      }

      // Construct the limiter candidates for the interior edges, i.e.
      // a linear function that has the cell average as value at the centroid
      // of the cell and the cell average of two neighboring cells as value at
      // the centroids of those cells. These candidates again have three
      // coefficients since they are at most linear functions
      std::vector<std::array<double, 3>>coeffs_candidates(n_edges);
      for (int joint_j = 0; joint_j < n_edges; ++joint_j)
      {
        // The computation of the candidate depends on two edges and the
        // respective adjacent cells. Therefore we fix the corresponding joint
        // number and cell number.
        auto joint_nr_next = (joint_j + 1) % n_edges;
        auto joint_nr_2next = (joint_j + 2) % n_edges;
        auto cell_nr_next = cell_nr_neighs[joint_nr_next];
        auto cell_nr_2next = cell_nr_neighs[joint_nr_2next];

        // A linear system has to be solved to compute the coefficients of the
        // candidates. This is done here, where the linear system is solved by
        // hand.
        // For each candidate the problem can be written as
        // A_ja_j = b_j
        // where
        // A_j = {  1, centroid_cell
        //          1, centroids_neigh(joint_nr_next)
        //          1, centroids_neigh(joint_nr_2next)  }
        // a_j = {  a_j,0
        //          a_j,1
        //          a_j,2 }
        // b_j = {  features[cell_nr][0]
        //          features[cell_nr_next][0]
        //          features[cell_nr_2next][0] }
        // The notation used here is adapted from
        // https://en.wikipedia.org/wiki/Inverse_matrix#Inversion_of_3_%C3%97_3_matrices
        // Since the template parameter is called d, the "d" from the wikipedia
        // page is here denoted by D
        double a = 1;
        double b = centroid_cell.x;
        double c = centroid_cell.y;
        double D = 1;
        double e = centroids_neigh[joint_nr_next].x;
        double f = centroids_neigh[joint_nr_next].y;
        double g = 1;
        double h = centroids_neigh[joint_nr_2next].x;
        double i = centroids_neigh[joint_nr_2next].y;
        double inv_A_j [3][3] = {
          {e*i - f*h, c*h - b*i, b*f - c*e},
          {f*g - D*i, a*i - c*g, c*D - a*f},
          {D*h - e*g, b*g - a*h, a*e - b*D}
        };  // this is the inverse matrix x determinant
        auto det = a * inv_A_j[0][0] + b * inv_A_j[1][0] + c * inv_A_j[2][0];
        double b_j [3] = {features[cell_i][0], features[cell_nr_next][0],
          features[cell_nr_2next][0]};
        for (int i = 0; i < 3; ++i)
        {
          coeffs_candidates[joint_j][i] = 0;
          for (int k = 0; k < 3; ++k)
          {
            coeffs_candidates[joint_j][i] += inv_A_j[i][k] * b_j[k];
          }
          coeffs_candidates[joint_j][i] /= det; // scale to correct value
        }
      }

      // Now we've got candidates at hand that have to be ordered by
      // decreasing Euclidean norm of their linear coefficients. This is
      // stored in rating_limiter
      std::vector<double> rating_limiter(n_edges);
      for (int joint_j = 0; joint_j < n_edges; ++joint_j)
      {
        rating_limiter[joint_j] = std::sqrt(
            coeffs_candidates[joint_j][1]*coeffs_candidates[joint_j][1] +
            coeffs_candidates[joint_j][2]*coeffs_candidates[joint_j][2] );
      }

      // ordering is a container in which the numbers
      // 1, 2, ..., n_interior_edges are ordered such that
      // rating_limiter[ordering[i]] > rating_limiter[ordering[j]], for i < j.
      // It is an index vector pointing in rating_limiter such that the
      // values of rating_limter are ordered in descending order
      std::vector<int> ordering(n_edges);
      int iterator = 0;
      std::iota(ordering.begin(), ordering.end(), iterator++);
      std::sort( ordering.begin(), ordering.end(),
          [&](int i, int j) {return rating_limiter[i] > rating_limiter[j];}
          );

      // Using this ordering, we can check the candidates if their evaluation
      // at the edge midpoints is between the cell averages of the solution in
      // the adjacent cells.
      // The number optimal candidate is stored in the following variable. We
      // initialize it with a value that is impossible to later verify the
      // programs choice
      auto optimal_candidate = n_edges;

      for (auto candidate_nr : ordering)
      {
        // Check if the candidate satisfies the marking criterion, i.e. is the
        // value at the edge midpoint between the cell averages of u_h in the
        // adjacent cells.
        bool has_passed = true;
        for (int joint_j = 0; joint_j < n_edges; ++joint_j)
        {
          if ( joint_j == (candidate_nr + 1) % n_edges ||
              joint_j == (candidate_nr + 2) % n_edges )
          {
            // For these joints the candidate satisfy the criterion by
            // construction
            continue;
          }

          // compute the value at the midpoint of the edge
          auto p = cell->ComputeMidOfJoint(joint_j);
          auto val_midpoint = coeffs_candidates[candidate_nr][0]
            + coeffs_candidates[candidate_nr][1] * p.x
            + coeffs_candidates[candidate_nr][2] * p.y;

          // check whether value at midpoint lies between the value at the
          // centroid of the cell and the centroid of the neighbor
          double tol = 1e-11; // tolerance for rounding errors
          auto neigh = cell->GetJoint(joint_j)->GetNeighbour(cell);
          auto cell_nr_neigh = neigh->GetCellIndex();
          auto low = std::min(features[cell_i][0],
              features[cell_nr_neigh][0]);
          auto high = std::max(features[cell_i][0],
              features[cell_nr_neigh][0]);
          if ( val_midpoint < low-tol ||
              val_midpoint > high+tol)
          {
            // discard candidate since it does not satisfy the criterion
            has_passed = false;
            break;
          }
        }
        if (has_passed)
        {
          // We found an acceptable limiter and can stop further searching
          optimal_candidate = candidate_nr;
          break;
        }
        else if (candidate_nr == ordering.back())
        {
          // none of the candidates is appropriate -> use integral mean as
          // reconstruction. We set therefore the optimal candidate to some
          // number, say -1
          optimal_candidate = -1;
        }
      }
      if ( optimal_candidate < -1 || optimal_candidate >= n_edges)
      {
        ErrThrow("Something went wrong with choosing a appropriate function ",
            "for reconstruction.");
      }
      else if (optimal_candidate == -1)
      {
        // Replace by integral mean
        coeffs[0] = features[cell_i][0];
        coeffs[1] = coeffs[2] = 0;
      }
      else
      {
        // Replace by chosen limiter
        for (auto i = 0; i < 3; ++i)
        {
          coeffs[i] = coeffs_candidates[optimal_candidate][i];
        }
      }
    }
    else if (lim_name_mod == "LinQuadDeriv")
    {
      for (auto i = 0; i < 3; ++i)
      {
        coeffs[i] = features[cell_i][i];
      }
    }
    else if (lim_name_mod == "ConstTriaReco" || lim_name_mod == "ConstQuadDeriv")
    {
      // Replace by integral mean, i.e. change the constant coefficient but
      // leave the remaining coefficients untouched to be 0.
      coeffs[0] = features[cell_i][0];
    }
    else if (lim_name_mod == "ConstJump" || lim_name_mod == "ConstJumpMod"
        || lim_name_mod == "ConstJumpL1Norm" || lim_name_mod ==
        "ConstJumpL2Norm" || lim_name_mod == "ConstJumpLinftyNorm")
    {
      // Replace by integral mean, i.e. change the constant coefficient but
      // leave the remaining coefficients untouched to be 0.
#ifdef __2D__
      coeffs[0] = fe_function->compute_cell_average(cell_i);
#endif
    }

    // Replace the solution
    // Now we are ready to change the FE Function to the newly constructed
    // function or the cell average.
    // For this, we have to evaluate the nodal functionals of the FE at the
    // new function to change values accordingly

    // Find out the points of the nodal functions where to evaluate the
    // limiter
    auto fe = fe_space->get_fe(cell_i);
    auto nodal_functionals = fe.GetNodalFunctional();
    int n_points_NF;
    const double* xi_nf; // points on reference cell
    const double* eta_nf;  // points on reference cell
    nodal_functionals->GetPointsForAll(n_points_NF, xi_nf, eta_nf);

    // Get points where to evaluate the limiter
    std::vector<double> x(n_points_NF);  // points on original cell
    std::vector<double> y(n_points_NF);  // points on original cell
    if (lim_name_mod == "LinTriaReco")
    {
      auto ref_trans = fe_space->get_fe(cell_i).GetRefTransID();
      FEDatabase::SetCellForRefTrans(cell, ref_trans);
      FEDatabase::GetOrigFromRef(ref_trans, n_points_NF, xi_nf, eta_nf,
          x.data(), y.data());
    }

    // Evaluate the limiter at the points
    std::vector<double> val_limiter(n_points_NF, coeffs[0]);
    if (lim_name_mod == "LinTriaReco")
    {
      for (int pt_i = 0; pt_i < n_points_NF; ++pt_i)
      {
        val_limiter[pt_i] = coeffs[0] + coeffs[1] * x[pt_i]
          + coeffs[2] * y[pt_i];
      } // endfor
    }
    else if (lim_name_mod == "LinQuadDeriv")
    {
      // Evaluate the limiter at the points of the nodal functionals
      for (int pt_i = 0; pt_i < n_points_NF; ++pt_i)
      {
        // These are exactly the basis functions in the decomposition given
        // above, i.e. the basis functions that satisfy the Kronecker delta
        // property of the functionals.
        double val_psi = xi_nf[pt_i] / 2;
        double val_xi = eta_nf[pt_i] / 2;
        val_limiter[pt_i] = coeffs[0] + coeffs[1] * val_psi
          + coeffs[2] * val_xi;
      } // endfor
    }

    // Compute the new coefficients
    std::vector<double> new_coeffs(fe.GetN_DOF());
    nodal_functionals->GetAllFunctionals(nullptr, nullptr, val_limiter.data(),
        new_coeffs.data());

    // Replace the coefficients of this FE function for this cell by the new
    // coefficients
    for (unsigned int dof_i = 0; dof_i < fe_space->get_n_local_dof(cell_i);
        ++dof_i)
    {
      values[fe_space->get_global_dof(cell_i, dof_i)] = new_coeffs[dof_i];
    }
  }
}
