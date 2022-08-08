// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================
#define __HMM_1986__

double DIFFUSION;

void ExampleFile()
{
  Output::print("Example: HMM1986.h with diffusion set to ", DIFFUSION);
}
// exact solution (this is the solution for eps = 0)
auto& Exact = unknown_solution_2d;

// kind of boundary condition
void BoundCondition(int, double, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
bool modified_HMM;

void BoundValue(int BdComp, double Param, double &value)
{
  double discontinuity_point = (modified_HMM) ? 0.25 : 0.3;
    switch(BdComp)
    {
	case 0:
	case 1:
	    value = 0;
	    break;
	case 2:
	    if (Param < 1e-6)
		value = 0;
	    else
		value = 1;
	    break;
	case 3:
	    if ( Param < discontinuity_point )
		value = 1;
	    else
		value = 0;
	    break;
       	    
    }
}

void BilinearCoeffs(int n_points, const double*, const double*,
                    const double*const*, double **coeffs)
{
  static double eps=DIFFUSION;
  int i;
  double *coeff, arg;

  arg = -M_PI/3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = std::cos(arg);
    coeff[2] = std::sin(arg);
    coeff[3] = 0;
    coeff[4] = 0;
  }
}

double compute_layer_width(const TFEFunction2D& u, double y_cut)
{
  // note that this is heavily inspired by 'ComputeBdLayer' in Hemker1996.h
  // Consider a line starting at (xcut, y0) going in the positive y-direction.
  // This should be behind the cylinder where the solution is one and then 
  // (going up) zero.
  // We will go along that line to find the first point where ufct is below 
  // valmax, then we go further along that line to find the first point where 
  // ufct is below valmin. Remembering where these two transitions appeared, we
  // are able to compute the "layer width". This is described also in 
  // Augustin, Caiazzo, Fiebach, Fuhrmann, John, Linke, Umla: "An assessment of 
  // discretizations for convection-dominated convection–diffusion equations",
  // 2011.
  // We avoid having to search all cells for every point, which would be very
  // expensive, by only searching the neighboring cells from where the previous
  // point was found.
  constexpr int bound_points = 100001;
  constexpr double valmin = 0.1;
  constexpr double valmax = 0.9;
  constexpr double x0 = 0.;
  constexpr double Width = 1;

  double xstart = x0;
  double xend = Width;
  double h = Width / (bound_points-1);
  parmoon::Point p(xstart, y_cut);

  // First, we have to determine the cell in which the starting point (0, y_cut)
  // lies. Therefore, we have to do a loop over all cells and determine the
  // cell. This cell is called starting_cell_nr.
  auto collection = u.GetFESpace()->GetCollection();
  auto n_cells = collection->GetN_Cells();
  int starting_cell_nr = 0;
  for(int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    auto current_cell = collection->GetCell(cell_i);
    bool is_point_in_cell = current_cell->PointInCell(p);
    if (is_point_in_cell)
    {
      starting_cell_nr = cell_i;
      break;
    }
  }
  
  // Now a loop over all x values can be performed. In each iteration we have to
  // determine the cell in which the current point lies. For this purpose the
  // current cell number and the current cell are determined. In the beginning
  // the current cell number is of course the previously computed
  // starting_cell_nr.
  int current_cell_nr = starting_cell_nr;
  int jstart = 0;
  for(int j = 0; j < bound_points; j++)
  {
    p.x = x0 + h * j;
    collection->find_cell(p, current_cell_nr);
    double val;
    auto current_cell = collection->GetCell(current_cell_nr);
    u.FindValueLocal(current_cell, current_cell_nr, p.x, p.y, &val);
    if(val > valmin)
    {
      xstart = p.x;
      jstart = j;
      break;
    }
  }

  for(int j=jstart+1;j<bound_points;j++)
  {
    p.x = x0 + h * j;
    collection->find_cell(p, current_cell_nr);
    double val;
    auto current_cell = collection->GetCell(current_cell_nr);
    u.FindValueLocal(current_cell, current_cell_nr, p.x, p.y, &val);
    if(val > valmax)
    {
      xend = p.x;
      break;
    }
  }
  auto lw = std::abs(xend - xstart);

  if (lw >= std::abs(Width - xstart))
  {
    Output::warn("HMM1986::compute_layer_width", "Failed to compute layer ",
        "width at y=",y_cut,". I did not find an end point.");
  }
  if (lw <= h)
  {
    Output::warn("HMM1986::compute_layer_width", "Layer width equals step ",
        "width at y=",y_cut,"!");
  }
  return lw;
}

/*******************************************************************************/
// computes the data for the evaluation of the solution 
/*******************************************************************************/
void hmm1986_postprocessing(ConvectionDiffusion<2> & cd2d)
{
  auto & uh = cd2d.get_function();

  // Compute cutlines
  std::array<double, 3> y_cut{0.1, 0.25, 0.5}; // add more as desired
  for(auto y : y_cut)
  {
    auto lw = compute_layer_width(uh, y);
    Output::print("layer width at y=", y, ": ", std::setprecision(12), lw);
  }

  // Compute (an approximation of) the global minimal and maximal value of uh
  uh.PrintMinMax("solution", false);
  // Compute the global minimal and maximal value of the FE-vector that is
  // cell wise transformed to Pk elements of the same order (if not Pk elements
  // are used already)
  uh.PrintMinMax("Pk solution vector", true);

  // Compute mean oscillation using cell wise minimum and maximum of uh
  constexpr double umin = 0.;
  constexpr double umax = 1.;
  double osc_mean = uh.compute_mean_oscillation(umin, umax, false);
  Output::print("mean oscillation ", std::setprecision(12), osc_mean);
  // Compute mean oscillation using cell wise minimum and maximum of the
  // FE-vector transformed to Pk elements of the same order
  osc_mean = uh.compute_mean_oscillation(umin, umax, true);
  Output::print("mean oscillation using Pk nodal functionals ",
      std::setprecision(12), osc_mean);
  // Compute mean oscillations of the FE-vector itself
  if (!cd2d.get_space()->is_discontinuous())
  {
    auto& entries = cd2d.get_solution().get_entries_vector();
    double osc_mean_vec = 0.;
    for(double e : entries)
    {
      osc_mean_vec += std::max(0., e-umax) + std::max(0., umin-e);
    }
    osc_mean_vec /= entries.size();
    Output::print("mean vector oscillation ", std::setprecision(12),
        osc_mean_vec);
  }
}
