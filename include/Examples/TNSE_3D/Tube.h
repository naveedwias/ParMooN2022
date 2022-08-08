/*********************************************************
 * Setup for simulation of flow through ascending aorta
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;
std::vector<size_t> windkessel_id;

double _TUBE_RADIUS_ = 0.015;
double _UMAX_ = 1.5;

void ExampleFile()
{
  Output::print("Example: Tube.h");
}


void InitialU1(double, double, double, double *values)
{
  
    values[0] = 0;
}


void InitialU2(double, double, double, double *values)
{
  
    values[0] = 0;
}


void InitialU3(double, double, double, double *values)
{

    values[0] = 0;
}


void InitialP(double, double, double, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP = unknown_solution_3d;

// ========================================================================
// boundary conditions
// ========================================================================


// kind of boundary condition (for FE space needed)
void BoundCondition(int reference, double , double , double , BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (reference == (int)neumann_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if (reference == (int)windkessel_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }
  

  // set Nitsche BC
  for (unsigned int j = 0; j < nitsche_id.size(); j++)
  {
    if (reference == (int)nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;
      return;
    }
  }
          
}

void U1BoundValue(int reference, double, double, double, double &value) // (int BdComp, double Param, double &value)
{
  switch (reference) {
  case 1:
    value = 0.;
    break;
  case 2:
  case 3:
    value = 0.;
    break;
  default:
     Output::print(" ** ERROR: no U1BoundValue for bd. reference = ",reference);
     break;
   }

  // reset value = 0 for Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (reference == (int)neumann_id[j])
    {
      value = 0.;
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if (reference == (int)windkessel_id[j])
    {
      value = 0.;
      return;
    }
  }
}

void U2BoundValue(int reference, double, double, double, double &value)
{
  switch (reference) {
  case 1:
    value = 0.;
    break;
  case 2:
  case 3:
    value = 0.;
    break;
  default:
     Output::print(" ** ERROR: no U2BoundValue for bd. reference = ",reference);
     break;
   }

  // reset value = 0 for Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (reference == (int)neumann_id[j])
    {
      value = 0.;
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if (reference == (int)windkessel_id[j])
    {
      value = 0.;
      return;
    }
  }
}

void U3BoundValue(int reference, double x, double y, double, double &value)
{
  double rMAX = _TUBE_RADIUS_;
  double r2 = x*x + y*y;
  double uMAX = _UMAX_;

  // time smoothing
  double t=TDatabase::TimeDB->CURRENTTIME;
  double T1 = 0.5;
  
  double sint=1.;
  if (t<T1) {
    sint = std::sin(2*M_PI/(4*T1)*t);
  }
  
  switch (reference) {
  case 1:
    // paraboloid
    value = sint*uMAX*(1 - r2/(rMAX*rMAX));
    
    break;
  case 2:
  case 3:
    value = sint*uMAX*(1 - r2/(rMAX*rMAX));
    break;
  default:
     Output::print(" ** ERROR: no U3BoundValue condition for bd. reference = ",reference);
     break;
   }

  // reset value = 0 for Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (reference == (int)neumann_id[j])
    {
      switch (reference) {
      case 1:
	value = 0.;
	break;
      case 3:
      case 2:
	value = 0.;
	break;
      }
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if (reference == (int)windkessel_id[j])
    {
      value = 0.;
      return;
    }
  }
}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][5] = 0.; // sigma (Brinkman)

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][4] = 0; 
  }
}


