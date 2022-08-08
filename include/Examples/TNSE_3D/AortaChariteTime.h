/*******************************************************************************
 * Setup for time-dependent simulation of flow through ascending aorta
 ******************************************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;
std::vector<size_t> windkessel_id;


void ExampleFile()
{
  Output::print<1>("Example: AortaChariteTime.h");
}

// =============================================================================
// exact solution
// =============================================================================
auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP  = unknown_solution_3d;

// =============================================================================
// initial condition
// =============================================================================
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

// =============================================================================
// boundary conditions
// =============================================================================
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

void U1BoundValue(int reference, double, double, double, double &value)
{
  switch (reference)
  {
    case 7: // inlet
      value = 0.;
      break;
    case 2: // outlets
    case 3:
    case 4:
    case 5:
      value = 0.;
      break;
    case 1: // wall
      value = 0.;
      break;
#ifdef _MPI
    case 0:
      break;
#endif
    default:
      Output::print(" ** ERROR: no condition for bd. reference = ",reference);
      break;
  }
}

void U2BoundValue(int reference, double, double, double, double &value)
{
  switch (reference)
  {
    case 7: // inlet
      value = 0.;
      break;
    case 2: // outlets
    case 3:
    case 4:
    case 5:
      value = 0.;
      break;
    case 1: // wall
      value = 0.;
      break;
#ifdef _MPI
    case 0:
      break;
#endif
    default:
      Output::print(" ** ERROR: no condition for bd. reference = ",reference);
      break;
  }
}

void U3BoundValue(int reference, double , double , double, double &value)
{
  switch (reference)
  {
    case 7: // inlet
      value = 0.;
      break;
    case 2: // outlets
    case 3:
    case 4:
    case 5:
      value = 0.;
      break;
    case 1: // wall
      value = 0.;
      break;
#ifdef _MPI
    case 0:
       break;
#endif
    default:
      Output::print(" ** ERROR: no condition for bd. reference = ",reference);
      break;
  }
}


// =============================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// =============================================================================
void LinCoeffs(int n_points, const double*, const double*, const double*,
               const double*const*, double** coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameter
    coeffs[i][0] = effective_viscosity;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][4] = 0; 
  }
}

