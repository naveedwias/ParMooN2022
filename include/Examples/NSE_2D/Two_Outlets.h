// Navier-Stokes problem, Poiseuille-like with two outlets

double effective_viscosity = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;
std::vector<size_t> windkessel_id;

double UMAX = 1;

void ExampleFile()
{
  Output::print<1>("Example: Two_Outlets.h");
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_2d;
auto& ExactU2 = unknown_solution_2d;
auto& ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double , BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (i == (int)neumann_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if (i == (int)windkessel_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }
  

  // set Nitsche BC
  for (unsigned int j = 0; j < nitsche_id.size(); j++)
  {
    if (i == (int)nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;
      
      return;
    }
  }
}

void U1BoundValue(int BdComp, double, double &value)
{
  // loop to impose Neumann boundary conditions
  ///@attention we set =0 here, as Neumann BC are imposed using boundary assembling
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if ( BdComp == (int)neumann_id[j])
    {
      value = 0.;
      return;
    }
  }

  for (unsigned int j = 0; j < windkessel_id.size(); j++)
  {
    if ( BdComp == (int)windkessel_id[j])
    {
      value = 0.;
      return;
    }
  }

  // if we do not have Neumann or Windkessel, then we impose the Dirichlet value
  switch(BdComp)
  {
  case 0:
  case 2:
  case 3:
  case 5:
    value = 0;
    break;
  case 1:
    value = 0.333 * UMAX / 0.2;
    break; 
  case 4:
    value = 0.667 * UMAX / 0.2;
    break;
  case 6:
    value = UMAX; //4*UMAX*Param*(1-Param);
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void U2BoundValue(int, double, double &value)
{
  value = 0.;
  
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{

  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = effective_viscosity;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0.;//g
    // additional coefficient (used only in the Brinkman problem)
    coeff[4] = 0.;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void NonLinCoeffs(int n_points, double *, double *,
                  double **parameters, double **coeffs)
{
  
  double *coeff, *param;

  for( int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = effective_viscosity;
    coeff[1] = param[0];
    coeff[2] = param[1];

  }
}

