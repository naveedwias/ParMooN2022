/*********************************************************
 * Setup for simulation of flow through ascending aorta
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;
std::vector<size_t> windkessel_id;

// inlet normal for the mesh we are using
double nx = -0.105282895457120;
double ny = 0.232829669272479;
double nz = 0.966801870618088;

void ExampleFile()
{
  Output::print<1>("Example: Aorta.h");
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
  case 4:
  case 5:
  case 6:
  case 1:
     value = 0.;
     break;
     
     
   case 2:
     value = 0.;
     break;
     
   default:
     Output::print(" ** ERROR: no condition for bd. reference = ",reference);
     break;
   }
}

void U2BoundValue(int reference, double, double, double, double &value)
{
  switch (reference) {
   case 4:
   case 5:
   case 6:
     value = 0.;
     break;
     
   case 1:
     value = 1.*ny;
     break;
     
   case 2:
     value = 0.;
     break;
     
   default:
     Output::print(" ** ERROR: no condition for bd. reference = ",reference);
     break;
   }
}

void U3BoundValue(int reference, double , double , double, double &value)
{
  switch (reference) {
   case 4:
   case 5:
   case 6:
     value = 0.;
     break;
     
   case 1:
     value = 1.*nz;
     break;
     
   case 2:
     value = 0.;
     break;
     
   default:
     Output::print(" ** ERROR: no condition for bd. reference = ",reference);
     break;
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


