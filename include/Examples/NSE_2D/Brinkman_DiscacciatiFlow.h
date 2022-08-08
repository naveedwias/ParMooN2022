/******************************************************
 * Brinkman 2D problem:
 * Poiseuille-Problem from the thesis of Discacciati
 * Domain =  [0,12]x[0,8]
 * unknown solution
 * authors: Alfonso Caiazzo and Laura Blank
 ******************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1.;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Discacciati_Flow.h");
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

void BoundCondition(int i, double, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  if (i == 1)
    cond = NEUMANN;


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

void U1BoundValue(int BdComp, double Param, double &value)
{
  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be always zero!!!!
 /*
  for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if ( BdComp == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      switch(BdComp)
      {
      case 1:
        value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
        break;
      default:
        Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
        exit(1);
        break;
      }
      return;
    }
  }
  */

  double UMAX;
  // loop to impose (strong or weak) Dirichlet
  switch(BdComp)
  {
  case 0: value = 0.;
  break;
  case 1: value = 0.;
  break;
  case 2: value = 0.;
  break;
  case 3:
    UMAX = 0.1; // m/s
    value = 4. * UMAX * (1.-Param) * (Param);
    break;
  default:
    Output::print("No boundary component with this number.");
    exit(1);
    break;
  }
}

void U2BoundValue(int, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Brinkman problem:
// mueff, sigma := mu/permeability, f1, f2, g,
// with:
// -mueff Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// (lhs and rhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================
void LinCoeffs(int n_points, const double *x, const double *y,
               const double *const*, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][4] = sigma; // = mu/K

    // Discacciati region parameters
    if (   ( ( ( (x[i] >=7  ) && (x[i] <= 8  ) )  && ( (y[i] >= 4) && (y[i] <= 6) ) ) )
        || ( ( ( (x[i] >= 9) && (x[i] <= 10) )  && ( (y[i] >= 2) && (y[i] <= 4.2) ) ) )  )
    {
      coeffs[i][0]= 0.;
      coeffs[i][4]= 1000000.;
    }
    else if ( ( (x[i] >=5  ) && ( (y[i] >= 2) && (y[i] <= 6) ) ) )
    {
      coeffs[i][0] = 0.;
      coeffs[i][4]= 100.;
    }
    else
    {
      coeffs[i][0] = 0.01;
      coeffs[i][4]= 0.;
    }

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = 0.;
  }
}



