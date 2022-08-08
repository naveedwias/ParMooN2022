// Bechmark for Darcy problem, exact solution is
// 

void ExampleFile()
{
  Output::print<1>("Example: SinSinSolution.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  values[1] = 4*M_PI*M_PI*std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
  values[2] = -4*M_PI*M_PI*std::cos(2*M_PI*x)*std::cos(2*M_PI*y);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
  values[1] = -4*M_PI*M_PI*std::cos(2*M_PI*x)*std::cos(2*M_PI*y);
  values[2] = 4*M_PI*M_PI*std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
  values[1] = 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  values[2] = 2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
  values[3] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double, BoundCond &cond)  
{
  cond = (bdComp == 0 || bdComp == 3) ? NEUMANN : DIRICHLET;
}

// u \cdot n
void FluxBoundValue(int bdComp, double t, double &value)
{
  switch(bdComp)
  {
    case 0: 
    {
      value = 0; // Neumann
      //value = 2*M_PI*std::sin(2*M_PI*t); // Dirichlet
      break;
    }
    case 1: 
    {
      //value = 0; // Neumann
      value = -2*M_PI*std::sin(2*M_PI*t); // Drichlet
      break;
    }
    case 2: 
    {
      //value = 0; // Neumann
      value = -2*M_PI*std::sin(2*M_PI*(1-t)); // Dirichlet
      break;
    }
    case 3:
    {
      value = 0; // Neumann
      //value = 2*M_PI*std::sin(2*M_PI*(1-t)); // Dirichlet
      break;
    }
    default: cout << "wrong boundary part number" << endl;
    
    break;
  }
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  const double eps = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 8*M_PI*M_PI*std::sin(2*M_PI*X[i])*std::sin(2*M_PI*Y[i]); // g
  }
}


