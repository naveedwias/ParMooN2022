//Simple Example on the cube [0,1]^3
//with Dirichlet boundary condition

void ExampleFile()
{
  Output::print("A simple SinCos Darcy program on [0,1]^3 with Dirichlet ",
                "boundary condition using mixed finite elements") ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[1] = M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[2] = M_PI*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[3] = M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[4] = -3*M_PI*M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
}

void ExactU2(double x, double y,double z, double *values)
{
  values[0] = -M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[1] = M_PI*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[2] = M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[3] = -M_PI*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z);
  values[4] = 3*M_PI*M_PI*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[1] = M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[2] = -M_PI*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z);
  values[3] = M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[4] = 3*M_PI*M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
}


void ExactP(double x, double y, double z, double *values)
{
  values[0] = std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[1] = -M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[2] = M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[3] = M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[4] = -3*M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void FluxBoundValue(int, double x, double y, double z, double &value)
{
  double tol = 1e-10;
   bool pointOnBoundary = true;
   if(std::abs(x) <= tol) // (x == 0.0)
   {
     if( (y>=0 && y<=1) && ( z>=0 && z<=1 ))
       value = 0; // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(x-1.0) <= tol) //(x == 1.0)
   {
     if( (y>=0 && y<=1) && ( z>=0 && z<=1 ))
       value = 0; // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(y) <= tol) // (y == 0.0)
   {
     if( (x>=0 && x<=1) && ( z>=0 && z<=1 ))
       value = M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(y-1.0) <= tol) // (y == 1.0)
   {
     if( (x>=0 && x<=1) && ( z>=0 && z<=1 ))
       value = -M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(z) <= tol) // (z == 0.0)
   {
     if( (x>=0 && x<=1) && ( y>=0 && y<=1 ))
       value = M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(z-1.0) <= tol) // (z == 1.0)
   {
     if( (x>=0.0 && x<=1.0) && ( y>=0.0 && y<=1.0 ))
       value = -M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else
     pointOnBoundary = false;

   if(!pointOnBoundary)
   {
     ErrThrow("trying to evaluate boundary data at a point not belonging to "
              "the boundary, (", x, ",", y, ",", z, ")");
   }
}

void BoundConditionPressure(int, double, double, double, BoundCond &cond)
{
  cond = NEUMANN; 
}

void PressureBoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// coefficients in the pde
void LinCoeffs(int n_points, const double *X, const double *Y, const double *Z,
               const double *const*, double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = eps;      // permeability
    coeffs[i][1] = 0;        // f1, forces
    coeffs[i][2] = 0;        // f2, forces
    coeffs[i][3] = 0;        // f3, forces
    //source terms:
    coeffs[i][4] = 3 * M_PI * M_PI * std::cos(M_PI*X[i]) * std::sin(M_PI*Y[i]) * std::sin(M_PI*Z[i]);
  }
}


