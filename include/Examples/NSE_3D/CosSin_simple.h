/**
 * @file A Navier--Stokes test problem with sine-cosine solution.
 *
 * The boundary data is adapted to the [0.1]^3 unit cube example, which is
 * availabe as default geometry in ParMooN. It will throw an error if you
 * try running it on any other domain - just to make you aware of that fact.

 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","Simple example with sin and cos solution and p=0");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =          std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[1] =      -M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[2] =       M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[3] =       M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[4] = -3*M_PI*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z); //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =          std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[1] =       M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[2] =      -M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[3] =       M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z);
  values[4] = -3*M_PI*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); //Laplacien
}
void ExactU3(double x, double y,  double z, double *values)
{
  values[0] =      -2*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[1] =   -2*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[2] =   -2*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z);
  values[3] =    2*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[4] = 6*M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z); //Laplacien
}

void ExactP(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}
// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU1(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U2BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU2(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U3BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU3(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}

void LinCoeffs(int n_points, const double * X, const double * Y,
               const double * Z, const double *const*, double **coeffs)
{
  const double eps = DIMENSIONLESS_VISCOSITY;
  double u1[5], u2[5], u3[5], p[5];
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] =  eps;

    ExactU1(X[i], Y[i], Z[i], u1);
    ExactU2(X[i], Y[i], Z[i], u2);
    ExactU3(X[i], Y[i], Z[i], u3);
    ExactP( X[i], Y[i], Z[i], p);

    coeffs[i][1] = -eps * u1[4] + p[1]; //Stokes: diffusion and pressure gradient
    coeffs[i][2] = -eps * u2[4] + p[2];
    coeffs[i][3] = -eps * u3[4] + p[3];
    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) //Navier--Stokes:
    {                                              // add convective terms
      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3];
      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3];
      coeffs[i][3] += u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3];
    }
    coeffs[i][4] = 0; // g
    coeffs[i][5] = 0; // sigma
  }
}

// From here it's old stuff: explicitely formulated dirichlet and neumann
// boundary conditions - if you want to reuse them, take care of the code!
// // kind of boundary condition (for FE space needed)
//void BoundCondition(double x, double y, double z, BoundCond &cond)
//{
//  double tol = 1e-10;
//  if((std::abs(1+x) < tol) || (std::abs(1+y) < tol) || (std::abs(1+z) < tol)
//       || (std::abs(1-z) < tol))
//    cond = DIRICHLET;
//  else
//    cond = NEUMANN;
//
//}
//
// // former values of boundary condition
//void U1BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((std::abs(1+x) < tol) || (std::abs(1+y) < tol) || (std::abs(1+z) < tol)
//       || (std::abs(1-z) < tol))
//    value = std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z); //Dirichlet
//  else{
//    if(std::abs(1-x) < tol)
//    {
//      value = -eps*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += -eps*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
//    }
//    else
//    {
//      value = eps*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
//    }
//  }
//}
//void U2BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((std::abs(1+x) < tol) || (std::abs(1+y) < tol) || (std::abs(1+z) < tol)
//       || (std::abs(1-z) < tol))
//    value = std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); //Dirichlet
//  else{
//    if(std::abs(1-x) < tol)
//    {
//      value = eps*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
//    }
//    else
//    {
//      value = -eps*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += -eps*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
//    }
//  }
//}
//void U3BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((std::abs(1+x) < tol) || (std::abs(1+y) < tol) || (std::abs(1+z) < tol)
//       || (std::abs(1-z) < tol))
//    value = -2*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z); //Dirichlet
//  else{
//    if(std::abs(1-x) < tol)
//    {
//      value = -eps*2*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
//    }
//    else
//      {
//      value = -eps*2*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::cos(M_PI*z);
//      }
//  }
//}
