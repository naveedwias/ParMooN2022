/**
 * Simple example for Navier-Stokes 3D, used for development, debugging and testing.
 * This example has two nice features: it calculates the boundary values
 * (Dirichlet only!) and the right hand side matching the intended exact
 * solution.
 *
 * Pressure solution is adapted to Default_UnitCube ([0,1]^3) geometry.
 *
 * @author Clemens Bartsch
 * @date 2016/04/04
 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","test_u_2_p_1.");
}

// exact solution
void ExactU1(double, double y, double, double *values)
{
  values[0] =  y*y;
  values[1] =  0;
  values[2] =  2*y;
  values[3] =  0;
  values[4] =  2; //Laplacien
}
void ExactU2(double x, double, double z, double *values)
{
  values[0] =  x*z;
  values[1] =  z;
  values[2] =  0;
  values[3] =  x;
  values[4] =  0; //Laplacien
}
void ExactU3(double, double y, double, double *values)
{
  values[0] =  -y*y;
  values[1] =  0;
  values[2] =  -2*y;
  values[3] =  0;
  values[4] =  -2; //Laplacien
}

void ExactP(double x, double y, double z, double *values)
{
  values[0] = x+y+z - 1.5; //adapted to be L^2_0 on unit cube
  values[1] = 1;
  values[2] = 1;
  values[3] = 1;
  values[4] = 0;
}

/* ****
 From here it's the same for all NSE3D test Examples.
 * **** */
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
  std::vector<double> u1(5,0.0);
  std::vector<double> u2(5,0.0);
  std::vector<double> u3(5,0.0);
  std::vector<double> p(5,0.0);

  for(int i=0;i<n_points;i++)
  {
    // get the space coordinates
    double x = X[i];
    double y = Y[i];
    double z = Z[i];

    // load the values vectors with exact values
    ExactU1(x,y,z,&u1.at(0));
    ExactU2(x,y,z,&u2.at(0));
    ExactU3(x,y,z,&u3.at(0));
    ExactP(x,y,z,&p.at(0));

    coeffs[i][0] = eps;
    coeffs[i][1] = -eps*u1[4] + ( u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3] ) + p.at(1) ; // f1
    coeffs[i][2] = -eps*u2[4] + ( u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3] ) + p.at(2) ; // f2
    coeffs[i][3] = -eps*u3[4] + ( u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3] ) + p.at(3) ; // f3
    coeffs[i][4] = 0; // g watch out that u is divergence free!
    coeffs[i][5] = 0; // sigma
  }
}
