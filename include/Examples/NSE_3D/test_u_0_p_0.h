/**
 * Simple example for Navier-Stokes 3D, used for development, debugging and testing.
 * Should work for each geometry (when Dirichlet conditions only!).
 * Note: Got equipped with Neumann boundary conditions (assuming the unit cube
 * [0,1]^3) on the front face (x==1).
 *
 * @note: There is a note in the ParMooN to do directory on what to regard when
 * coupling Neumann and Dirichlet boundary conditions in a ParMooN 3D examples.
 *
 * Constant velocity solution u = (1,0,0).
 * Constant pressure solution p = 0.
 * Dirichlet boundary conditions, and one cube facet with Neumann boundaries .
 * 
 * @author ???, Clemens Bartsch imported this from MooNMD
 * @date 2016/03/11 Import to ParMooN.
 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","test_u_0_p_0. Constant velocity solution u=(1,0,0), "
	 "constant pressure solution p=0. \n");
}

// exact solution
void ExactU1(double, double, double, double *values)
{
  values[0] =  1;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacian
}
void ExactU2(double, double, double, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacian
}
void ExactU3(double, double, double, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacian
}

void ExactP(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// Helper function - finds out if (x,y,z) is on the interior of the face
// with x==1 (the "front face", if one imagines it in Cartesian system).
bool on_front_face_interior(double x, double y, double z)
{
  double tol = 1e-10;
  return std::abs(1-x) < tol
      && !(std::abs(y -1) < tol)  && !(std::abs(y) < tol)
      && !(std::abs(z -1) < tol)  && !(std::abs(z) < tol);
}

/* ****
 From here it's the same for all NSE3D test Examples.
 * **** */
/**
 * Helpful for adjusting NSE Neumann bdry conditions to a given analytical
 * solution - ONLY LAPLACETYPE 0 implemented!
 * Calculates S*n, where S is the Cauchy stress tensor. Neumann bdry conditions
 * in NSE mean to prescribe the value of S*n.
 * For Laplacetype 0 we have S = nu * grad(u) - p*I, and n is the outer normal
 * of the boundary part.
 * @param[in] x, y, z The boundary point.
 * @param[in] n0, n1, n2 The outer normal on the bdry part in (x,y,z).
 * @param[in] i The column of the result we are interested in.
 */

double S_times_n(
    double x, double y, double z,
    double n0, double n1, double n2, size_t i)
{
  //load the desired exact values of u and p
  double u1[5] = {0,0,0,0,0};
  double u2[5] = {0,0,0,0,0};
  double u3[5] = {0,0,0,0,0};
  double p[5] = {0,0,0,0,0};
  ExactU1(x,y,z, u1);
  ExactU2(x,y,z, u2);
  ExactU3(x,y,z, u3);
  ExactP(x,y,z,p);

  //set up the gradient matrix
  double du[3][3];
  du[0][0] = u1[1];
  du[0][1] = u1[2];
  du[0][2] = u1[3];
  du[1][0] = u2[1];
  du[1][1] = u2[2];
  du[1][2] = u2[3];
  du[2][0] = u3[1];
  du[2][1] = u3[2];
  du[2][2] = u3[3];

  // the outer normal entries as a vector
  double n[3] = {n0, n1, n2};

  double nu = DIMENSIONLESS_VISCOSITY;

  double sum = 0;
  for (int k=0; k<3; ++k)
  {
    //TODO Remove this global Database dependency!
    if(TDatabase::ParamDB->LAPLACETYPE == 0)
      sum += nu * du[i][k] * n[k];
    else if(TDatabase::ParamDB->LAPLACETYPE == 1)
      sum += nu * (du[i][k] + du[k][i])* n[k];
  }
  sum -= p[0]*n[i]; //subtract pressure part
  return sum;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{

  if(on_front_face_interior(x,y,z)) //Neumann only on interior of the plane with x == 1
    cond = NEUMANN;
  else
  {
    cond = DIRICHLET;
  }
}
// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  if(on_front_face_interior(x,y,z))
    value = S_times_n(x,y,z, 1,0,0, 0);
  else
  {
  double diri[5];
  ExactU1(x,y,z,diri);
  value = diri[0]; //Dirichlet value
  }
}
void U2BoundValue(int, double x, double y, double z, double &value)
{
  if(on_front_face_interior(x,y,z))
    value = S_times_n(x,y,z, 1,0,0, 1);
  else
  {
  double diri[5];
  ExactU2(x,y,z,diri);
  value = diri[0]; //Dirichlet value
  }
}
void U3BoundValue(int, double x, double y, double z, double &value)
{
  if(on_front_face_interior(x,y,z))
    value = S_times_n(x,y,z, 1,0,0, 2);
  else
  {
    double diri[5];
    ExactU3(x,y,z,diri);
    value = diri[0]; //Dirichlet value
  }
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
