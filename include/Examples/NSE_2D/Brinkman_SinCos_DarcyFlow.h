/************************************************
 *  Test for the Darcy limit of the Brinkman problem (mueff = 0)
 *  with sine-cosine solution taken from: DOI: https://doi.org/10.1137/08072632X.
 *
 *  exact solution:
 *  u(x,y) = [u1,u2]^T = [-2 pi std::cos(2 pi x) * std::sin(2 pi y),  -2 pi std::sin(2 pi x) * std::cos(2 pi y) ]^T
 *  p(x,y) = sigma * std::sin(2 Pi*x) * std::sin(2 Pi*y)
 *  u \cdot n = 0 on \partial \Omega
 *
 *  solves the Darcy Problem:
 *  sigma u + p = 0
 *  div u = 0
 *  u \cdot n = ... \partial Omega
 *
 *  The boundary condition is treated as an essential boundary condition.
 **************************************************/

// initialize physical parameters
double effective_viscosity = -1;
double sigma = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_SinCos_DarcyFlow.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -2 * M_PI * std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y);              //u1
  values[1] =  4 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);         //u1_x
  values[2] = -4 * M_PI * M_PI * std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y);         //u1_y
  values[3] = 16 * M_PI * M_PI * M_PI * std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y);    //Delta u1=u1_xx+u1_yy
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -2 * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);            //u2
  values[1] = -4 * M_PI * M_PI * std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y);       //u2_x
  values[2] =  4 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);       //u2_y
  values[3] = 16 * M_PI * M_PI * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);  //Delta u2=u2_xx + u2_yy
}

void ExactP(double x, double y, double *values)
{
  values[0] = sigma * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);                    //p
  values[1] = sigma * 2 * M_PI * std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y);           //p_x
  values[2] = sigma * 2 * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);           //p_y
  values[3] = sigma * (-8) * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);   //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions (Parametrisation of the boundary); Param \in [0,1]
// ========================================================================
void BoundCondition(int i, double, BoundCond &cond)
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
  switch(BdComp)
  {
  case 0: value = 0.;
  break;
  case 1: value = -2 * M_PI * std::sin( 2 * M_PI * Param );  // u \cdot n = u_1; x = 1 --> std::cos(2 pi x) = 1
  break;
  case 2: value = 0.;
  break;
  case 3: value = -2 * M_PI * std::sin( 2 * M_PI * (1-Param) );  // u \cdot n = - u_1; x = 0 --> std::cos(2 pi x) = 1
  break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: value = -2 * M_PI * std::sin( 2 * M_PI * Param);  // u \cdot n = - u_2; y = 0 --> std::cos(2 pi y) = 1
  break;
  case 1: value = 0.;
  break;
  case 2: value = -2 * M_PI * std::sin( 2 * M_PI * (1-Param));  // u \cdot n = u_2; y = 1 --> std::cos(2 pi y) = 1
  break;
  case 3: value = 0.;
  break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
}

// ========================================================================
// coefficients for Brinkman problem:
// mu, f1, f2, g, sigma = mu/permeability
// with:
// -mueff Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double *x, const double *y,
               const double *const*, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  //double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    //ExactP(x[i], y[i], val_p);

    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][4] = sigma;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0.; //-coeffs[i][0] * val_u1[3] + val_p[1] + coeffs[i][4] * val_u1[0];   // f1
    coeffs[i][2] = 0.; //-coeffs[i][0] * val_u2[3] + val_p[2] + coeffs[i][4] * val_u2[0];   // f2

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2]; // g
  }
}

