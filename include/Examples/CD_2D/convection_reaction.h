// ======================================================================
// Pure convection-reaction, no diffusion
// Houston, Schwab, Sueli, SIAM 39(6), 2002, doi: 10.1137/S0036142900374111
// ======================================================================

void ExampleFile()
{
  Output::print<1>("Example: convection_reaction.h") ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;
}

// exact solution
void Exact(double x, double y, double *values)
{
  const double p = M_PI; // Pi;

  double x_hat = 2*x-1; // transform from (0,1) to (-1,1)
  double y_hat = 2*y-1; // transform from (0,1) to (-1,1)

  double arg = p*(1+x_hat)*(1+y_hat)*(1+y_hat)/8; // argument of sin
  values[0] = 1+std::sin(arg);  //original solution

  double arg_x = p*(1+y_hat)*(1+y_hat)/4; // inner derivative wrt x
  double arg_y = p*(1+x_hat)*(1+y_hat)/2; // inner derivative wrt y

  values[1] = std::cos(arg)*arg_x; //x-derivative
  values[2] = std::cos(arg)*arg_y; //y-derivative
  double dxx = -std::sin(arg)*arg_x*arg_x;
  double dyy_1 = -std::sin(arg)*arg_y*arg_y;
  double dyy_2 = std::cos(arg)*p*(1+x_hat)/2;
  values[3] = dxx + dyy_1 + dyy_2 ; //laplacian
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double, BoundCond &cond)
{
  if(BdComp%3 == 0)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

//to get the boundary lengths of the example
void transform(const int BdComp, const double Param, double& x, double& y)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      break;
    case 1:
      x = 1.;
      y = Param;
      break;
    case 2:
      x = 1.-Param;
      y = 1.;
      break;
    case 3:
      x = 0.;
      y = 1.-Param;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  //coordinates for the part of boundary
  double x,y;
  transform(BdComp, Param, x, y);

  //type of boundary condition (DIRICHLET or NEUMANN)
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);

  //exact solution at the point (x,y), required for BC
  double u[4];
  Exact(x, y, u);

  if(cond == DIRICHLET)
  {
    value = u[0];
  }
  if(cond == NEUMANN)
  {
    value = 0;
  }
}

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *const*, double **coeffs)
{
  double *coeff;
  double exact[4];

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //param = parameters[i];

    double x_hat = 2*x[i]-1; // transform from (0,1) to (-1,1)
    double y_hat = 2*y[i]-1; // transform from (0,1) to (-1,1)

    coeff[0] = 0;
    coeff[1] = 2-y_hat*y_hat; //x-coordiante of convection b
    coeff[2] = 2-x_hat; //y-coordiante of convection b
    coeff[3] = 1+(1+x_hat)*(1+y_hat)*(1+y_hat); //reaction coeffecient c

    Exact(x[i], y[i], exact);

    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
  }
}
