// ======================================================================
// Sine problem
// ======================================================================

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: SineLaplace.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
  const double p = M_PI; // 2*M_PI;
  values[0] = std::sin(p*x)*std::sin(p*y);
  values[1] = p*std::cos(p*x)*std::sin(p*y);
  values[2] = p*std::sin(p*x)*std::cos(p*y);
  values[3] = -2*p*p*std::sin(p*x)*std::sin(p*y);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double, BoundCond &cond)
{
  if(BdComp==1)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

void transform(const int BdComp, const double Param, double& x, double& y, 
               double& nx, double& ny)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      nx = 0.;
      ny = -1.;
      break;
    case 1:
      x = 1.;
      y = Param;
      nx = 1.;
      ny = 0.;
      break;
    case 2:
      x = 1. - Param;
      y = 1.;
      nx = 0.;
      ny = 1.;
      break;
    case 3:
      x = 0.;
      y = 1. - Param;
      nx = -1.;
      ny = 0.;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  double exact[4];
  Exact(x, y, exact);
  if(cond == NEUMANN)
  {
    const double eps = DIFFUSION;
    value = eps * (nx * exact[1] + ny * exact[2]);
  }
  else // Dirichlet
    value = exact[0];
}

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *const*, double **coeffs)
{
  const double eps = DIFFUSION;
  double exact[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    coeffs[i][1] = 0;// parameters[i][0]; 
    coeffs[i][2] = 0;//parameters[i][1];
    coeffs[i][3] = 0;
    
    Exact(x[i], y[i], exact);

    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction

  // coeffs[i][4] = 0;
  }
}

