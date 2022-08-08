// ======================================================================
// Polynomial problem without diffusion for the UnitSquare with a jump
// discontinuity.
// This MAYBE works on other rectangles as well as long as the sides are ordered
// counterclockwise, the first side is parallel to the x-axis and below the
// third side.
// ======================================================================


unsigned int deg;
// Compute for a given degree the number of coefficients of polynom and set
// the coefficients as 1, -2, 3, -4, ..., -8, 9, -1, 2, -3, ...
std::vector<int> coef;
unsigned int n_coef ;

void ExampleFile()
{
  Output::print<1>("Example: jump-polynomial.h with polynomial degree ", deg);
  auto n_coef = (deg*deg + 3*deg)/2 + 1;
  coef.resize(n_coef, 0);
  for (unsigned int i = 0; i < n_coef; ++i) {
    coef[i]= 1;
    if (i%2 == 1)
    {
      coef[i] *= -1;  // this may lead to a slower blow up of the solution
    } // end if
  } // end for
}


// exact solution
void Exact(double x, double y, double *values)
{
  // Compute for a given degree the number of coefficients of polynom and set
  // the coefficients as 1, -2, 3, -4, ..., -8, 9, -1, 2, -3, ...

  // values[0] = function at (x,y)
  unsigned int index = 0;
  values[0] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      auto val = coef[index] * std::pow(x,i) * std::pow(y,j);
      if (y < 0.5)
      {
        values[0] -= val;
      }
      else
      {
        values[0] += val;
      }

      index++;
    }
  }

  // values[1] = derivative wrt. x of function at (x,y)
  index = 0;
  values[1] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (i == 0)
      {
        index++;
      }
      else
      {
        auto val = coef[index] * (int) i * std::pow(x,i-1) * std::pow(y,j);
        if (y < 0.5)
        {
          values[1] -= val;
        }
        else
        {
          values[1] += val;
        }
        index++;
      }
    }
  }

  // values[2] = derivative wrt. y of function at (x,y)
  index = 0;
  values[2] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (j == 0)
      {
        index++;
      }
      else
      {
        auto val = coef[index] * (int) j * std::pow(x,i) * std::pow(y,j-1);
        if (y < 0.5)
        {
          values[2] -= val;
        }
        else
        {
          values[2] += val;
        }
        index++;
      }
    }
  }

  // values[3] = Laplacian of function at (x,y)
  index = 0;
  values[3] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. x
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (i <= 1)
      {
        index++;
      }
      else
      {
        auto val = coef[index] * (int) i * (int) (i-1) * std::pow(x,i-2) * std::pow(y,j);
        if (y < 0.5)
        {
          values[3] -= val;
        }
        else
        {
          values[3] += val;
        }
        index++;
      }
    }
  }
  index = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. y
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (j <= 1)
      {
        index++;
      }
      else
      {
        auto val = coef[index] * (int) j * (int) (j-1) * std::pow(x,i) * std::pow(y,j-2);
        if (y < 0.5)
        {
          values[3] -= val;
        }
        else
        {
          values[3] += val;
        }
        index++;
      }
    }
  }
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double, BoundCond &cond)
{
  if(BdComp == 3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

void transform(const int BdComp, const double Param, double& x, double& y, 
               double& nx, double& ny)
{ // Here the UnitSquare is implicitly assumed
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
    value = 0;
  }
  else // Dirichlet
    value = exact[0];
}

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *const*, double **coeffs)
{
  double exact[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = 0; // diffusion parameter
    coeffs[i][1] = 1; // parameters[i][0], 1st component of advection function b
    coeffs[i][2] = 0; // parameters[i][1], 2nd component of advection function b
    coeffs[i][3] = 2; // reaction function c

    Exact(x[i], y[i], exact);

    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction
  }
}


