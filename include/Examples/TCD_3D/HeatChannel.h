// ======================================================================
// Sine problem 3D
// ======================================================================
#ifndef HEATCHANNEL_H
#define HEATCHANNEL_H



void ExampleFile()
{
  Output::print<1>("Example: HeatChannel.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  if (std::abs(z)<1e-8)
  { cond = DIRICHLET; }
  else
  { cond  = NEUMANN; }
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  value = 0;
  if (std::abs(z)<1e-8)
  {
    value =  1;
  }
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  if (std::abs(z)<1e-8)
  {
    values[0] = 1;
  }
  else
  {
    values[0] = 0;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double u = TDatabase::ParamDB->P0;
  static double eps = TDatabase::ParamDB->P1;//1.9e-5;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = u;
    coeff[4] = 0;
    coeff[5] = 0.;
  }
}

#endif
