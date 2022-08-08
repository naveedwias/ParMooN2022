// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  Output::print<1>("Example: HeatChanel.h");
}

// exact solution
auto& Exact = unknown_solution_3d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double z, BoundCond &cond)
{
  if (std::abs(z)<1e-8)
  { cond = DIRICHLET; }
  else
  { cond  = NEUMANN; }
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
   double r2;
   
   value = 0;
       
   r2 = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
   if (std::abs(z)<1e-8 && r2<0.25)
   {
      value =  -4.*(0.25-r2);
   }
   
    r2 = (x+0.5)*(x+0.5)+(y+0.5)*(y+0.5);
   if (std::abs(z)<1e-8 && r2<0.25)
   {
      value =  -4.*(0.25-r2);
   }
   
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *, const double *const*, double **coeffs)
{
  static double eps=1.9e-5, u=TDatabase::ParamDB->P0;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = u;
    coeff[4] = 0;
    coeff[5] = 0.;
  }
}

// initial conditon
void InitialCondition(double, double, double, double *values)
{ 
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;  
}
