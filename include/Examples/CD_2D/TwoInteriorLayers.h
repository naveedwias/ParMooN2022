// ======================================================================
// sharp characteristic interior layer
// Knopp, Lube, Rapin, CMAME 2002
// ======================================================================
#define __TWO_INTERIOR_LAYERS__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: TwoInteriorLayers.h");
}
// exact solution
auto& Exact = unknown_solution_2d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double, BoundCond &cond)
{
    if (i==3)
	cond = NEUMANN;
    else
	cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   if (BdComp==0)
   {
      if ((Param>1.0/3.0)&& (Param<2.0/3.0))
         value = 1;
      else
         value = 0;
   }
   else
      value = 0;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = DIFFUSION;
  for(int i = 0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    coeffs[i][0] = eps;
    coeffs[i][1] = -y;
    coeffs[i][2] = x;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = std::sqrt(x*x + y*y);
  }
}


 
