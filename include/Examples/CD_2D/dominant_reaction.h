// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================
#define __HMM_1986__

double DIFFUSION;

void ExampleFile()
{
  Output::print("Example: JK21_dominant_reaction.h");
}
// exact solution (this is the solution for eps = 0)
auto& Exact = unknown_solution_2d;

// kind of boundary condition
void BoundCondition(int, double, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
bool modified_HMM;

void BoundValue(int BdComp, double Param, double &value)
{
    value  = 0.0;
}

void BilinearCoeffs(int n_points, const double*, const double*,
                    const double*const*, double **coeffs)
{
  static double eps=DIFFUSION;
  int i;
  double *coeff;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0.004;
    coeff[2] = 0.012;
    coeff[3] = 1.0;
    coeff[4] = 1.0;
  }
}

