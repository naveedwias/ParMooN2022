// ======================================================================
// Rotating Convection Field from ABR17
// ======================================================================
#define __rotating_convection_field__

double DIFFUSION;

void ExampleFile()
{
  Output::print("Example: rotating_convection_field.h");
}
// exact solution (this is the solution for eps = 0)
auto& Exact = unknown_solution_2d;

// kind of boundary condition
void BoundCondition(int, double, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    //diffusion coeffecient
    double eps;
    eps = DIFFUSION;
    switch(BdComp)
    {
	case 0:
        if(Param > 0 && Param < 0.0001)
            value = Param/eps;
        else if(Param >= 0.0001 && Param <= 0.4999)
            value = 1;
        else if(Param > 0.4999 && Param <0.5)
            value = (0.5-Param)/eps;
        else
            value = 0;
        break;
	case 1:
	    value = 0;
	    break;
	case 2:
		value = 0;
	    break;
	case 3:
		value = 0;
	    break;
    }
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=DIFFUSION;
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    double x = X[i];
    double y = Y[i];
    coeff = coeffs[i];
    //param = parameters[i];

    coeff[0] = eps;
    coeff[1] = -y; //x-coordiante of convection b
    coeff[2] = x; //y-coordiante of convection b
    coeff[3] = 0; //reaction coeffecient c
    coeff[4] = 0; //rhs
  }
}
