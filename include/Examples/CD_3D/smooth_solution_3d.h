//Example 5 from Allendes, Barrenechea and Rankin, SIAM-SCI17

double PECLET_NUMBER;

void ExampleFile()
{
  Output::root_info("Example", "smooth_solution_3d.h");
}

// exact solution is not known
void Exact(double x, double y, double z, double *values)
{
  values[0] = x*y*z*(1-x)*(1-y)*(1-z); //original solution
  values[1] = y*z*(1-x)*(1-y)*(1-z)-x*y*z*(1-y)*(1-z);  //x-derivative
  values[2] = x*z*(1-x)*(1-y)*(1-z)-x*y*z*(1-x)*(1-z);  //y-derivative
  values[3] = x*y*(1-x)*(1-y)*(1-z)-x*y*z*(1-x)*(1-y);  //z-derivative
  values[4] = -2*y*z*(1-y)*(1-z)-2*x*z*(1-x)*(1-z)-2*x*y*(1-x)*(1-y);    //laplacian
}

// kind of boundary condition (needed for FE space)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double, double, double, double &value)
{
  value=0;
}

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *z, const double *const*, double **coeffs)
{
  double exact[5];
  for(int i = 0; i < n_points; ++i)
  {
    coeffs[i][0] = PECLET_NUMBER; //diffusion coefficient
    coeffs[i][1] = 1;   //x coordinate of convection b
    coeffs[i][2] = 1;   //y coordinate of convection b
    coeffs[i][3] = 1;   //z coordinate of convection b
    coeffs[i][4] = 1;   //reaction coefficient c
    
    
    Exact(x[i],y[i],z[i],exact);
    coeffs[i][5] = -coeffs[i][0]*exact[4]; // diffusion
    coeffs[i][5] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2] + coeffs[i][3]*exact[3]; // convection
    coeffs[i][5] += coeffs[i][4]*exact[0]; // reaction
  }
}

