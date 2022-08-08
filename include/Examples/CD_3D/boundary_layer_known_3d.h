//Example 6 from Allendes, Barrenechea and Rankin, SIAM-SCI17

double PECLET_NUMBER;

void ExampleFile()
{
  Output::root_info("Example", "boundary_layer_known_3d.h");
}

// exact solution is not known
void Exact(double x, double y, double z, double *values)
{
  double p = PECLET_NUMBER; 
  double fx = std::exp( (x - 1.0) / p );
  double a = std::exp( - 1.0 / p );
  double fxTerm = x - (fx - a) / (1.0 - a); 
  double denominator = (1.0 - a)*p;
  
  values[0] = y*z*(1.0 - y)*(1.0 - z)*fxTerm;
  values[1] = y*z*(1.0 - y)*(1.0 - z)*( 1.0 - fx / denominator );
  values[2] = (1.0 - 2.0*y)*z*(1.0 - z)*fxTerm;
  values[3] = (1.0 - 2.0*z)*y*(1.0 - y)*fxTerm;
  values[4] = - y*z*(1.0 - y)*(1.0 - z)*( fx / (denominator*p) )
              - 2.0*( y*(1.0 - y) + z*(1.0 - z) )*fxTerm; 
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
    coeffs[i][2] = 0;   //y coordinate of convection b
    coeffs[i][3] = 0;   //z coordinate of convection b
    coeffs[i][4] = 0;   //reaction coefficient c
    
    
    Exact(x[i],y[i],z[i],exact);
    coeffs[i][5] = -coeffs[i][0]*exact[4]; // diffusion
    coeffs[i][5] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2] + coeffs[i][3]*exact[3]; // convection
    coeffs[i][5] += coeffs[i][4]*exact[0]; // reaction
  }
}

