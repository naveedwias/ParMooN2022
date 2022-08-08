// Time Convection Diffusion 2D problem for a hump changing its height
// circular inner layer
// John, Maubach, Tobiska, Num. Math. 1997

#ifndef __JOHNMAUBACHTOBISKA__
#define __JOHNMAUBACHTOBISKA__

namespace john_maubach_tobiska_inst
{
  double DiffCoeff = 1.;
  const double r0 = 0.25;
  const double x0 = 0.5;
  const double y0 = 0.5;

  // ===========================================================================
  // example file
  // ===========================================================================
  void ExampleFile()
  {
    Output::print("Example: JohnMaubachTobiska1997inst.h "
                  "(hump changing its height)") ;
    if(TDatabase::ParamDB->INTERNAL_QUAD_RULE != 99)
    {
      Output::print("INTERNAL_QUAD_RULE set to 99");
    }
  }

  // ===========================================================================
  // exact solution (this is the solution for eps = 0)
  // ===========================================================================
  void Exact(double x, double y, double* values)
  {
    double t = TDatabase::TimeDB->CURRENTTIME;
    double d = 2.0 * std::sqrt(1./DiffCoeff);

    double rxy = (x-x0)*(x-x0) + (y-y0)*(y-y0);
    double cxy = d*(r0*r0 - rxy);
    double hxy = 2.0*d / (1 + cxy*cxy);

    double g1xy    = std::atan(cxy) + M_PI/2.0;
    double g1xy_x  = -hxy*(x-x0);
    double g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
    double g1xy_y  = -hxy*(y-y0);
    double g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
    g1xy    = g1xy / M_PI; 
    g1xy_x  = g1xy_x / M_PI;
    g1xy_xx = g1xy_xx / M_PI;
    g1xy_y  = g1xy_y / M_PI;
    g1xy_yy = g1xy_yy / M_PI;

    double scale   = 16.0;
    double g2xy    = scale*x*(x - 1)*y*(y - 1);
    double g2xy_x  = scale*(2*x - 1)*y*(y - 1);
    double g2xy_xx = scale*2*y*(y - 1);
    double g2xy_y  = scale*(2*y - 1)*x*(x - 1);
    double g2xy_yy = scale*2*x*(x - 1);

    double gxy_xx = g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx;
    double gxy_yy = g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy;

    values[0] = std::sin(M_PI*t)*g1xy*g2xy;                     // u
    values[1] = std::sin(M_PI*t) * (g1xy_x*g2xy + g1xy*g2xy_x); // u_x
    values[2] = std::sin(M_PI*t) * (g1xy_y*g2xy + g1xy*g2xy_y); // u_y
    values[3] = std::sin(M_PI*t) * (gxy_xx + gxy_yy);           // u_xx + u_yy
  }

  // ===========================================================================
  // initial conditon
  // ===========================================================================
  void InitialCondition(double x, double y, double* values)
  {
    double val[4];
    Exact(x, y, val);
    values[0] = val[0];
  }

  // ===========================================================================
  // kind of boundary condition
  // ===========================================================================
  void BoundCondition(int, double, BoundCond& cond)
  {
    cond = DIRICHLET;
  }

  // ===========================================================================
  // value of boundary condition
  // ===========================================================================
  void BoundValue(int, double, double& value)
  {
    value = 0.;
  }

  // ===========================================================================
  // coefficients
  // ===========================================================================
  void BilinearCoeffs(int n_points, const double* X, const double* Y,
                      const double* const*, double** coeffs)
  {
    double* coeff;
    double t = TDatabase::TimeDB->CURRENTTIME;
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    // previous discrete time
    double t_m1 = t-tau;

    double b1 = 2.;
    double b2 = 3.;
    double c = 1.;

    double d  = 2.0 * std::sqrt(1.0/DiffCoeff);

    for(int i=0 ; i<n_points ; i++)
    {
      coeff = coeffs[i];
      double x = X[i];
      double y = Y[i];

      double rxy =  (x-x0)*(x-x0) + (y-y0)*(y-y0);
      double cxy =  d*(r0*r0 - rxy);
      double hxy =  2.0*d / (1 + cxy*cxy);

      double g1xy    = std::atan(cxy) + M_PI/2.0;
      double g1xy_x  = -hxy*(x-x0);
      double g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
      double g1xy_y  = -hxy*(y-y0);
      double g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
      g1xy    = g1xy / M_PI; 
      g1xy_x  = g1xy_x / M_PI;
      g1xy_xx = g1xy_xx / M_PI;
      g1xy_y  = g1xy_y / M_PI;
      g1xy_yy = g1xy_yy / M_PI;

      double scale   = 16.0;
      double g2xy    = scale*x*(x - 1)*y*(y - 1);
      double g2xy_x  = scale*(2*x - 1)*y*(y - 1);
      double g2xy_xx = scale*2*y*(y - 1);
      double g2xy_y  = scale*(2*y - 1)*x*(x - 1);
      double g2xy_yy = scale*2*x*(x - 1);

      double gxy    = std::sin(M_PI*t)*g1xy*g2xy;
      double gxy_x  = std::sin(M_PI*t)*(g1xy_x*g2xy + g1xy*g2xy_x);
      double gxy_xx = std::sin(M_PI*t)
                      * (g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx);
      double gxy_y  = std::sin(M_PI*t)*(g1xy_y*g2xy + g1xy*g2xy_y);
      double gxy_yy = std::sin(M_PI*t)
                      * (g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy);
      double gxy_l  = gxy_xx + gxy_yy;

      coeff[0] = DiffCoeff;
      coeff[1] = b1;
      coeff[2] = b2;
      coeff[3] = c;
      coeff[4] = M_PI*std::cos(M_PI*t)*g1xy*g2xy - coeff[0]*gxy_l
                 + coeff[1]*gxy_x + coeff[2]*gxy_y + coeff[3]*gxy;

      // for schemes needing the residual, residual of previous time,
      // without time derivative
      gxy    = std::sin(M_PI*t_m1)*g1xy*g2xy;
      gxy_x  = std::sin(M_PI*t_m1)*(g1xy_x*g2xy + g1xy*g2xy_x);
      gxy_xx = std::sin(M_PI*t_m1)
                      * (g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx);
      gxy_y  = std::sin(M_PI*t_m1)*(g1xy_y*g2xy + g1xy*g2xy_y);
      gxy_yy = std::sin(M_PI*t_m1)
                      * (g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy);
      gxy_l  = gxy_xx + gxy_yy;
      coeff[5] = M_PI*std::cos(M_PI*t_m1)*g1xy*g2xy - coeff[0]*gxy_l
                 + coeff[1]*gxy_x + coeff[2]*gxy_y + coeff[3]*gxy; 

      // time derivative
      coeff[6] = (std::sin(M_PI*t)*g1xy*g2xy - std::sin(M_PI*t_m1)*g1xy*g2xy)
                 / t_m1;
    }
  }

} // namespace john_maubach_tobiska_inst

#endif
