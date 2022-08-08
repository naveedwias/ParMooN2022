// Time Convection Diffusion 2D problem in a unit square domain [0,1]x[0,1],
// prescribed solution: u(x,y) = sin(PI*t) * sin(2*PI*x) * sin(2*PI*y)
// with convection field b = (1-y , -2/3+x) and reaction term c = 1

#ifndef __SIN3_0_PI__
#define __SIN3_0_PI__


namespace sin3_0_pi
{
  double DiffCoeff = 1.;

  // ===========================================================================
  // example file
  // ===========================================================================
  void ExampleFile()
  {
    int rank = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if(rank==0)
    {
      Output::print("\nExample: Sin3_0_pi.h");
      if(TDatabase::ParamDB->INTERNAL_QUAD_RULE != 99)
      {
        Output::print("INTERNAL_QUAD_RULE set to 99");
      }
    }
    TDatabase::ParamDB->INTERNAL_QUAD_RULE = 99;
  }

  // ===========================================================================
  // exact solution
  // ===========================================================================
  void Exact(double x, double y, double* values)
  {
    double t = TDatabase::TimeDB->CURRENTTIME;

    // u
    values[0] = std::sin(M_PI*t) *std::sin(2*M_PI*x) *std::sin(2*M_PI*y);
    // u_x
    values[1] = std::sin(M_PI*t) *2*M_PI*std::cos(2*M_PI*x) *std::sin(2*M_PI*y);
    // u_y
    values[2] = std::sin(M_PI*t) *std::sin(2*M_PI*x) *2*M_PI*std::cos(2*M_PI*y);
    // u_xx + u_yy
    values[3] = -8.*M_PI*M_PI * std::sin(M_PI*t)
                              * std::sin(2*M_PI*x) * std::sin(2*M_PI*y);
  }

  // ===========================================================================
  // initial conditon
  // ===========================================================================
  void InitialCondition(double x,  double y, double* values)
  {
    double exact[4];
    Exact(x, y, exact);
    values[0] = exact[0];
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
    double exact[4];
    double rhs;
    double t = TDatabase::TimeDB->CURRENTTIME;

    for(int i=0 ; i<n_points ; i++)
    {
      coeff = coeffs[i];
      double x = X[i];
      double y = Y[i];
      Exact(x, y, exact);

      // viscosity coefficient
      coeff[0] = DiffCoeff;

      // convection field
      coeff[1] = 1. - y;
      coeff[2] = -2./3. + x;

      // reaction coefficient
      coeff[3] = 1.0;

      // source term
      /// u_t
      rhs  = M_PI*std::cos(M_PI*t) * std::sin(2*M_PI*x) * std::sin(2*M_PI*y);
      /// - DiffCoeff * (u_xx+u_yy)
      rhs += -coeff[0] * exact[3];
      /// b1 * u_x + b2 * u_y
      rhs += coeff[1] * exact[1]
           + coeff[2] * exact[2];
      /// c * u
      rhs += coeff[3] * exact[0];
      coeff[4] = rhs;
    }
  }

} // namespace sin3_0_pi

#endif

