// Time Convection Diffusion 2D problem for a traveling wave in a unit square
// domain [0,1]x[0,1], used for the paper: SUPG reduced order models for
// convection-dominated convection–diffusion–reaction equations, S. Giere,
// T. Iliescu, V. John and D. Wells, 2015
// The thickness of the internal layer is sqrt(DiffCoeff)

#ifndef __TRAVELING_WAVE__
#define __TRAVELING_WAVE__


namespace traveling_wave
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
      Output::print("\nExample: traveling_wave_sqrt.h");
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
    double eps = std::sqrt(DiffCoeff); 
    double t = TDatabase::TimeDB->CURRENTTIME;
    double ch = std::cosh((-0.5 - t + x + y) / eps);
    double th = std::tanh((-0.5 - t + x + y) / eps);

    // u
    values[0]  = 0.5 * std::sin(M_PI*x) * std::sin(M_PI*y) * (1 + th);

    // u_x
    values[1]  = ( 0.5 * std::sin(M_PI*x) * std::sin(M_PI*y) ) / (eps * ch*ch);
    values[1] += 0.5 * M_PI * std::cos(M_PI*x) * std::sin(M_PI*y) * (1 + th);

    // u_y 
    values[2]  = 0.5 * std::sin(M_PI*x) * std::sin(M_PI*y) / (eps * ch*ch);
    values[2] += 0.5 * M_PI * std::cos(M_PI*y) * std::sin(M_PI*x) * (1 + th);

    // u_xx + u_yy
    values[3]  = - M_PI*M_PI * std::sin(M_PI*x) * std::sin(M_PI*y) * (1 + th);
    values[3] += M_PI * std::sin(M_PI*(x+y)) / (eps * ch*ch);
    values[3] -= 2 * std::sin(M_PI*x) * sin(M_PI*y) * th / ( eps*eps * ch*ch); 
  }

  // ===========================================================================
  // initial conditon
  // ===========================================================================
  void InitialCondition(double x,  double y, double* values)
  {
    double eps = std::sqrt(DiffCoeff); 
    double t = TDatabase::TimeDB->CURRENTTIME;

    double th = std::tanh((-0.5 - t + x + y) / eps);

    values[0]  = 0.5 * std::sin(M_PI*x) * std::sin(M_PI*y) * (1 + th);
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
    double rhs;
    double val[4];

    for(int i=0 ; i<n_points ; i++)
    {
      coeff = coeffs[i];
      double x = X[i];
      double y = Y[i];
      double eps = std::sqrt(DiffCoeff);

      // viscosity coefficient
      coeff[0] = DiffCoeff;

      // convection: b1 
      coeff[1] = std::cos(M_PI/3.);
      // convection: b2
      coeff[2] = std::sin(M_PI/3.);

      // reaction coefficient
      coeff[3] = 1.0;

      // source term
      Exact(x, y, val);
      /// u_t
      double ch = std::cosh((-0.5 - t + x + y) / eps);
      rhs  = - 0.5 * sin(M_PI*x) * sin(M_PI*y) / (eps * ch*ch);
      /// c * u
      rhs += coeff[3] * val[0];
      /// b1 * u_x + b2 * u_y
      rhs += coeff[1] * val[1] + coeff[2] * val[2];
      /// - DiffCoeff * (u_xx+u_yy)
      rhs += -coeff[0] *  val[3];
      coeff[4] = rhs;
    }
  }

} // namespace traveling_wave

#endif
