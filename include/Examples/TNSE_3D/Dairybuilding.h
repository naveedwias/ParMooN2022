// Navier-Stokes problem
// Dairybuilding in 3D
//

double DIMENSIONLESS_VISCOSITY;
// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  Output::print("Example: Dairybuilding.h");
}

void InitialU1(double , double , double z, double *values)
{
  // interpolate experimental data (50m downstream)
  std::vector<double> zval(19),uval(19);
  zval[0] = 0.; uval[0] = 0.;
  zval[1] = 1.e-2; uval[1] = 3.52;
  zval[2] = 2.e-2; uval[2] = 3.77;
  zval[3] = 3.e-2; uval[3] = 3.94;
  zval[4] = 4.e-2; uval[4] = 4.01;
  zval[5] = 5.e-2; uval[5] = 4.10;
  zval[6] = 7.5e-2; uval[6] = 4.24;
  zval[7] = 10e-2; uval[7] = 4.41;
  zval[8] = 12.5e-2; uval[8] = 4.55;
  zval[9] = 15.e-2; uval[9] = 4.67;
  zval[10] = 20.e-2; uval[10] = 4.82;
  zval[11] = 30.e-2; uval[11] = 5.18;
  zval[12] = 40.e-2; uval[12] = 5.46;
  zval[13] = 50.e-2; uval[13] = 5.71;
  zval[14] = 60.e-2; uval[14] = 5.96;
  zval[15] = 70.e-2; uval[15] = 6.23;
  zval[16] = 80.e-2; uval[16] = 6.38;
  zval[17] = 90.e-2; uval[17] = 6.57;
  zval[18] = 100.e-2; uval[18] = 6.81;

  int i_range = -1;
  for (unsigned int i=0; i<zval.size()-1; i++)
  {
    if ( (z>=zval[i]) && (z<zval[i+1]) )
      i_range = i;
  }
  if (i_range>=0)
  {
    values[0] = uval[i_range] +
      (uval[i_range+1]-uval[i_range])*
      (z - zval[i_range])/(zval[i_range+1] - zval[i_range]);
  }
  else
    values[0] = uval[18]; // case z = 100;
  
  // separate treatment of the first interval to mimick the boundary layer
  if ( (z>zval[0]) && (z<zval[1]) )
      values[0] = uval[1];

  // if P5 positive, start with 0, else developed flow
  double t;
  double t_start = TDatabase::ParamDB->P5;
  double t_end = TDatabase::ParamDB->P6;
  if (t_start >= 0)
    t = 0;
  else if (t_end <= 0)
    t = 1;
  else // t_start negative and t_end positive, smoothed initial solution
    t = -t_start/(t_end-t_start);  
  
  values[0] *= t;
  //values[0] = 0;
  //values[0] *= 5e-2;
}

void InitialU2(double , double , double , double *values)
{
  
    values[0] = 0;
}

void InitialU3(double , double , double , double *values)
{

    values[0] = 0;
}

void InitialP(double , double ,  double , double *values)
{
  values[0] = 0;
}

auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP = unknown_solution_3d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  double height = 1, outflow = 2.5, inflow = -0.5, width = 0.342;
  double eps = 1.e-6;

  // DIRICHLET everywhere except on...
  cond = DIRICHLET;

  // ... and neumann for outflow
  if (std::abs(x-outflow) < eps )
    cond = NEUMANN; 
  
  // ... both sides and top surface (free slip bc)...
  if ( x > inflow + eps && x < outflow-eps)
  {
    if (std::abs(y) < eps || std::abs(y-width) < eps || std::abs(z-height)<eps)
    {
//       cond = DIRICHLET;
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 0;
      TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM = 0; //0=strong penetr.
      TDatabase::ParamDB->FRICTION_TYPE = 1;
      TDatabase::ParamDB->FRICTION_CONSTANT = 0;
      TDatabase::ParamDB->FRICTION_POWER = 1;
      TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e20;
      TDatabase::ParamDB->PENETRATION_POWER= -4;
    }
  }

  
  // set lower boundary extra
  if (std::abs(z) < eps)
    cond = DIRICHLET;
  if ((std::abs(x) < eps)  && (std::abs(z) <= 0.0048+eps))
    cond = DIRICHLET;
  if ((x>=-eps) && (x<=0.032+eps)  && (std::abs(z-0.0048) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.0032) < eps)  && (z <= 0.0048+eps))
    cond = DIRICHLET;
  if (( x>= 0.032 -eps) && (x<=0.1485+eps)  && (std::abs(z-0.0016) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.1485) < eps)  && (z <= 0.0016+eps))
    cond = DIRICHLET;
  if (( x>= 0.1485 -eps) && (x<=0.1503+eps)  && (std::abs(z-0.008) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.1503) < eps)  && (z <= 0.008+eps))
    cond = DIRICHLET;
  if (( x>= 0.1503 -eps) && (x<=0.1917+eps)  && (std::abs(z-0.0034) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.1917) < eps)  && (z <= 0.008+eps))
    cond = DIRICHLET;
  if (( x>= 0.1917 -eps) && (x<=0.1935+eps)  && (std::abs(z-0.008) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.1935) < eps)  && (z <= 0.008+eps))
    cond = DIRICHLET;
  if (( x>= 0.1935 -eps) && (x<=0.3385+eps)  && (std::abs(z-0.0016) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.3385) < eps)  && (z <= 0.008+eps))
    cond = DIRICHLET;
  if (( x>= 0.3385 -eps) && (x<=0.342+eps)  && (std::abs(z-0.008) <= eps))
    cond = DIRICHLET;
  if ((std::abs(x-0.342) < eps)  && (z <= 0.008+eps))
    cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double x, double, double z, double &value)
{
  double inflow = -0.5;
  //double noise = TDatabase::ParamDB->P8/100.0;
  // Homogeneous everywhere except on...
  value = 0;

  double t = TDatabase::TimeDB->CURRENTTIME;
  double t_start = TDatabase::ParamDB->P5;
  double t_end = TDatabase::ParamDB->P6;
  double eps = 1.e-6;

  if (t <= t_start)
    t = 0;
  else if (t >= t_start && t <= t_end)
    t = (t-t_start)/(t_end-t_start);
  else // (if t >= t_end_force)
    t = 1;

  // interpolate experimental data (50m downstream)
  std::vector<double> zval(19),uval(19);
  zval[0] = 0.; uval[0] = 0.;
  zval[1] = 1.e-2; uval[1] = 3.52;
  zval[2] = 2.e-2; uval[2] = 3.77;
  zval[3] = 3.e-2; uval[3] = 3.94;
  zval[4] = 4.e-2; uval[4] = 4.01;
  zval[5] = 5.e-2; uval[5] = 4.10;
  zval[6] = 7.5e-2; uval[6] = 4.24;
  zval[7] = 10.e-2; uval[7] = 4.41;
  zval[8] = 12.5e-2; uval[8] = 4.55;
  zval[9] = 15.e-2; uval[9] = 4.67;
  zval[10] = 20.e-2; uval[10] = 4.82;
  zval[11] = 30.e-2; uval[11] = 5.18;
  zval[12] = 40.e-2; uval[12] = 5.46;
  zval[13] = 50.e-2; uval[13] = 5.71;
  zval[14] = 60.e-2; uval[14] = 5.96;
  zval[15] = 70.e-2; uval[15] = 6.23;
  zval[16] = 80.e-2; uval[16] = 6.38;
  zval[17] = 90.e-2; uval[17] = 6.57;
  zval[18] = 100.e-2; uval[18] = 6.81;

  // ... Inflow in the surface 0yz ...
  if ( x < inflow + eps ) // || std::abs(x-outflow) < eps ) // be careful that at it is 0 at y=0
  {
    int i_range = -1;
    for (unsigned int i=0; i<zval.size()-1; i++)
    {
      if ( (z>=zval[i]) && (z<zval[i+1]) )
        i_range = i;
    }
    if (i_range>=0)
    {
      value = uval[i_range] +
        (uval[i_range+1]-uval[i_range])*
        (z - zval[i_range])/(zval[i_range+1] - zval[i_range]);
    }
    else
      value = uval[18]; // case z = 100;

    if ( (z>zval[0]) && (z<zval[1]) )
      value = uval[1];
    
    
    //if (z<0.2)
    //{
    //    value +=  noise * ((double)rand()/RAND_MAX-0.5);
    //}
    value *= t;
  }
  else
    value =0;
  
  
  // set lower boundary extra
  /*if (std::abs(z) < eps)
    value = 0;
  if ((std::abs(x) < eps)  && (std::abs(z) <= 0.0048+eps))
    value = 0;
  if ((x>=-eps) && (x<=0.032+eps)  && (std::abs(z-0.0048) <= eps))
    value = 0;
  if ((std::abs(x-0.0032) < eps)  && (std::abs(z) <= 0.0048+eps))
    value = 0;
  if (( x>= 0.032 -eps) && (x<=0.1485+eps)  && (std::abs(z-0.0016) <= eps))
    value = 0;
  if ((std::abs(x-0.1485) < eps)  && (std::abs(z) <= 0.0016+eps))
    value = 0;
  if (( x>= 0.1485 -eps) && (x<=0.1503+eps)  && (std::abs(z-0.008) <= eps))
    value = 0;
  if ((std::abs(x-0.1503) < eps)  && (std::abs(z) <= 0.008+eps))
    value = 0;
  if (( x>= 0.1503 -eps) && (x<=0.1917+eps)  && (std::abs(z-0.0034) <= eps))
    value = 0;
  if ((std::abs(x-0.1917) < eps)  && (std::abs(z) <= 0.008+eps))
    value = 0;
  if (( x>= 0.1917 -eps) && (x<=0.1935+eps)  && (std::abs(z-0.008) <= eps))
    value = 0;
  if ((std::abs(x-0.1935) < eps)  && (std::abs(z) <= 0.008+eps))
    value = 0;
  if (( x>= 0.1935 -eps) && (x<=0.3385+eps)  && (std::abs(z-0.0016) <= eps))
    value = 0;
  if ((std::abs(x-0.3385) < eps)  && (std::abs(z) <= 0.008+eps))
    value = 0;
  if (( x>= 0.3385 -eps) && (x<=0.342+eps)  && (std::abs(z-0.008) <= eps))
    value = 0;
  if ((std::abs(x-0.342) < eps)  && (std::abs(z) <= 0.008+eps))
    value = 0;*/
  
  //value *= 5e-2;
}

// value of boundary condition
void U2BoundValue(int, double , double , double , double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  double nu = DIMENSIONLESS_VISCOSITY;

  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f3
    coeffs[i][4] = 0; // g
  }
}

