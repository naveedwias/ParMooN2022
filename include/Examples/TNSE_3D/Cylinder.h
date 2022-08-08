// Navier-Stokes problem
// channel flow in 3D
//
#include <MooNMD_Io.h>



double DIMENSIONLESS_VISCOSITY;
// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  Output::print("Example: Cylinder.h, Re ");
}


void InitialU1(double, double, double, double *values)
{
  
    values[0] = 0;
}


void InitialU2(double, double, double, double *values)
{
  
    values[0] = 0;
}


void InitialU3(double, double, double, double *values)
{

    values[0] = 0;
}


void InitialP(double, double, double, double *values)
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
  cond = DIRICHLET;

  // Neumann boundary condition on outlet  (top)
  double _R_CYLINDER = 2.;
  double _HEIGHT = 10.;
  double r = std::sqrt(x*x+y*y);
  if (std::abs(z-_HEIGHT)<1e-8) {
    if ( std::abs(r - _R_CYLINDER)>1e-8) {
      cond = NEUMANN;
      //Output::print("NEUMANN on top");
    }
  }

  /*
  // Neumann boundary condition on inlet  (bottom)
   if (std::abs(z)<1e-8) {
    if ( std::abs(r2 - _R_CYLINDER*_R_CYLINDER)>1e-8) {
      cond = NEUMANN;
      //Output::print("NEUMANN on bottom");
    }
  }
  */

}


// value of boundary condition
void U1BoundValue(int, double, double, double, double &value)
{
  value = 0;
}


// value of boundary condition
void U2BoundValue(int, double, double, double, double &value)
{
  value = 0;
}


// value of boundary condition
void U3BoundValue(int, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  double _R_CYLINDER = 2.;
  double _UMAX = 1.;
//  double _HEIGHT = 10.;
//  double _DELTA_P = 5.;
//  double _OMEGA = 2*3.1415;

  // set value to 0 for the side walls (Dirichlet)
  // and outlet on the top (Neumann)
  value = 0;

  // parabolic profile at inlet (Dirichlet)
  double r2 = x*x+y*y;
  double RC2 = _R_CYLINDER*_R_CYLINDER;
  if (std::abs(z)<1e-8) {
  value = _UMAX*(1-r2/RC2)*t;
  }

  /*
  // plug profile
  if (abs(r2-_R_CYLINDER*_R_CYLINDER)<1e-8) {
    //
    value = 0.;
  } else {
    value = _UMAX*t;
  }
  */
}




// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  double nu = DIMENSIONLESS_VISCOSITY;

//  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  //Output::print("bulk velocity: mean ", u1, " sim ", u2);
  
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f3
  }
}
