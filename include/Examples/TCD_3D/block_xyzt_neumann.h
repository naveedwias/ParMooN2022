#include "Constants.h"

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = true;
double diffCoeff;

void print_information_about_example()
{
  Output::print("Example: block_xyzt_neumann - solution x*y*z*t in domain of form of block 1x2x3 with Dirichlet and Neumann boundary conditions");
}
  
  
// void get_c_exact_solution_values(double x, double y, double z, double t, double *values);
void get_c_exact_solution_values(double x, double y, double z, double t, double* values)
{
  values[0] = x*y*z*t;
  values[1] = y*z*t;
  values[2] = x*z*t;
  values[3] = x*y*t;
  values[4] = 0;
}


// double get_c_exact_solution_time_derivative_value(double x, double y, double z, double t)
double get_c_exact_solution_time_derivative_value(double x, double y, double z, double )
{
  return x*y*z;
}


// void get_epsilon_values(double x, double y, double z, double t, double* values)
void get_epsilon_values(double , double , double , double , double* values)
{
  values[0] = diffCoeff;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
}


// void get_b1_values(double x, double y, double z, double t, double* values)
void get_b1_values(double , double , double , double , double* values)
{
  values[0] = 1.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// void get_b2_values(double x, double y, double z, double t, double* values)
void get_b2_values(double , double , double , double , double* values)
{
  values[0] = 1.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// void get_b3_values(double x, double y, double z, double t, double* values)
void get_b3_values(double , double , double , double , double* values)
{
  values[0] = 1.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// double get_coeff_c_value(double x, double y, double z, double t)
double get_coeff_c_value(double , double , double , double )
{
  return 1.0;
}


// double get_f_value(double x, double y, double z, double t, const double* solution_values)
double get_f_value(double x, double y, double z, double t, const double* )
{
  double exact[5], b1[5], b2[5], b3[5], epsilon[5];
  get_c_exact_solution_values(x, y, z, t, exact);
  get_b1_values(x, y, z, t, b1);
  get_b2_values(x, y, z, t, b2);
  get_b3_values(x, y, z, t, b3);
  get_epsilon_values(x, y, z, t, epsilon);
  
  return get_c_exact_solution_time_derivative_value(x, y, z, t) - epsilon[0]*exact[4]
         - ( epsilon[1]*exact[1] + epsilon[2]*exact[2] + epsilon[3]*exact[3] )
         + b1[0]*exact[1] + b2[0]*exact[2] + b3[0]*exact[3]
         + get_coeff_c_value(x, y, z, t)*exact[0];
}

//
// boundary conditions
//


// void get_boundary_condition_tcd(int, double x, double y, double z, BoundCond &cond)
void get_boundary_condition(int, double x, double y, double z, BoundCond &cond)
{
  if( (x == 1.0 && y > 0.0 && z > 0.0) || (y == 2.0 && x > 0.0 && z > 0.0) || (z == 3.0 && x > 0.0 && y > 0.0) )
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}


// double get_c_boundary_value(int BdComp, double x, double y, double z, double t)
double get_c_boundary_value(int, double x, double y, double z, double t)
{
  double val[5], epsion_val[5];
  get_c_exact_solution_values(x, y, z, t, val);
  get_epsilon_values(x, y, z, t, epsion_val);
  
  if(x == 1.0 && y > 0.0 && z > 0.0)
    return epsion_val[0]*val[1];
  if(y == 2.0 && x > 0.0 && z > 0.0)
    return epsion_val[0]*val[2];
  if(z == 3.0 && x > 0.0 && y > 0.0)
    return epsion_val[0]*val[3];
    
  return val[0];
}


//
// routines defined solely by means of other routines
//

void get_c_exact_solution_values_no_t(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  get_c_exact_solution_values(x, y, z, t, values);
}


void get_c_initial_values(double x, double y, double z, double *values)
{
  get_c_exact_solution_values(x, y, z, 0.0, values);
}


void get_c_boundary_value_no_t(int BdComp, double x, double y, double z, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  value = get_c_boundary_value(BdComp, x, y, z, t);
}


// void get_c_boundary_time_derivative_no_t(int BdComp, double x, double y, double z, double &value)
void get_c_boundary_time_derivative_no_t(int, double x, double y, double z, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  value = get_c_exact_solution_time_derivative_value(x, y, z, t);
}


void bilinear_coeffs(int n_points, const double* X, const double* Y, const double* Z, const double*const* , double** coeffs)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  double b1_val[5], b2_val[5], b3_val[5];
  double epsilon_val[5];
  
  for(int i=0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    double z = Z[i];
    
    get_b1_values(x, y, z, t, b1_val);
    get_b2_values(x, y, z, t, b2_val);
    get_b3_values(x, y, z, t, b3_val);
    get_epsilon_values(x, y, z, t, epsilon_val);

    coeffs[i][0] = epsilon_val[0];
    coeffs[i][1] = b1_val[0];
    coeffs[i][2] = b2_val[0];
    coeffs[i][3] = b3_val[0];
    coeffs[i][4] = get_coeff_c_value(x, y, z, t);
        
    coeffs[i][5] = get_f_value(x, y, z, t, b1_val);
  }
}
