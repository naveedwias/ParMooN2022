#include "Constants.h"

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = true;
double diffCoeff;

void print_information_about_example()
{
  Output::print("Example: rectangle_xyt_neumann - solution x*y*t with Dirichlet (parts 0, 3) and Neumann (parts 1, 2) boundary conditions and domain 3x1");
}
  

// void get_c_exact_solution_values(double x, double y, double t, double* values)
void get_c_exact_solution_values(double x, double y, double t, double* values)
{
  values[0] = x*y*t;
  values[1] = y*t;
  values[2] = x*t;
  values[3] = 0;
}


// double get_c_exact_solution_time_derivative_value(double x, double y, double t)
double get_c_exact_solution_time_derivative_value(double x, double y, double )
{
  return x*y;
}


// void get_epsilon_values(double x, double y, double t, double* values)
void get_epsilon_values(double , double , double , double* values)
{
  values[0] = diffCoeff;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// void get_b1_values(double x, double y, double t, double* values)
void get_b1_values(double , double , double , double* values)
{
  values[0] = 1.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// void get_b2_values(double x, double y, double t, double* values)
void get_b2_values(double , double , double , double* values)
{
  values[0] = 1.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}


// double get_coeff_c_value(double x, double y, double t)
double get_coeff_c_value(double , double , double )
{
  return 1.0;
}


// double get_f_value(double x, double y, double t, const double* solution_values)
double get_f_value(double x, double y, double t, const double* )
{
  double exact[4], b1[4], b2[4], epsilon[4];
  get_c_exact_solution_values(x, y, t, exact);
  get_b1_values(x, y, t, b1);
  get_b2_values(x, y, t, b2);
  get_epsilon_values(x, y, t, epsilon);
  
  return get_c_exact_solution_time_derivative_value(x, y, t) - epsilon[0]*exact[3]
         - ( epsilon[1]*exact[1] + epsilon[2]*exact[2] )
         + b1[0]*exact[1] + b2[0]*exact[2] + get_coeff_c_value(x, y, t)*exact[0];
}


//
// boundary conditions
//


void get_boundary_condition(int BdComp, double, BoundCond &cond)
{
  if(BdComp == 0 || BdComp == 3)
    cond = DIRICHLET;
  else // BdComp = 0, 2
    cond = NEUMANN;
}


void get_xy_of_boundary_parametrization(int BdComp, double Param,
                                        double& x, double& y)
{
  switch(BdComp)
  {
    case 0:
      x = 3.0*Param;
      y = 0.0;
    break;
    case 1:
      x = 3.0;
      y = Param;
    break;
    case 2:
      x = 3.0*(1.0 - Param);
      y = 1.0;
    break;
    case 3:
      x = 0.0;
      y = 1-Param;
    break;
    default:
      ErrThrow("wrong boundary component ", BdComp);
  }
}


double get_c_boundary_value(int BdComp, double Param, double t)
{
  double x, y;
  get_xy_of_boundary_parametrization(BdComp, Param, x, y);
  
  double val[4];
  get_c_exact_solution_values(x, y, t, val);
  
  if(BdComp == 0 || BdComp == 3)
  {
    return val[0];
  }
  else
  {
    double epsion_val[4];
    get_epsilon_values(x, y, t, epsion_val);
    
    if(BdComp == 2)
      return epsion_val[0]*val[2];
    else // BdComp = 1
      return epsion_val[0]*val[1];
  }
}


//
// routines defined solely by means of other routines
//

void get_c_exact_solution_values_no_t(double x, double y, double* values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  get_c_exact_solution_values(x, y, t, values);
}


void get_c_initial_values(double x, double y, double* values)
{
  get_c_exact_solution_values(x, y, 0.0, values);
}


void get_c_boundary_value_no_t(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  value = get_c_boundary_value(BdComp, Param, t);
}


void get_c_boundary_time_derivative_no_t(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  double x, y;
  get_xy_of_boundary_parametrization(BdComp, Param, x, y);
  
  value = get_c_exact_solution_time_derivative_value(x, y, t);
}


// void bilinear_coeffs(int n_points, const double* X, const double* Y, const double*const* par_val, double** coeffs)
void bilinear_coeffs(int n_points, const double* X, const double* Y, const double*const* , double** coeffs)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  double b1_val[4], b2_val[4];
  double epsilon_val[4];
  
  for(int i=0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    
    get_b1_values(x, y, t, b1_val);
    get_b2_values(x, y, t, b2_val);
    get_epsilon_values(x, y, t, epsilon_val);

    coeffs[i][0] = epsilon_val[0];
    coeffs[i][1] = b1_val[0];
    coeffs[i][2] = b2_val[0];
    coeffs[i][3] = get_coeff_c_value(x, y, t);
        
    coeffs[i][4] = get_f_value(x, y, t, b1_val);
  }
}
