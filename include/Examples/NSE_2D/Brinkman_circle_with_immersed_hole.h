/*********************************************************
 * Analytical solution for the Brinkman 2D problem
 * 
 * Consider a circular domain with radius r_1 and center point (0,0).
 * Inscribed to that circle is another on with smaller radius r_0 and the same center.
 * The inner circle is not considered as boundary, but the boundary condition is taken
 * into account using a singular term in the mass conservation equation (immersed boundary method)
 * (The appropriate mesh is a circle)
 * The radial velocity is at the inner boundary is fixed by the strength
 * of the singular force term. The exact solution  (outside of the inner ring) should be
 * u(r=r0) = u0
 * p(r=R1) = 0
 * and then the solution is 
 * u(r) = u0 * r_0/r_1
 * p(r) = -sigma * u0 * r_0 * std::log(r/r_1) 
 *
 * r := std::sqrt(x^2 + y^2)
 * Hence,
 * u_1 = u0 * r_0 * x/(x^2 + y^2);
 * u_2 = u0 * r_0 * y/(x^2 + y^2);
 * p = -sigma * u0 * r_0 * 0.5 * std::log( (x^2 + y^2)/r_1^2 );
 * du_1/dx = u0 * r_0 * ( (-x^2 +y^2)/(x^2 + y^2)^2 );
 * d^2u_1/dx^2 = u0 * r_0 * ( 2 * x * (x^2-3y^2)/((x^2+y^2)^3))
 * du_1/dy = u0 * r_0 * ( -2xy/((x^2+y^2)^2) ); 
 * d^2u_1/dy^2 = u0 * r_0 * ( 2*x*(3y^2-x^2)/((x^2+y^2)^3) )
 * 
 * Delta (u_1,u_2) = 0
 * (f_1,f_2) = (0,0);
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1.;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;


double u0 = 1e-5;
double r_1 = 1000;
double r_0 = 10;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_circle_with_immersed_hole.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double r2 = x*x + y*y;
  if (r2<=r_0) {
    values[0]=0.;
    values[1]=0.;
    values[2]=0.;
    values[3]=0.;
  } else {
    values[0] = u0 * r_0 * x/r2;        
    values[1] = u0 * r_0 * ( (-x*x+y*y)/(r2*r2) );
    values[2] = u0 * r_0 * ( -2*x*y/(r2*r2) );    
    values[3] = 0.;
  }
}

void ExactU2(double x, double y, double *values)
{
  double r2 = x*x + y*y;
  if (r2<=r_0) {
    values[0]=0.;
    values[1]=0.;
    values[2]=0.;
    values[3]=0.;
  } else {
    values[0] = u0 * r_0 * y/r2;
    values[1] = u0 * r_0 * (-2*y*x/(r2*r2));
    values[2] = u0 * r_0 * (-y*y + x*x)/( r2*r2);
    values[3] = 0.;
  }
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */

void ExactP(double x, double y, double *values)
{
  double r2 = x*x + y*y;
  if (r2<=r_0) {
    values[0]=0.;
    values[1]=0.;
    values[2]=0.;
    values[3]=0.;
  } else {
    values[0] = -sigma * u0 * r_0 * 0.5 * std::log( r2/(r_1*r_1) );     
    values[1] = -sigma * u0 * r_0 * x/r2;
    values[2] = -sigma * u0 * r_0 * y/r2;
    values[3] = 0.;
  }
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = NEUMANN; 
}

void U1BoundValue(int BdComp, double, double &value)
{
  
  switch(BdComp)
  {
  case 0:
    {
      value = 0.;  //p(r=R1) = 0
      break;
    }
  default: cout << "No boundary component with this number." << endl;
    break;
  }
}



void U2BoundValue(int BdComp, double, double &value)
{
  switch(BdComp)
  {
  case 0:
    {
      value = 0.;  //p(r=R1) = 0
      break;
    }
  default: cout << "No boundary component with this number." << endl;
    break;
  }
}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double *x, const double *y,
               const double *const*, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][4] = sigma;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;

    // source data
    
    double epsilon = 5*r_0; //0.05;
    
    double Q_in = u0*2*M_PI*r_0;

    
    // approximated delta function centered in (0,0)
    double x_center_source = 0.;
    double y_center_source = 0.;
    double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
    double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
    bool at_source = (x_distance_to_source < epsilon*epsilon) *
      (y_distance_to_source < epsilon*epsilon);

    coeffs[i][3] = 0.;
    if(at_source)
    {
      double delta_h = std::cos(M_PI*(x[i] - x_center_source)/epsilon) + 1;
      delta_h *= std::cos(M_PI*(y[i] - y_center_source)/epsilon) + 1;
      delta_h /= 4.*epsilon*epsilon;
      //g(x,y):  RHS for mass conservation equation
      coeffs[i][3] = delta_h * Q_in;
    }
  }
}


