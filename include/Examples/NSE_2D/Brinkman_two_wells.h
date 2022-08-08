/*********************************************************
 * Model for geothermal flow (geothermal plant with 2 wells)
 * 
 * Consider a rectangular domain (L1 x L2) with two small circles
 * of radius r_well at positions (xi, yi) and (xe, ye) (i=injection, e=extraction)
 * Both circles are physical boundaries but not resolved by the mesh
 * Instead, they are modeled as singular source and sink
 * Parameters/BC
 * Q = incoming/outgoing flow at wells
 * |u| at wells = Q/(2*Pi*r_well)
 * p(outer boundary) = 0
 *
 * For this problem we do not have an exact solution
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1.;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

// we consider a domain [0,10000]x[0,6000] (m)
double Q = 150./3600.; // 150m^3/h
double r_well = 0.2; //m 
double xi = 4500.;
double yi = 3000.;
double xe = 5500.;
double ye = 3000.;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_two_wells.h");
}

// ========================================================================
// exact solution (in this case, this is only an approximated solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double r2 = (x-xi)*(x-xi) + (y-yi)*(y-yi);

  values[0] = Q/(2*M_PI) * (x-xi)/r2;        
  values[1] = Q/(2*M_PI) * ( (-(x-xi)*(x-xi)+(y-yi)*(y-yi))/(r2*r2) );
  values[2] = Q/(2*M_PI) * ( -2*(x-xi)*(y-yi)/(r2*r2) );

  r2 = (x-xe)*(x-xe) + (y-ye)*(y-ye);
  values[0] -= Q/(2*M_PI) * (x-xe)/r2;        
  values[1] -= Q/(2*M_PI) * ( (-(x-xe)*(x-xe)+(y-ye)*(y-ye))/(r2*r2) );
  values[2] -= Q/(2*M_PI) * ( -2*(x-xe)*(y-ye)/(r2*r2) );
  
  values[3] = 0.; 
}

void ExactU2(double x, double y, double *values)
{
  double r2 = (x-xi)*(x-xi) + (y-yi)*(y-yi);

  values[0] = Q/(2*M_PI) * (y-yi)/r2;
  values[1] = Q/(2*M_PI) * (-2*(y-yi)*(x-xi)/(r2*r2));
  values[2] = Q/(2*M_PI) * (-(y-yi)*(y-yi) + (x-xi)*(x-xi))/( r2*r2);

  r2 = (x-xe)*(x-xe) + (y-ye)*(y-ye);
  values[0] -= Q/(2*M_PI) * (y-ye)/r2;
  values[1] -= Q/(2*M_PI) * (-2*(y-ye)*(x-xe)/(r2*r2));
  values[2] -= Q/(2*M_PI) * (-(y-ye)*(y-ye) + (x-xe)*(x-xe))/( r2*r2);
  
  values[3] = 0.; 
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */

void ExactP(double x, double y, double *values)
{
  double r2 = (x-xi)*(x-xi) + (y-yi)*(y-yi);
  double r_1 = 6000.; // max(L1-xi,L2-yi)/2

  values[0] = -sigma * Q/(2*M_PI) * 0.5 * std::log( r2/(r_1*r_1) );     
  values[1] = -sigma * Q/(2*M_PI) * (x-xi)/r2;
  values[2] = -sigma * Q/(2*M_PI) * (y-yi)/r2;

  r2 = (x-xe)*(x-xe) + (y-ye)*(y-ye);
  values[0] -= -sigma * Q/(2*M_PI) * 0.5 * std::log( r2/(r_1*r_1) );     
  values[1] -= -sigma * Q/(2*M_PI) * (x-xe)/r2;
  values[2] -= -sigma * Q/(2*M_PI) * (y-ye)/r2;

  values[3] = 0.;
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
  case 1:
  case 2:
  case 3:
    {
      value = 0.;  
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
  case 1:
  case 2:
  case 3:
    {
      value = 0.;  
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

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = 0.;

    
    std::vector<double> singular_x,singular_y,singular_sign;
    // source in (xi,yi)
    singular_x.push_back(xi);
    singular_y.push_back(yi);
    singular_sign.push_back(1.);
    // sink in (xe,ye)    
    singular_x.push_back(xe);
    singular_y.push_back(ye);
    singular_sign.push_back(-1.);

    // parameter for approximated delta function
    double epsilon = 25*r_well; 

    // approximated delta function centered in (x[m],y[m])
    for (unsigned int m=0; m<singular_x.size(); m++) {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
      
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      bool at_source = (x_distance_to_source < epsilon*epsilon) *
	(y_distance_to_source < epsilon*epsilon);
      
      if(at_source)
	{
	  Output::print<4>(" adding a singular source - point ", m);
	  double delta_h = std::cos(M_PI*(x[i] - x_center_source)/epsilon) + 1;
	  delta_h *= std::cos(M_PI*(y[i] - y_center_source)/epsilon) + 1;
	  delta_h /= 4.*epsilon*epsilon;
	  //g(x,y):  RHS for mass conservation equation
	  coeffs[i][3] += singular_sign[m] * delta_h * Q;

	  // other option: add a constant value so that int delta = Q
	  //coeffs[i][3] += singular_sign[m] * Q/(4.*epsilon*epsilon);
	    
	}
      
    }

    
  }
}


