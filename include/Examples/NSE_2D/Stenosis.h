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


void ExampleFile()
{
  Output::print<1>("Example: Stenosis.h");
}

// ========================================================================
// exact solution (in this case, this is only an approximated solution
// ========================================================================
void ExactU1(double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
  values[3] = 0.; 
}

void ExactU2(double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
  values[3] = 0.; 
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */

void ExactP(double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
  values[3] = 0.; 
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (i == (int)neumann_id[j])
    {
      cond = NEUMANN;
      cout << " NEUMANN ON BD " << i << endl;
      return;
    }
  }

  // set Nitsche BC
  for (unsigned int j = 0; j < nitsche_id.size(); j++)
  {
    if (i == (int)nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;
      
      return;
    }
  }
}

void U1BoundValue(int, double, double &value)
{

   value = 0;
   
}



void U2BoundValue(int BdComp, double t, double &value)
{
  double UMAX = 10.;
  if (BdComp == 21)
  {
    value = UMAX * 4*t*(1-t);
  } else {
    value = 0;
  }

}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double*, const double*,
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
    
  }
}


