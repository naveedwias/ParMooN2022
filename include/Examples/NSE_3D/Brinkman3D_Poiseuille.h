/*********************************************************
 * Analytical solution for the Brinkman 3D problem
 * - analogous to a Poiseuille flow
 * Stokes regime (sigma = 0):
   ux = Poiseuille solution (unitary pressure drop)
 * Darcy regime (effective_viscosity = 0):
   ux = 1/sigma
 * Brinkman regime:
   ux = exponential boundary layer
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
double sigma = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

//@attention !! we assume a cylinder of radius 2 and heigth 10 !!
double _R_CYLINDER = 2.;
double _HEIGHT = 10.;
double _DELTA_P = 2.;
void ExampleFile()
{
  Output::print<1>("Example: Brinkman3D_Poiseuille.h");
}

// ========================================================================
// exact solution
// Note: computed assuming a pressure drop = 2, and mean pressure = 0
//       radius = 1 and height = 10
// ========================================================================


void ExactU1(double, double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.; // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU2(double, double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.; // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU3(double x, double y, double, double *values)
{
  // ratio between Stokes and Darcy terms
  double t = std::abs(std::sqrt( effective_viscosity/sigma));


  // DeltaP = 4 mu  L Umax / R ^2
  // Umax = R^2 * DeltaP / (4 mu L) = 1/(10 mu)
  double _UMAX = _DELTA_P * _R_CYLINDER*_R_CYLINDER /( 4. * _HEIGHT *  effective_viscosity);
  double _R2 = _R_CYLINDER*_R_CYLINDER;

  Output::print("sigma = ", sigma);
  if ( std::abs(sigma) < 1e-10 )
  {
    // Stokes regime
    values[0] = _UMAX * (1 - (x*x + y*y)/_R2 ); 
    values[1] = _UMAX * (- 2. * x /_R2 );
    values[2] = _UMAX * (- 2. * y /_R2);
    values[3] = 0.;
    values[4] = -4 * _UMAX/_R2;
  }
  else if  (t == 0)
  {
    // Darcy regime
    values[0] = _DELTA_P/_HEIGHT * 1./sigma;  
    values[1] = 0.;
    values[2] = 0.;
    values[3] = 0.;
    values[4] = 0.;
  }
  else
  {
    ///@todo check the solution for  t>0
    values[0] = (1./sigma) * (1+std::exp(2/t)-std::exp( (2- std::sqrt(x*x+y*y)  )/t) - std::exp( std::sqrt(x*x+y*y)/t)) / (1+std::exp(2/t));
        //(1./sigma) * (1+std::exp(1/t)-std::exp(1- ( (1/4)*x+(1/2)  )/t) - std::exp( ((1/4)*x+(1/2) )/t)) / (1+std::exp(1/t))+
                 //(1./sigma) * (1+std::exp(1/t)-std::exp(1- (  (1/4)*y+(1/2) )/t) - std::exp( ( (1/4)*y+(1/2) )/t)) / (1+std::exp(1/t));
    //(1./sigma) * (1+std::exp(1/t)-std::exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - std::exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (1+std::exp(1/t));
    values[1] = (1./sigma) * (std::exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - std::exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (4*t*(1+std::exp(1/t)));
    values[2] = (1./sigma) * (std::exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - std::exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (4*t*(1+std::exp(1/t)));
    values[3] = 0.;
    values[4] = 2* (1./sigma) * (std::exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - std::exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (16*t*t*(1+std::exp(1/t)));
  }
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */
void ExactP(double, double, double z, double *values)
{
  values[0] = 1-(_DELTA_P/_HEIGHT)*z;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = -(_DELTA_P/_HEIGHT);
  values[4] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================


// kind of boundary condition (for FE space needed)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{

  double _R_CYLINDER = 2.;
  double _HEIGHT = 10.;
  double r = std::sqrt(x*x+y*y);
  
  if (nitsche_id.size()==1)
  {
    if (nitsche_id[0] == 2)
    {
      cond = DIRICHLET_WEAK;
    }
    else
    {
      Output::print(" Warning: Nitsche BC on component ", nitsche_id[0],
		    " are not what I expect. Setting cond = DIRICHLET");
      cond = DIRICHLET;
    }
    
  }
  else
  {
    cond = DIRICHLET;
  }

  if (neumann_id.size()==1)
    Output::print(" Warning: neumann_id.size = 1. This example should not support this case");
    
  if (neumann_id.size()==2)
  {
    // set Neumann BC on top and bottom
    if ((std::abs(z-_HEIGHT)<1e-14) || (std::abs(z-0)<1e-14) )
    {
      if( std::abs(r - _R_CYLINDER)>1e-14  )
      {
	cond = NEUMANN;
      }
    }
  }
}

void U1BoundValue(int, double, double, double, double &value) // (int BdComp, double Param, double &value)
{
    value = 0.;

}

void U2BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

void U3BoundValue(int, double x, double y, double, double &value)
{
  if (neumann_id.size())
  {
    value = 0.; // we only need Dirichlet on the wall, i.e., u = 0
  }
  else
  {
    // all dirichlet BC
    double _UMAX = _DELTA_P * _R_CYLINDER*_R_CYLINDER /( 4. * _HEIGHT *  effective_viscosity);
    double _R2 = _R_CYLINDER*_R_CYLINDER;
    value = _UMAX * (1 - (x*x + y*y)/_R2 );
  }
}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][5] = sigma;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][4] = 0; 
  }
}


