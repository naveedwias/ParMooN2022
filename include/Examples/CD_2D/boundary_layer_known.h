// ======================================================================
// boundary layer
// John, Knobloch, Savescu, CMAME 2011
// ======================================================================
#define __cd2d_hump__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: boundary_layer_known.h") ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;
}

// exact solution
void Exact(double x, double y, double *values)
{
    //intermediate calculation variables
    double exp_x, exp_y, exp_x_y, u_xx, u_yy;
    
    //diffusion coeffecient
    double eps;
    eps = 1. / DIFFUSION;
     
    //intermediate calculations
    exp_x=std::exp(2*(x-1)*eps);
    exp_y=std::exp(3*(y-1)*eps);
    exp_x_y=std::exp((2*(x-1)+3*(y-1))*eps);
    
    values[0] = x*y*y-y*y*exp_x-x*exp_y+exp_x_y;  //original solution
    values[1] = y*y-y*y*2*eps*exp_x-exp_y+ 2*eps*exp_x_y; //x-derivative
    values[2] = 2*x*y-2*y*exp_x-x*3*eps*exp_y+ 3*eps*exp_x_y; //y-derivative
    
    u_yy = 2*x- 2*exp_x- 9*x*eps*eps*exp_y+ 9*eps*eps*exp_x_y;
    u_xx = -y*y*4*eps*eps*exp_x+ 4*eps*eps*exp_x_y;
    
    values[3] = u_xx+u_yy; //laplacian
 }

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

//to get the boundary lengths of the example
void transform(const int BdComp, const double Param, double& x, double& y)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      break;
    case 1:
      x = 1.;
      y = Param;
      break;
    case 2:
      x = 1.-Param;
      y = 1.;
      break;
    case 3:
      x = 0.;
      y = 1.-Param;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  //coordinates for the part of boundary
  double x,y;
  transform(BdComp, Param, x, y);
  
  //type of boundary condition (DIRICHLET or NEUMANN)
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  
  //exact solution at the point (x,y), required for BC
  double u[4];
  Exact(x, y, u);
  
  if(cond == DIRICHLET)
  {
    value = u[0];
  }
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
    Exact(x,y,values);
}


void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = DIFFUSION;
  double *coeff;
  double exact[4]; 
  
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //param = parameters[i];
  
    coeff[0] = eps;
    coeff[1] = 2; //x-coordiante of convection b
    coeff[2] = 3; //y-coordiante of convection b
    coeff[3] = 1; //reaction coeffecient c
    
    Exact(X[i], Y[i], exact);
    
    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
  }
}
