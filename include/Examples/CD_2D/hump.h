// ======================================================================
// circular inner layer
// John, Maubach, Tobiska, Num. Math. 1997
// ======================================================================
#define __cd2d_hump__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: hump.h") ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;
}

// exact solution
void Exact(double x, double y, double *values)
{
    //intermediate calculation variables
    double rxy,cxy,hxy,g1xy,g1xy_x,g1xy_xx,g1xy_y,g1xy_yy;
    double g2xy,g2xy_x,g2xy_xx,g2xy_y,g2xy_yy,gxy_xx,gxy_yy;
    
    //diffusion coeffecient
    double eps;
    eps = 1. / DIFFUSION;
     
    //intermediate calculations
    rxy =  (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
    cxy =  2.0 * std::sqrt(eps)*(0.25*0.25 - rxy);
    hxy =  2.0*2.0 * std::sqrt(eps)/ (1 + cxy*cxy);
    g1xy    =  atan(cxy) + M_PI/2.0;
    g1xy_x  = -hxy*(x-0.5);
    g1xy_xx = -2*cxy*g1xy_x*g1xy_x - hxy;
    g1xy_y  = -hxy*(y-0.5);
    g1xy_yy = -2*cxy*g1xy_y*g1xy_y - hxy;
    g1xy    =  g1xy / M_PI; 
    g1xy_x  =  g1xy_x / M_PI;
    g1xy_xx =  g1xy_xx / M_PI;
    g1xy_y  =  g1xy_y / M_PI;
    g1xy_yy =  g1xy_yy / M_PI;
    
    g2xy    =   16.0*x*(x - 1)*y*(y - 1); 
    g2xy_x  =   16.0*(2*x - 1)*y*(y - 1);
    g2xy_xx =   16.0*2*y*(y - 1);
    g2xy_y  =   16.0*(2*y - 1)*x*(x - 1);
    g2xy_yy =   16.0*2*x*(x - 1);
    
    values[0] = g1xy*g2xy;  //original solution
    values[1] = (g1xy_x*g2xy + g1xy*g2xy_x); //x-derivative
    values[2] = (g1xy_y*g2xy + g1xy*g2xy_y); //y-derivative
    
    //intermediate calculations for finding the laplacian
    gxy_xx  =  g1xy_xx*g2xy + 2*g1xy_x*g2xy_x + g1xy*g2xy_xx;  
    gxy_yy  =  g1xy_yy*g2xy + 2*g1xy_y*g2xy_y + g1xy*g2xy_yy;
    
    values[3]  = (gxy_xx + gxy_yy); //laplacian
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
      x = 0.+Param;
      y = 0.;
      break;
    case 1:
      x = 1.;
      y = 0.+Param;
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
  double eps=DIFFUSION;
  double *coeff;
  double exact[4]; 
  
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
  
    coeff[0] = eps;
    coeff[1] = 2; //x-coordiante of convection b
    coeff[2] = 3; //y-coordiante of convection b
    coeff[3] = 2; //reaction coeffecient c
    
    Exact(X[i], Y[i], exact);
    
    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
  }
}
