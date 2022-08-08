#include "math.h"

const double distanceTol = 1.0e-7;

/**
* @brief Rotating body problem as described in Ku20
* @author Abhinav Jha
*/

double DIFFUSION;

void ExampleFile()
{
  Output::root_info("Example", "steady_circular_advection.h");
}

bool isIn0Part(double x, double y)
{
  double radiusSq = x*x + y*y;
  
  if( radiusSq < 0.15*0.15 || 
    ( radiusSq > 0.45*0.45 && radiusSq < 0.55*0.55 ) || 
      radiusSq > 0.85*0.85)
    return true;
  else
    return false;
}

bool isIn1Part(double x, double y)
{
  double radiusSq = x*x + y*y;
  
  if(radiusSq >= 0.15*0.15 && radiusSq <= 0.45*0.45)
    return true;
  else
    return false;
}

bool isInCosPart(double x, double y)
{
  double radiusSq = x*x + y*y;
  
  if(radiusSq >= 0.55*0.55 && radiusSq <= 0.85*0.85)
    return true;
  else
    return false;
}

double valueIn0Part(double , double )
{
  return 0.0;
}

double valueIn1Part(double , double )
{
  return 1.0;
}

double valueInCosPart(double x, double y)
{
  double radius = std::sqrt(x*x + y*y);
  double cosPart = std::cos( 10.0 * M_PI * (radius - 0.7)/3.0 );
  
  return pow(cosPart, 2);
}

double partsValue(double x, double y)
{
  if(isIn0Part(x,y))
    return valueIn0Part(x,y);
  
  if(isIn1Part(x,y))
    return valueIn1Part(x,y);
  
  if(isInCosPart(x,y))
    return valueInCosPart(x,y);
    
  return 0.0; 
}


void Exact(double x, double y , double *values)
{
  
  
  if(isInCosPart(x,y))
  {
    double radius = std::sqrt(x*x + y*y);
    double sinCosArg = (20.0 * M_PI * (radius - 0.7))/3.0;
    double sinPart = std::sin( sinCosArg );
    double frontPart = 10.0 * M_PI / (3.0 * radius);
    double xyPart = frontPart * x*y;
    
    double cosPart = std::cos( sinCosArg );
    double secDerCommonFrontPart = xyPart * 20.0 * M_PI / (3.0 * radius);
    
    values[1] = - frontPart* xyPart * sinPart;
    values[2] = - frontPart* xyPart * sinPart;
    values[3] = secDerCommonFrontPart * cosPart * (y - x)
                + frontPart * sinPart * (y - x);
  }
  else
  {
    values[1] = 0.0;
    values[2] = 0.0;
    values[3] = 0.0;
  }
  values[0] = partsValue(x, y);
}

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
      x = 1. - Param;
      y = 1.;
      break;
    case 3:
      x = 0.;
      y = 1. - Param;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}


void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  // find coordinate point on the boundary
   double x, y;
   transform(BdComp, Param, x, y);
   if(std::abs(y - 0.0) < distanceTol)
     cond = NEUMANN;
    else
     cond = DIRICHLET;
}

void BoundValue(int BdComp, double Param, double &value)
{
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinate point on the boundary
  double x, y;
  transform(BdComp, Param, x, y);
  
  if( std::abs(x-0.0) <distanceTol )
    value = partsValue(x,y);
  else
    value = 0.0;
}


// kind of boundary condition (for FE space needed)
void BoundCondition_Poisson(int BdComp, double Param, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue_Poisson(int, double, double &value)
{
  value = 0.0;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double *coeff;
  double exact[4];
  
  for(int i = 0; i < n_points; ++i)
  {
    coeff = coeffs[i];
    //param = parameters[i];
    coeff[0] = 0.0;
    coeff[1] = Y[i];   //ux
    coeff[2] = -X[i];   //uy
    coeff[3] = 0.001;   //reaction coefficient

    Exact(X[i], Y[i], exact);
    
    coeff[4] = coeff[3] * exact[0];
  }
}

