double diffCoeff;

void print_information_about_example()
{
  Output::print("Example: rotating_shapes - meant as 3D analogy of rotating_bodies");
}

bool isInBox(double x, double y, double z, double t)
{
  const double xCenterOfRotation = 0.5;
  const double yCenterOfRotation = 0.5;
//   const double zCenterOfRotation = 0.5;
  
  const double xCenterOfBoxInitial = 0.5;
  const double yCenterOfBoxInitial = 0.25;
  const double zCenterOfBoxInitial = 0.5;
  
  double xCenterOfBox = ( xCenterOfBoxInitial - xCenterOfRotation )*std::cos(t)
                      - ( yCenterOfBoxInitial - yCenterOfRotation )*std::sin(t)
                      + xCenterOfRotation;
                      
  double yCenterOfBox = ( xCenterOfBoxInitial - xCenterOfRotation )*std::sin(t)
                      + ( yCenterOfBoxInitial - yCenterOfRotation )*std::cos(t)
                      + yCenterOfRotation;
                      
  double zCenterOfBox = zCenterOfBoxInitial;
  
  const double halfBoxEdge = 0.125;
  
  if(   std::abs(xCenterOfBox - x) <= halfBoxEdge
     && std::abs(yCenterOfBox - y) <= halfBoxEdge
     && std::abs(zCenterOfBox - z) <= halfBoxEdge )
    return true;
  else
    return false;
}

bool isInCylinder(double x, double y, double z, double t)
{
  const double xCenterOfRotation = 0.5; 
  const double yCenterOfRotation = 0.5;
//   const double zCenterOfRotation = 0.5;
  
  const double xCenterOfCylInitial = 0.5;
  const double yCenterOfCylInitial = 0.75;
  const double zCenterOfCylInitial = 0.5;
  
  double xCenterOfCyl = ( xCenterOfCylInitial - xCenterOfRotation )*std::cos(t)
                      - ( yCenterOfCylInitial - yCenterOfRotation )*std::sin(t)
                      + xCenterOfRotation;
                      
  double yCenterOfCyl = ( xCenterOfCylInitial - xCenterOfRotation )*std::sin(t)
                      + ( yCenterOfCylInitial - yCenterOfRotation )*std::cos(t)
                      + yCenterOfRotation;
                      
  double zCenterOfCyl = zCenterOfCylInitial;
  
  const double halfCylHeight = 0.25;
  const double CylRadius = 0.125;
  
  if(   (xCenterOfCyl - x)*(xCenterOfCyl - x)
      + (yCenterOfCyl - y)*(yCenterOfCyl - y) <= CylRadius*CylRadius
     && std::abs(zCenterOfCyl - z) <= halfCylHeight )
    return true;
  else
    return false;
}

bool isInHollowCylinder(double x, double y, double z, double t)
{
  const double xCenterOfRotation = 0.5; 
  const double yCenterOfRotation = 0.5;
//   const double zCenterOfRotation = 0.5;
  
  const double xCenterOfCylInitial = 0.5;
  const double yCenterOfCylInitial = 0.75;
  const double zCenterOfCylInitial = 0.5;
  
  double xCenterOfCyl = ( xCenterOfCylInitial - xCenterOfRotation )*std::cos(t)
                      - ( yCenterOfCylInitial - yCenterOfRotation )*std::sin(t)
                      + xCenterOfRotation;
                      
  double yCenterOfCyl = ( xCenterOfCylInitial - xCenterOfRotation )*std::sin(t)
                      + ( yCenterOfCylInitial - yCenterOfRotation )*std::cos(t)
                      + yCenterOfRotation;
                      
  double zCenterOfCyl = zCenterOfCylInitial;
  
  const double halfCylHeight = 0.25;
  const double CylRadius = 0.125;
  const double InnerCylRadius = 0.0625;
  
  double distanceSquared = (xCenterOfCyl - x)*(xCenterOfCyl - x)
                         + (yCenterOfCyl - y)*(yCenterOfCyl - y);
  
  if(   distanceSquared <= CylRadius*CylRadius
     && distanceSquared >= InnerCylRadius*InnerCylRadius
     && std::abs(zCenterOfCyl - z) <= halfCylHeight )
    return true;
  else
    return false;
}

bool isInCone(double x, double y, double z, double t)
{
  const double xCenterOfRotation = 0.5;
  const double yCenterOfRotation = 0.5;
//   const double zCenterOfRotation = 0.5;
  
   const double xCenterOfConInitial = 0.75;
   const double yCenterOfConInitial = 0.5;
   const double zCenterOfConInitial = 0.5;
  
  double xCenterOfCon = ( xCenterOfConInitial - xCenterOfRotation )*std::cos(t)
                      - ( yCenterOfConInitial - yCenterOfRotation )*std::sin(t)
                      + xCenterOfRotation;
                      
  double yCenterOfCon = ( xCenterOfConInitial - xCenterOfRotation )*std::sin(t)
                      + ( yCenterOfConInitial - yCenterOfRotation )*std::cos(t)
                      + yCenterOfRotation;
                      
  double zCenterOfCon = zCenterOfConInitial;
  
  const double ConHeight = 0.5;
  const double halfConHeight = 0.5*ConHeight;
  const double ConRadius = 0.125;
  double zBottomFaceDist = z - zCenterOfConInitial - halfConHeight;
  
  if(   std::abs(zCenterOfCon - z) <= halfConHeight
     && (xCenterOfCon - x)*(xCenterOfCon - x)+ (yCenterOfCon - y)*(yCenterOfCon - y)
        <= (zBottomFaceDist*ConRadius/ConHeight)*(zBottomFaceDist*ConRadius/ConHeight) )
    return true;
  else
    return false;
}

bool isInBall(double x, double y, double z, double t)
{
  const double xCenterOfRotation = 0.5;
  const double yCenterOfRotation = 0.5;
//   const double zCenterOfRotation = 0.5;
  
  const double xCenterOfConInitial = 0.25;
  const double yCenterOfConInitial = 0.5;
  const double zCenterOfConInitial = 0.5;
  
  double xCenterOfCon = ( xCenterOfConInitial - xCenterOfRotation )*std::cos(t)
                      - ( yCenterOfConInitial - yCenterOfRotation )*std::sin(t)
                      + xCenterOfRotation;
                      
  double yCenterOfCon = ( xCenterOfConInitial - xCenterOfRotation )*std::sin(t)
                      + ( yCenterOfConInitial - yCenterOfRotation )*std::cos(t)
                      + yCenterOfRotation;
                      
  double zCenterOfCon = zCenterOfConInitial;
  
  const double radius = 0.125;
  const double radiusSquared = radius*radius;
  
  double distanceSquared = (xCenterOfCon - x)*(xCenterOfCon - x)
                         + (yCenterOfCon - y)*(yCenterOfCon - y)
                         + (zCenterOfCon - z)*(zCenterOfCon - z);
                         
  if( distanceSquared <= radiusSquared)
    return true;
  else
    return false;
}

void exactSolution(double x, double y, double z, double t, double *values)
{
   double val = ( isInBox(x, y, z, t) || isInHollowCylinder(x, y, z, t) || isInCone(x, y, z, t) ) ? 1.0 : 0.0;
//  double val = ( isInBox(x, y, z, t) || isInCylinder(x, y, z, t) || isInCone(x, y, z, t) ) ? 1.0 : 0.0;
//  double val = ( isInBox(x, y, z, t) || isInCone(x, y, z, t) ) ? 1.0 : 0.0;
   
  values[0] = val;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
}

void Exact(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  exactSolution(x, y, z, t, values);
}

void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void InitialCondition(double x, double y, double z, double *values)
{
  exactSolution(x, y, z, 0.0, values);
}

void BoundValue(int, double , double , double , double &value)
{
  value = 0.0;
}

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *, const double *const*, double **coeffs)
{
//   double diffCoeff = 1.0e-20;
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    double xVelocity = 0.5 - y[i];
    double yVelocity = x[i] - 0.5;
    double zVelocity = 0.0;
    
    coeff[0] = diffCoeff;
    coeff[1] = xVelocity;
    coeff[2] = yVelocity;
    coeff[3] = zVelocity;
    coeff[4] = 0.0;
    coeff[5] = 0.0; 
  }
}

