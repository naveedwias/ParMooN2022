// Bechmark for Darcy problem, exact solution of the flux is linear 
// 

// choose different obstacles:
// 0 - circle in center
// 1 - square in center
// 2 - diamond in center
// 3 - two boxes, one at the top, one at the bottom
// 4 - two half circles, one at the top, one at the bottom
// 5 - two boxes, one at the top, one at the bottom, boxes have an offset
// 6 - two half ellipsoids with an offset, one at the top, one at the bottom,
// 7 - continuous change in permeability in form of a hump
int obstacle = 0;
// inside the obstacle the porousity sigma is multiplied by this factor
double factor = 10.;


void ExampleFile()
{
  obstacle = TDatabase::ParamDB->P0;
  Output::print<1>("Example: Obstacle.h.");
  Output::print<1>("The obstacle is ", factor, " times less conductive and ",
                   "has the shape of");
  if(factor == 1.)
  {
    Output::warn("Obstacle.h",
                 "The permeability factor is 1.0, which effectively means "
                 "there is no obstacle. Set 'darcy_permeability_jump' to a "
                 "value larger than one.");
  }
  switch(obstacle)
  {
    case 0:
      Output::print<1>(" a circle in the center");
      break;
    case 1:
      Output::print<1>(" a square in the center");
      break;
    case 2:
      Output::print<1>(" a diamond in the center");
      break;
    case 3:
      Output::print<1>(" two boxes on one vertical line, one at the top, ",
                       "the other at the bottom\n");
      break;
    case 4:
      Output::print<1>(" two half circles on one vertical line, one at the ",
                       "top, one at the bottom");
      break;
    case 5:
      Output::print<1>(" two boxes with a vertical offset, on at the top, one ",
                       "at the bottom");
      break;
    case 6:
      Output::print<1>("two half ellipsoids with a vertical offset, on at the ",
                       "top, one at the bottom");
      break;
    case 7:
      Output::print<1>("hump in the center (continuous change in permeability)");
      break;
    default:
      ErrThrow("Unknown obstacle. Set P0 in the input file to 0,...,7");
      break;
  }
}

// ========================================================================
// multi indices used for various things
// ========================================================================
// no exact solution known
auto ExactU1 = unknown_solution_2d;
auto ExactU2 = unknown_solution_2d;
auto ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)  
{
  cond = DIRICHLET;
}

// u \cdot n
void FluxBoundValue(int bdComp, double, double &value)
{
  switch(bdComp)
  {
    case 0: 
      value = 0.0;
      break;
    case 1: 
      value = 1.0;
      break;
    case 2: 
      value = 0.0;
      break;
    case 3:
      value = -1.0;
      break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
}

// ========================================================================
// coefficient and right hand side
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  double r;
  for(int i = 0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    
    coeffs[i][0] = eps;
    switch (obstacle)
    {
      case 0:
        r = std::sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
        if( r < 0.2 )
          coeffs[i][0] *= factor;
        break;
      case 1:
        if( x>0.45 && x<0.55 && y>0.35 && y<0.65)
          coeffs[i][0] *= factor;
        break;
      case 2:
        if (x+y>0.8 && x+y<1.2 && x-y<0.2 && x-y>-0.2)
          coeffs[i][0] *= factor;
        break;
      case 3:
        if(x>0.35 && x<0.65 && (y<0.2 || y>0.8))
          coeffs[i][0] *= factor;
        break;
      case 4:
        r = std::sqrt((x-0.5)*(x-0.5)+y*y);
        if( r < 0.2 )
          coeffs[i][0] *= factor;
        r = std::sqrt((x-0.5)*(x-0.5)+(y-1)*(y-1));
        if( r < 0.2 )
          coeffs[i][0] *= factor;
        break;
      case 5:
        if((x>0.65 && x<0.75 && y<0.6) || (x>0.25 && x<0.35 && y>0.4))
          coeffs[i][0] *= factor;
        break;
      case 6:
        r = std::sqrt((x-0.75)*(x-0.75)+0.1*y*y);
        if( r < 0.2 )
          coeffs[i][0] *= factor;
        r = std::sqrt((x-0.25)*(x-0.25)+0.1*(y-1)*(y-1));
        if( r < 0.2 )
          coeffs[i][0] *= factor;
        break;
      case 7:
        r = std::sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
          if( r < 0.2 )
            coeffs[i][0] *= 1+0.5*factor*(std::cos(5*M_PI*r)+1);
        break;
    }
    
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // g
  }
}
