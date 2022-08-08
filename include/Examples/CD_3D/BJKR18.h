//Example from [BJKR18]

double PECLET_NUMBER;
bool dirichlet_on_sides = false;
double eps=1e-7;

void ExampleFile()
{
  Output::root_info("Example", "BJKR18.h");
}

// exact solution is not known
auto& Exact = unknown_solution_3d;

// find out if the point (x,y,z) is on the rectangle
bool boundary_rectangle(double x, double y, double z)
{
  // this is fact returns true for points on the outer boundary
  
  //xy plane with z=0
  if(std::abs(z)<eps)
    return true;
  //xy plane with z=2
  else if(std::abs(z-2)<eps)
    return true;
  //yz plane with x=0
  else if(std::abs(x)<eps)
    return true;
  //yz plane with x=5
  else if(std::abs(x-5)<eps)
    return true;
  //xz plane with y=0
  else if(std::abs(y)<eps)
    return true;
  //xz plane with y=2
  else if(std::abs(y-2)<eps)
    return true;
  else
    return false;
}

// kind of boundary condition (needed for FE space)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  if(std::abs(x-5)<eps && std::abs(y-2)>eps && std::abs(z-2)>eps && std::abs(y)>eps && std::abs(z)>eps)
    cond=NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  bool outer_box = boundary_rectangle(x, y, z);
  if(!outer_box)
  {
    value = 0.; // Dirichlet
  }
  else
  {
    BoundCond bc;
    BoundCondition(-1, x,y,z, bc);
    if(bc == NEUMANN)
      value = 0;
    else
      value = 1;
  }
}

void BilinearCoeffs(int n_points, const double *x, const double *,
                    const double *, const double *const*, double **coeffs)
{
  double l_x;
  l_x=(0.19*x[0]*x[0]*x[0]-1.42*x[0]*x[0]+2.38*x[0])/4;
  for(int i = 0; i < n_points; ++i)
  {
    coeffs[i][0] = PECLET_NUMBER; //diffusion coefficient
    coeffs[i][1] = 1;   //ux
    coeffs[i][2] = l_x;   //uy
    coeffs[i][3] = l_x;   //uz
    coeffs[i][4] = 0;   //reaction coefficient
    coeffs[i][5] = 0;   //rhs
  }
}

