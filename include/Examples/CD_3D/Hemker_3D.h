/**
 * This is the originally 2D Hemker example adopted to the 3D problem
 *
 * The domain is \f$\Omega = [-3,9] x [-3,3] x [0,6] \setminus C\f$ where 
 * \f$C = \{(x,y,z) | x^2+y^2 < 1\}\f$ is a cylinder. It is [0,6] in the 
 * z-direction because the sandwich meshes currently allow positive z, 
 * otherwise we would choose \f$[-3,3]\f$. Also we chose this interval so that
 * for each fixed x the cut \f$\{x\}\cap\Omega\f$ is a square.
 * 
 * We impose a constant flow field \f$b=(1,0,0)\f$ and no reaction. Also the 
 * right hand side is set to zero. The only non-zero contribution to the 
 * right hand side of the resulting system stems from the boundary condition 
 * at the cylinder.
 * 
 * The boundary conditions are Dirichlet at the inflow \f$x=-3\f$ and at the 
 * cylinder \f$x^2+y^2=1\f$. If the `bool dirichlet_on_sides` is true, Dirichlet
 * conditions are also set on the sides \f$y=-3\f$ and \f$y=3\f$. All other 
 * boundaries are of type Neumann.
 * 
 * All boundary values are zero (Dirichlet and Neumann) except on the cylinder
 * where they are set to 1.
 */
double PECLET_NUMBER;
double tolerance = 1e-5; // to check coordinates on boundary
bool dirichlet_on_sides = false;

void ExampleFile()
{
  Output::root_info("Example", "Hemker_3D.h");
}

// exact solution is not known
auto& Exact = unknown_solution_3d;

// find out if the point (x,y,z) is on the cylinder
bool on_cylinder(double x, double y, double)
{
  // this is fact returns true for points inside the cylinder. This is needed 
  // here because the fe space evaluates at the center of face to determine the 
  // boundary condition. This center is inside the cylinder.
  return (x*x + y*y - 1. < tolerance);
}

// kind of boundary condition (needed for FE space)
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  cond = NEUMANN;

  // Dirichlet at the inflow boundary and on cylinder
  bool inflow = (std::abs(3.0 + x) < tolerance);
  if(inflow || on_cylinder(x, y, z))
    cond = DIRICHLET;
  else if(dirichlet_on_sides && 
          (std::abs(y-3.0) < tolerance || std::abs(y+3.0) < tolerance))
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  value = on_cylinder(x, y, z) ? 1. : 0.;
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *, const double *const*, double **coeffs)
{
  for(int i = 0; i < n_points; ++i)
  {
    coeffs[i][0] = PECLET_NUMBER; //diffusion coefficient
    coeffs[i][1] = 1;   //ux
    coeffs[i][2] = 0;   //uy
    coeffs[i][3] = 0;   //uz
    coeffs[i][4] = 0;   //reaction coefficient
    coeffs[i][5] = 0;   //rhs
  }
}

