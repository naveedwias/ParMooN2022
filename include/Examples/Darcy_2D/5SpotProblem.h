// Bechmark for Darcy problem, exact solution is known
// see Arif Masud, Thomas J.R. Hughes: "A stabilized mixed finite element
// method for Darcy flow", chapter 4.2
//
// Problem 1: The source term consists of two point sources at (0,0) and (1,1). 
// These are Dirac delta functions. To model these we use equivalent 
// prescribed flows at the first edges touching the corner points (0,0) and 
// (1,1). However for this to work we need to know the edge length and the  
// order of the used finite elements. That's why one has to modify the 
// function "FluxBoundValue" every time the grid is refined.
//
// Problem 2: On all other parts of the boundary we prescribe u.n=0 for 
// symmetry reasons. The analytical solution can be derived but it's rather
// complicated, see for example 
// Malik Mamode: "Fundamental solution of the Laplacian on flat tori and
//                boundary value problems for the planar Poisson equation in
//                rectangles"
// http://www.boundaryvalueproblems.com/content/2014/1/221
// Here we set the exact solution to zero for simplicity.

double mesh_width_at_corners = 1./32.;
unsigned int polynomial_order = 1;
double factor = 1.;

void ExampleFile()
{
  Output::print<1>("Example: 5SpotProblem.h with mesh width ",
                   mesh_width_at_corners, " and polynomial degree ",
                   polynomial_order, ". The hydraulic conductivity ", 
                   (factor == 1. ? "is uniform" : " has a checkerboard pattern"),
                   "(", (factor == 1. ? 2.0 : factor), ").");
  if(polynomial_order > 2)
  {
    Output::warn("5SpotProblem.h", "higher order polynomials are currently "
                 "not supported. Use with care!");
    polynomial_order = 2;
  } 
}

// ========================================================================
// exact solution
// ========================================================================
auto ExactU1 = unknown_solution_2d;
auto ExactU2 = unknown_solution_2d;
auto ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)  
{ cond = DIRICHLET; }
// u \cdot n
void FluxBoundValue(int bdComp, double Param, double &value)
{
  value = 0;
  double h = mesh_width_at_corners;
  int order = polynomial_order;
  switch(bdComp)
  {
    case 0: 
      if (Param<h)
      {
        value=-(1-Param/h)/(4*h);
        if( order == 2)
          value = -(3/(8*h))*(1-2*Param/h + Param*Param/(h*h));
      }
      break;
    case 1:
      if(Param>1-h)
      {
        value=(1-(1-Param)/h)/(4*h);
        if( order == 2)
          value = (3/(8*h))*(1-2*(1-Param)/h + (1-Param)*(1-Param)/(h*h));
      }
      break;
    case 2:
      if (Param<h)
      {
        value=(1-Param/h)/(4*h);
        if( order == 2)
          value = (3/(8*h))*(1-2*Param/h + Param*Param/(h*h));
      }
      break;
    case 3:
      if(Param>1-h)
      {
        value=-(1-(1-Param)/h)/(4*h);
        if( order == 2)
          value = -(3/(8*h))*(1-2*(1-Param)/h + (1-Param)*(1-Param)/(h*h));
      }
      break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
}
// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    if ( (X[i] < 0.5 && Y[i]<0.5) || (X[i] > 0.5 && Y[i]>0.5) )
      coeffs[i][0] = 2.0;
    else 
      coeffs[i][0] = 2.0 * factor;
    // RHS for exact solution
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // g
  }
}


