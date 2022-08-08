// Brinkman 2D problem, 
/* from the Paper:
   @Article{Codina:
   A stabilized finite element method for generalized stationary 
   incompressible flows; 2001
   }
 */

// author: Laura Blank

/*
   - Domain =  [-4,4] x [0,6] \cup [-7,7] x [6,9]
   - unknown solution
   - bc: parabolic velocity profile
         u = (u1(y), 0)^T,
         u1(y):= 1-(y/3 - 3)^2 at inlet x=-7
   - bc: free flow at x=7
   - bc: u = (0,1)^T at y=9
   - bc: u = (0,0) else
 */


// physical parameter
// These should be reset when constructing the Example class


double effective_viscosity = -1;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Porous_Cavity_Tshape.h");
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_2d;
auto& ExactU2 = unknown_solution_2d;
auto& ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================

void BoundCondition(int i, double, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (i == (int) neumann_id[j])
    {
      cond = NEUMANN;
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

void U1BoundValue(int BdComp, double Param, double &value)
{
  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be alway zero!!!!
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (BdComp == (int)neumann_id[j])
    {
      switch(BdComp)
      {
        case 3:
          value = 0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }


// loop to impose (strong or weak) Dirichlet
switch(BdComp)
{
  case 0: value = 0.;
           break;
  case 1: value = 0.;
          break;
  case 2: value = 0.;
          break;
  //case 3: value = neumann;
  //        break;
  case 4: value = 1.;
           break;
  case 5: value = (-(1-Param)*(1-Param) + 2*(1-Param));
          break;
  case 6: value = 0.;
          break;
  case 7: value = 0.;
          break;
  default: cout << "No boundary component with this number." << endl;
           break;
}
}

void U2BoundValue(int BdComp, double, double &value)
{
  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be alway zero!!!!
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (BdComp == (int)neumann_id[j])
    {
      switch(BdComp)
      {
        case 3:
          value = 0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }

switch(BdComp)
{
  case 0: value = 0.;
           break;
  case 1: value = 0.;
          break;
  case 2: value = 0.;
          break;
  //case 3: value = neumann;
  //        break;
  case 4: value = 0.;
           break;
  case 5: value = 0.;
          break;
  case 6: value = 0.;
          break;
  case 7: value = 0.;
          break;
  default: cout << "No boundary component with this number." << endl;
           break;
}
}
// ========================================================================
// coefficients for Brinkman problem: viscosity, effective viscosity, permeability, f1, f2, g
// (lhs and rhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================
void LinCoeffs(int n_points, const double *x, const double *y,
               const double *const*, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = effective_viscosity;

    // (f1,f2)(x,y): RHS for momentum equation
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    ExactP(x[i], y[i], val_p);

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2];
    coeffs[i][4] = sigma;

    coeffs[i][1] = -coeffs[i][0] * val_u2[3]  + val_p[2] + coeffs[i][4] * val_u2[0];  // f1
    coeffs[i][2] = -coeffs[i][0] * val_u2[3]  + val_p[2] + coeffs[i][4] * val_u2[0];  // f2
  }
}


