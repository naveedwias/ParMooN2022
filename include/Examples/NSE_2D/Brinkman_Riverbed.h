// Brinkman 2D problem: Riverbed
// author: Laura Blank
/*
 * domain Omega = (0, 2) x (0, 2)
 * subdomains pure flow = (0, 2) x (1.5+-eps, 2)
 *            porous = (0, 2) x (0, 1.5+-eps)
 * Interface Gamma = {0, 2} x {1.5+-eps}
 * So the interface starts at (0, 1.5) following some function through
 * the domain and ending at (2, 1.5)

 Boundary conditions:
      Dirichlet on velocity on (0, 2) x {0, 2}
      periodic & Neumann on {0, 2} x (0, 2)
*/

// physical parameter
// These should be reset when constructing the Example class


double effective_viscosity = -1;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Riverbed.h");
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

void BoundCondition(int bdComp, double, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (bdComp == (int) neumann_id[j]) //(possibly 0 (left bottom), 2, 3, 5)
    {
      cond = NEUMANN;
      return;
    }
  }

  // cond = (bdComp==0 || bdComp==2 || bdComp==3 || bdComp==5) ? NEUMANN : DIRICHLET;

  // set Nitsche BC
  for (unsigned int j = 0; j < nitsche_id.size(); j++)
  {
    if (bdComp == (int)nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;
      return;
    }
  }
}

// ========================================================================
// boundary values
// ========================================================================

void U1BoundValue(int BdComp, double, double &value)
{
  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be alway zero!!!!
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (BdComp == (int)neumann_id[j])
    {
       switch(BdComp)
       {
         case 0:
           value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
           break;
         case 2:
           value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
           break;
         case 3:
           value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
           break;
         case 5:
           value = 0.; //0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
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
     case 1:
       value = 0; // 1.; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
       break;
     case 4:
       value = 0; // 1.; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
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
          case 0:
            value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
            break;
          case 2:
            value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
            break;
          case 3:
            value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
            break;
          case 5:
            value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
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
      case 1:
        value = 0; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
        break;
      case 4:
        value = 0; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
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
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][4] = sigma;
    coeffs[i][0] = effective_viscosity;

    if (sigma == -2 && effective_viscosity == -2)
    {
      // porous region
      if (  ( ( x[i] >= 0. ) && (x[i] <= 0.9) && (y[i] <= 1./9. * x[i] + 1.5) ) ||
            ( ( x[i] > 0.9 ) && (x[i] <= 1.) && (y[i] <= (-1.) * x[i] + 2.5) )  ||
            ( ( x[i] > 1. ) && (x[i] <= 1.9) && (y[i] <= 1./9. * x[i] + 25./18.) ) ||
            ( ( x[i] > 1.9 ) && ( x[i] <= 2.) && (y[i] <= (-1.) * x[i] + 3.5) ) )
      {
        coeffs[i][0] = 0.; //0.0001;
        coeffs[i][4] = 1.e7; //1./0.000001;
      }
      else  // flow region
      {
        coeffs[i][0] = 0.001; //0.000001;
        coeffs[i][4] = 0.;
        //coeffs[i][6] = 1000000.0;
      }
    }
    else
      ErrThrow("In the input file one has to set inverse_permeability: -2"
              " and effective_viscosity: -2 in order to use the riverbed "
              "example appropriately. "
              "(The actual values are set in Brinkman_Riverbed.h and in "
              "BoundaryAssembling2D.C for Nitsche in a fixed scenario.)");


    /* ***************************************** */
    // if sources and sinks via analytic_coefficient_function (delta distr.), then \neq 0
    coeffs[i][9] = 0;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0.;   // f1 (rhs of Brinkman problem for u1)
    coeffs[i][2] = 0.; // f2 (rhs of Brinkman problem for u2)

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = 0;
  }
}


