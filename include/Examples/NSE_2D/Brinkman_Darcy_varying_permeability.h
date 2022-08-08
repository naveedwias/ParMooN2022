/*********************************************************
 * Region with different permeabilities for the Brinkman 2D problem
 * - analogous to a Poiseuille flow
 * Stokes regime (sigma = 0):
   ux = Poiseuille solution (unitary pressure drop)
 * Darcy regime (effective_viscosity = 0):
   ux = 1/sigma
 * Brinkman regime:
   ux = exponential boundary layer
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
double sigma = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Darcy_varying_permeability.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double, double *values)
{
  // ratio between Stokes and Darcy terms
  double t = std::abs(std::sqrt( effective_viscosity/sigma));

  if ( std::abs(sigma) < 1e-10 )
  {
    // Stokes regime
    values[0] = 0.;// 1./(2.*effective_viscosity) * y * (1.-y);      
    values[1] = 0.;
    values[2] = 0.;// 1./(2.*effective_viscosity) * (1. - 2.*y);    
    values[3] = 0.;//-1./(effective_viscosity);
  }
  else if  (t == 0)
  {
    // Darcy regime
    values[0] = 0.;//1./sigma;  
    values[1] = 0.;
    values[2] = 0.;
    values[3] = 0.;
  }
  else
  {
    values[0] = 0.;//(1./sigma) * (1+std::exp(1/t)-std::exp((1-y)/t) - std::exp(y/t)) / (1+std::exp(1/t));      
    values[1] = 0.;
    values[2] = 0.;//(1./sigma) * (std::exp((1-y)/t)-std::exp(y/t))/(t * (1+std::exp(1/t)));    
    values[3] = 0.;//(1./sigma) * (-std::exp((1-y)/t)-std::exp(y/t))/(t * t * (1+std::exp(1/t))); 
  }
}

void ExactU2(double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */
void ExactP(double, double, double *values)
{
  values[0] = 0.;  //0.5-x;   
  values[1] = 0.;  //-1.;
  values[2] = 0.;
  values[3] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (i == (int)neumann_id[j])
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
  //double t = std::abs(std::sqrt(effective_viscosity/sigma));

  // loop to impose Neumann boundary conditions
  ///@attention we set =0 here, as Neumann BC are imposed using boundary assembling
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if ( BdComp == (int)neumann_id[j])
    {
      switch(BdComp)
      {
      case 1:
        value = 0.; 
        break;
      case 3:
        value = 0.; 
        break;
     case 0:
        value = 0.; 
        break;
      case 2:
        value = 0.; 
        break;
      default:
        Output::print("I cannot impose Neumann boundary condition on component ",
            BdComp);
        exit(1);
        break;
      }
      return;
    }
  }

  // loop to impose (strong or weak) Dirichlet
  switch(BdComp)
  {
  case 0: value = 0;
  break;
  case 2: value = 0;
  break;
  case 3:
        if (Param <= 0.4 || Param >= 0.6) 
         {
           value = 0.;
         } 
        else 
      {
value = -750 * 0.1326 * Param * Param + 750 * 0.1326 * Param - 180 * 0.1326; 
//value = -750 * 0.265 * Param * Param + 750 * 0.265 * Param - 180 * 0.265; 
//value = -50 * Param * Param + 50 * Param - 12;  //-50y^2+50y-12 
// value = -100 * Param * Param + 100 * Param -24;  //-100y^2+100y-24 
    }
    break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
}

void U2BoundValue(int, double, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, const double *x, const double *y,
               const double *const*, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  //double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    //ExactP(x[i], y[i], val_p);

    // physical parameters
    //coeffs[i][0] = effective_viscosity;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = 0.; //val_u1[1] + val_u2[2];
    //coeffs[i][4] = sigma;


// short channel
/*if (x[i] <= 0.75)
  {
  if (y[i] <= 0.4 || y[i] >= 0.6)
    {
    coeffs[i][4] = 10000.0; //100000000.0; //1000000.0; //  100000000000.0; //1000000.0; //1000000.0; // 1.; //
    coeffs[i][0] = 0.;
    //cout << "upper/lower layer: " <<  coeffs[i][4] << endl;
}
  else  //middle layer
    {
    coeffs[i][4] = 0.; //100000.0;
    coeffs[i][0] = 0.1;
    //cout << "middle layer: " <<  coeffs[i][4] << endl;
    }
  }
else 
  {
  coeffs[i][4] = 10000.0; //100000000.0; //1000000.0; // 100000000000.0;//1000000.0; //10000000000.0; //0.; //
  coeffs[i][0] = 0.; //0.01; //
  //cout << "right layer: " <<  coeffs[i][4] << endl;
  }
*/


// long channel
/*if (x[i] <= 20)
  {
  if (y[i] <= 0.4 || y[i] >= 0.6)
    {
    coeffs[i][4] = 10000.0; //100000000000.0; //1000000.0; //  100000000000.0; //1000000.0; //1000000.0; // 1.; //
    coeffs[i][0] = 0.;
    //cout << "upper/lower layer: " <<  coeffs[i][4] << endl;
}
  else  //middle layer
    {
    coeffs[i][4] = 0.; //100000.0;
    coeffs[i][0] = 0.1;
    //cout << "middle layer: " <<  coeffs[i][4] << endl;
    }
  }
else 
  {
  coeffs[i][4] = 10000.0; //100000000000.0; //1000000.0; // 100000000000.0;//1000000.0; //10000000000.0; //0.; //
  coeffs[i][0] = 0.; //0.01; //
  //cout << "right layer: " <<  coeffs[i][4] << endl;
  }
*/

// long long channel
if (x[i] <= 100)
  {
  if (y[i] <= 0.4 || y[i] >= 0.6)
    {
    coeffs[i][4] = 100000000000.0; //10000.0; //1000000.0; //  100000000000.0; //1000000.0; //1000000.0; // 1.; //
    coeffs[i][0] = 0.;
    //cout << "upper/lower layer: " <<  coeffs[i][4] << endl;
}
  else  //middle layer
    {
    coeffs[i][4] = 0.; //100000.0;
    coeffs[i][0] = 0.1;
    //cout << "middle layer: " <<  coeffs[i][4] << endl;
    }
  }
else 
  {
  coeffs[i][4] = 100000000000.0; //10000.0; //1000000.0; // 100000000000.0;//1000000.0; //10000000000.0; //0.; //
  coeffs[i][0] = 0.; //0.01; //
  //cout << "right layer: " <<  coeffs[i][4] << endl;
  }

 }
}


