// Brinkman 2D problem, 
/* Poiseuille-Problem from the Paper: 
   @Article{Hannukainen2011,
   author="Hannukainen, Antti
   and Juntunen, Mika
   and Stenberg, Rolf",
   title="Computations with finite element methods for the Brinkman problem",
   journal="Computational Geosciences",
   year="2011",
   month="Jan",
   day="01",
   volume="15",
   number="1",
   pages="155--166",
   abstract="Various finite element families for the Brinkman flow (or Stokes--Darcy flow) are tested numerically. Particularly, the effect of small permeability is studied. The tested finite elements are the MINI element, the Taylor--Hood element, and the stabilized equal order methods. The numerical tests include both a priori analysis and adaptive methods.",
   issn="1573-1499",
   doi="10.1007/s10596-010-9204-4",
   url="https://doi.org/10.1007/s10596-010-9204-4"
   }
 */

// authors: Alfonso Caiazzo and Laura Blank

/*
   Domain =  [0,1]x[0,1]
   t^2 = K * (mu_eff/mu)
   u(x,y) = [ K/mu * (1+std::exp(1/t)-std::exp((1-y)/t) - std::exp(y/t)) / (1+std::exp(1/t)  ,  0 ],   for t>0
   u(x,y) = [ K/mu , 0 ],    for t=0
   p(x,y) = -x + 1/2
 */


// physical parameter
// These should be reset when constructing the Example class


double effective_viscosity = -1;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Exponential_Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double, double y, double *values)
{
  double t = std::abs(std::sqrt(effective_viscosity/sigma));

  if (sigma == 0)
  {
    ErrThrow("The coefficient sigma is zero. Division by zero not allowed (for exact solution). ");
  }

  if (t == 0)
  {
    values[0] = 1/sigma;                                                               //u1
    values[1] = 0;                                                                  //u1_x
    values[2] = 0;                                                                  //u1_y
    values[3] = 0;                                                                  //Delta u1
  }
  else
  {
    values[0] = 1/sigma * (1+std::exp(1/t)-std::exp((1-y)/t) - std::exp(y/t)) / (1+std::exp(1/t));             //u1
    values[1] = 0;                                                                      //u1_x
    values[2] = 1/sigma * (std::exp((1-y)/t)-std::exp(y/t))/(t * (1+std::exp(1/t)));                      //u1_y
    values[3] = 1/sigma * (-std::exp((1-y)/t)-std::exp(y/t))/(t * t * (1+std::exp(1/t)));                 //Delta u1
  }
}

void ExactU2(double, double, double *values)
{
  values[0] = 0;            //u2
  values[1] = 0;            //u2_x
  values[2] = 0;            //u2_y
  values[3] = 0;            //Delta u2
}

void ExactP(double x, double, double *values)
{
  values[0] = 0.5-x;                      //p
  values[1] = -1;                         //p_x
  values[2] = 0;                          //p_y
  values[3] = 0;                          //Delta p=p_xx+p_yy
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
  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be alway zero!!!!
  for (unsigned int j = 0; j < neumann_id.size(); j++)
  {
    if (BdComp == (int)neumann_id[j])
    {
      switch(BdComp)
      {
        case 1:
          value = 0.;//TDatabase::ParamDB->neumann_boundary_value[j];
          break;
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

  double t = std::abs(std::sqrt(effective_viscosity/sigma));

  // loop to impose (strong or weak) Dirichlet
  switch(BdComp)
  {
    case 0: value = 0.;
            break;
    case 1: value = 1/sigma * (1+std::exp(1/t)-std::exp((1-Param)/t) - std::exp(Param/t)) / (1+std::exp(1/t));
            //value = 0; // TEST:sources/sinks,
            break;
    case 2: value = 0.;
            break;
    case 3: value = 1/sigma * (1+std::exp(1/t)-std::exp((Param)/t) - std::exp((1-Param)/t)) / (1+std::exp(1/t));
            //value = 0; // TEST: sources/sinks,
            break;
    default: cout << "No boundary component with this number." << endl;
             break;
  }
}

void U2BoundValue(int, double, double &value)
{
  value = 0.;
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
    /*
       if (x[i] < 0.5 && y[i] < 0.5)
       {
       coeffs[i][6] = 0.001;
       }
     */

    /* ***************************************** */
    // if sources and sinks via analytic_coefficient_function (delta distr.), then \neq 0
    coeffs[i][9]=0;
    /*
       double use_source_term = 1;
       if (use_source_term == 1)
       {
       coeffs[i][9]= parameters[i][0];
       }
     */

    // Adimensional parameter t^2 = mue/mu*K
    coeffs[i][0] = effective_viscosity;

    // (f1,f2)(x,y): RHS for momentum equation
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    ExactP(x[i], y[i], val_p);

    coeffs[i][1] = 0;
    //-coeff[5] * val_u1[3] - val_p[1] + (coeff[4]/coeff[6]) * val_u1[0];
    //(coeff[4]/coeff[6])-1;   // f1 (rhs of Brinkman problem for u1)

    coeffs[i][2] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2]; //0;
    coeffs[i][4] = sigma; //0;
  }
}


