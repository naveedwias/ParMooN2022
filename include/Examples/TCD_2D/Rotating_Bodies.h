/*!
 * Three rotating bodies - a hump, a cone and slotted cylinder, subject to a
 * counterclockwise flow.
 *
 * Use UnitSquare.GEO and UnitSquare.PRM for the geometry.
 *
 * Example file containing assembling information and boundary conditions
 * for a 2D time-dependent convection-diffusion example. This well-known
 * Benchmark problem originates from:
 *
 * R.J. Leveque (1996): High-Resolution Conservative Algorithms for Advection
 * in Incompressible Flow. SIAM J. Numer. Anal. 33(2), pp. 627-665
 *
 * @note: Some example specific functions for error measuring are
 * commented out so far, because they do not yet compile in ParMooN
 * and not needed at the moment.
 *
 * @author ??, import from MooNMD by Clemens Bartsch
 * @date 2015/11/22
 */

#ifndef __ROTATING_BODIES__
#define __ROTATING_BODIES__

//#include <LocalProjection.h>
//#include <MainRoutines2D.h>

/** Print out some information on the example file. */
void ExampleFile()
{
  Output::print("Example: Rotating_Bodies.h");
}

constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = false;
/** The exact solution */
void Exact(double x, double y, double *values)
{
  double hump, conical, x0, y0;
  double t = TDatabase::TimeDB->CURRENTTIME;

  //rotation
  x0 = std::cos(t)*x+std::sin(t)*y-0.5*std::cos(t)-0.5*std::sin(t)+0.5;
  y0 = -std::sin(t)*x+std::cos(t)*y+0.5*std::sin(t)-0.5*std::cos(t)+0.5;

  values[0] = 0;
  //cylinder
  if((x0-0.5)*(x0-0.5)+(y0-0.75)*(y0-0.75)<=0.15*0.15)
  {
    if(std::abs(x0-0.5) >= 0.0225 || y0 >= 0.85)
      values[0] = 1;
  }
  //conical body
  conical = 1.0/0.15*std::sqrt((x0-0.5)*(x0-0.5)+(y0-0.25)*(y0-0.25));
  if(conical <= 1.0)
    values[0] = 1.0 - conical;
  //hump
  hump = 1.0/0.15*std::sqrt((x0-0.25)*(x0-0.25)+(y0-0.5)*(y0-0.5));
  if(hump <= 1.0)
    values[0] = 0.25*(1.0+(std::cos(M_PI*hump)));

  values[1]=0;
  values[2]=0;
  values[3]=0;
}


/** The type of the boundary condition - Dirichlet on all components. */
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}


/** The value of boundary condition - 0 everywhere*/
void BoundValue(int, double, double &value)
{
  value = 0;
}


/** The initial condition - describes the initial shape of the
 *  bodies.*/
void InitialCondition(double x,  double y, double *values)
{
  double hump, conical;

  values[0] = 0;
  //cylinder
  if((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75)<= 0.15*0.15)
  {
    if(std::abs(x-0.5) >= 0.0225 || y >= 0.85)
      values[0] = 1.0;
  }
  //conical body
  conical = 1.0/0.15*std::sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25));
  if(conical <= 1.0)
    values[0] = 1.0 - conical;
  //hump
  hump = 1.0/0.15*std::sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5));
  if(hump <= 1.0)
    values[0] = 0.25*(1.0+(std::cos(M_PI*hump)));
}

/**
 * Coefficient function, used in the assembling process of the CDR problem.
 */
void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = 1.0e-20;

  int i;
  double *coeff;                                  // *param;
  double x, y;
  // double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0.5 - y;
    // convection in y direction
    coeff[2] = x - 0.5;
    // reaction
    coeff[3] = 0;
    // rhs
    coeff[4] = 0;
    // rhs from previous time step
    coeff[5] = 0;
  }
}

///****************************************************************/
///* computes errors in l1 and l2 norm                            */
///****************************************************************/
//void ComputeDiscreteErrors(TFESpace2D *fespace,
//TFEFunction2D *u,
//double *sol,
//double *errors,
//double *lump_mass)
//{
//  int i,j,index, N_U, N_Cells, *DOF;
//  double x, y, val[4], *err;
//  double l1_error, l2_error;
//  FE2D CurrentElement;
//  TCollection *coll;
//  TBaseCell *cell;
//
//  // get arrays with the numbering of the dof
//  N_U = fespace->get_n_dof();
//
//  // get collection and number of cells
//  coll = fespace->GetCollection();
//  N_Cells = coll->GetN_Cells();
//
//  // vector for error
//  err = new double[N_U];
//
//  // loop over cells
//  for(i=0;i<N_Cells;i++)
//  {
//    // next cell
//    cell = coll->GetCell(i);
//    // pointer to global indices of dof connected with this cell
//    DOF = fespace->GetGlobalDOF(i);
//    // finds finite element on cell i
//    CurrentElement = fespace->get_fe_type(i);
//    switch(CurrentElement)
//    {
//      // P_1
//      case C_P1_2D_T_A:
//	case C_UL1_2D_T_A:
//        for (j=0;j<3;j++)
//        {
//          // global index
//          index = DOF[j];
//          // coordinates of vertex
//          cell->GetVertex(j)->GetCoords(x, y);
//          Exact(x,y,val);
//          err[index]=val[0] - sol[index];
//        }
//        break;
//      case C_Q1_2D_Q_A:
//      case C_Q1_2D_Q_M:
//	case C_UL1_2D_Q_A:
//	case C_UL1_2D_Q_M:
//        for (j=0;j<4;j++)
//        {
//          // global index
//          if(j<2)
//            index = DOF[j];
//          if(j==2)
//            index = DOF[3];
//          if(j==3)
//            index = DOF[2];
//          // coordinates of vertex
//          cell->GetVertex(j)->GetCoords(x, y);
//          Exact(x,y,val);
//          err[index] = val[0] - sol[index];
//        }
//        break;
//      default:
//        ErrThrow("ComputeDiscreteErrors not implemented for element ",
//                 CurrentElement);
//        break;
//    }
//  }
//
//  l1_error = l2_error = 0.0;
//  for (i=0;i<N_U;i++)
//  {
//    l1_error += lump_mass[i] * std::abs(err[i]);
//    l2_error += lump_mass[i] * err[i]*err[i];
//
//  }
//  delete err;
//
//  errors[0] = l1_error;
//  errors[1] = std::sqrt(l2_error);
//}
//
//void ComputeExtremalValues(int N, double *sol, double  *values)
//{
//   int i;
//   double max, min;
//
//   min = 1e10;
//   max = -1e10;
//
//   for(i=0;i<N;i++)
//   {
//      if(sol[i] > max)
//         max = sol[i];
//      if(sol[i] < min)
//         min = sol[i];
//   }
//
//   values[0] = min;
//   values[1] = max;
//}
//
//void UltraLocalError_DG(TFEFunction2D *uh,
//			double *val)
//{
//  int i,j, N_Cells, N_Edges;
//  double *Values, min_uh = 1e10, max_uh = 1e-10, x, y, xs, ys;
//  double val_uh[4], val_u[4], eps=1e-8;
//  TCollection *coll;
//  TFESpace2D *fespace;
//  TBaseCell *cell;
//
//  // get coefficients of fe function
//  fespace = uh->GetFESpace2D();
//  Values = uh->GetValues();
//  // get collection
//  coll = fespace->GetCollection();
//  // get number of cells
//  N_Cells = coll->GetN_Cells();
//
//  // loop over the mesh cells
//  for(i=0;i<N_Cells;i++)
//  {
//    // get cell
//    cell = coll->GetCell(i);
//    // number of edges == number of vertices
//    N_Edges = cell->GetN_Edges();
//    xs = ys = 0;
//    // loop over the vertices
//    for (j=0;j< N_Edges;j++)
//    {
//	// get coordinates of the edges
//	cell->GetVertex(j)->GetCoords(x, y);
//	xs += x;
//	ys += y;
//    }
//    // bary center
//    xs /= N_Edges;
//    ys /= N_Edges;
//
//    // loop over the vertices
//    for (j=0;j< N_Edges;j++)
//    {
//	// get coordinates of the edges
//	cell->GetVertex(j)->GetCoords(x, y);
//	// go slighly away from the corner
//	x = x + eps * (xs-x);
//	y = y + eps * (ys-x);
//	// compute local function value
//	uh->FindGradientLocal(cell,i ,x, y, val_uh);
//	// check for maximum and minimum
//	if (val_uh[0] > max_uh)
//	    max_uh = val_uh[0];
//	if (val_uh[0] < min_uh)
//	    min_uh = val_uh[0];
//	// compute exact solutoin
//	//Exact(x,y,val_u);
//	// error in this node
//	//err[index]= val_u[0] - val_uh[0];
//    }
//  }
//  val[2] = min_uh;
//  val[3] = max_uh;
//} // UltraLocalProjection_DG
//
////**************************************************************
////  UltraLocalErrors_TCD
////  computes errors in the vertices of the mesh cells
////  for ultra local projection schemes for some TCD problems
////**************************************************************
//
//void UltraLocalError_TCD(TFEFunction2D *uh,
//        double *val, double *lumpmass)
//{
//  int i,j, N_Cells, N_Edges;
//  double *Values, min_uh = 1e10, max_uh = 1e-10, x, y, xs, ys;
//  double val_uh[4], val_u[4], eps=1e-8;
//  TCollection *coll;
//  TFESpace2D *fespace;
//  TBaseCell *cell;
//
//  // get coefficients of fe function
//  fespace = uh->GetFESpace2D();
//  Values = uh->GetValues();
//  // get collection
//  coll = fespace->GetCollection();
//  // get number of cells
//  N_Cells = coll->GetN_Cells();
//
//  // loop over the mesh cells
//  for(i=0;i<N_Cells;i++)
//  {
//    // get cell
//    cell = coll->GetCell(i);
//    // number of edges == number of vertices
//    N_Edges = cell->GetN_Edges();
//    // loop over the vertices
//    for (j=0;j< N_Edges;j++)
//    {
//	// get coordinates of the edges
//	cell->GetVertex(j)->GetCoords(x, y);
//	// compute local function value
//	uh->FindGradientLocal(cell,i ,x, y, val_uh);
//	// check for maximum and minimum
//	if (val_uh[0] > max_uh)
//	    max_uh = val_uh[0];
//	if (val_uh[0] < min_uh)
//	    min_uh = val_uh[0];
//	// compute exact solutoin
//	//Exact(x,y,val_u);
//	// error in this node
//	//err[index]= val_u[0] - val_uh[0];
//    }
//  }
//  val[2] = min_uh;
//  val[3] = max_uh;
//} // UltraLocalProjection_TCD
//
//void ComputeDataForEvaluationOfSolution(TFESpace2D *ConcentrationSpaces, TFEFunction2D *SolArray,
//					int N_Unknowns, double *sol, double *lump_mass, double *errors, bool is_localprojection)
//{
//    if (!is_localprojection &&(TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS!=2))
//    {
//      //ComputeDiscreteErrors(ConcentrationSpaces, SolArray,
//      //		      sol, errors, lump_mass);
//      //Output::print(TDatabase::TimeDB->CURRENTTIME, " l1: ", errors[0], " l2: ", errors[1]);
//      ComputeExtremalValues(N_Unknowns,sol,errors);
//      Output::print(TDatabase::TimeDB->CURRENTTIME, " min: ", errors[0],
//                    " max: ", errors[1], " ", errors[1] - errors[0]);
//	//exit(1);
//  }
//  else
//  {
//      if (is_localprojection)
//      {
//        UltraLocalError_TCD(SolArray, errors+5, lump_mass);
//        Output::print(TDatabase::TimeDB->CURRENTTIME, " min: ", errors[7],
//                      " max: ", errors[8], " ", errors[8] - errors[7]);
//      }
//      else
//      {
//	  // DG
//	  //ComputeExtremalValuesDG(N_Unknowns,sol,errors);
//	  UltraLocalError_DG(SolArray, errors+5);
//	  Output::print(TDatabase::TimeDB->CURRENTTIME, " min: ", errors[7],
//                  " max: ", errors[8], " ", errors[8] - errors[7]);
//      }
//  }
//  if(TDatabase::ParamDB->MEASURE_ERRORS)
//  {
//    //Output::print(TDatabase::TimeDB->CURRENTTIME, " mass of solution ",
//    //ComputeBilinearFormOfTwoFEFunctions2D_SameCells(ConcentrationSpaces->GetCollection(),
//    //                SolArray, SolArray, nullptr, 101));
//  }
//}
//

#endif
