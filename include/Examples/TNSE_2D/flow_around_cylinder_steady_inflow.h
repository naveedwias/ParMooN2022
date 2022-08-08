// Navier-Stokes problem, Benchmark problem
// 
// u(x,y) = unknown
// p(x,y) = unknown

#define __BENCH__

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  Output::print("Example: flow_around_cylinder_steady_inflow.h");

  ExampleOutput::DeclareOutputVariable("Drag");
  ExampleOutput::DeclareOutputVariable("Lift");
  ExampleOutput::DeclareOutputVariable("deltaP");
  ExampleOutput::DeclareOutputVariable("v_0.5_0.2_x");
  ExampleOutput::DeclareOutputVariable("v_0.5_0.2_y");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double y, double *values)
{
  values[0] = 6.0*y*(1-y);
}

void InitialU2(double, double, double *values)
{
  values[0] = 0;
}

void InitialP(double, double, double *values)
{
  values[0] = 0;
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
void BoundCondition(int bdcomp, double, BoundCond &cond)
{
  switch(bdcomp)
  {
    case 1:
      cond = NEUMANN;
      break;
    default:
      cond = DIRICHLET;
      break;
  }      
}

void BoundConditionPressure(int, double, BoundCond &cond)
{
     cond = NEUMANN;
}

void PressureBoundValue(int, double, double &value)
{
  value = 0;
}

void U1BoundValue(int BdComp, double y, double &value)
{

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;					
            break;
    case 2: value = 0;
            break;
    case 3: value = 6.0*y*(1-y);
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }

}

void U1BoundValue_diff(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value=1.0*M_PI/8*std::cos(M_PI*t/8)*6*Param*(1-Param);
            break;
    case 2: value = 0;
            break;
    case 3: value=1.0*M_PI/8*std::cos(M_PI*t/8)*6*Param*(1-Param);
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }

}

void U2BoundValue(int BdComp, double, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

void U2BoundValue_diff(int BdComp, double, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{
  double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // dot f1
    coeff[4] = 0; // dot f2
  }
}

/** calculate characteristic values */
void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,
             TFEFunction2D *pfct, 
	     TFEFunction2D *u1oldfct, TFEFunction2D *u2oldfct,
             double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Edges,comp;
  std::vector<const FiniteElement*> LocalUsedElements(2, nullptr);
  double **OrigFEValues, *Orig;
  bool SecondDer[2] = { false, false };
  double *u1, *u2, *p, *u1old, *u2old;
  int N_Cells;
  double value, value1, value2, value3, value4, value5;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValues2[MaxN_BaseFunctions2D];
  double FEFunctValues3[MaxN_BaseFunctions2D];
  double FEFunctValues4[MaxN_BaseFunctions2D];
  double FEFunctValues5[MaxN_BaseFunctions2D];
  int N_DerivativesU = 3;
  double *Derivatives[MaxN_BaseFunctions2D];
  MultiIndex2D NeededDerivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                        MultiIndex2D::D01 };
  TFEFunction2D *vfct;
  double *v, nu = DIMENSIONLESS_VISCOSITY;
  double *Der, *aux;
  TJoint *joint;
  TBoundEdge *boundedge;
  const TBoundComp *BoundComp;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";
  double dt = TDatabase::TimeDB->TIMESTEPLENGTH;
  if (dt < 1e-8)
  {
    ErrThrow("time step to small ");
  }
  
  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  u1old = u1oldfct->GetValues();
  u2old = u2oldfct->GetValues();
  p = pfct->GetValues();

  auto USpace = u1fct->GetFESpace2D();
  auto PSpace = pfct->GetFESpace2D();
  
  aux = new double [MaxN_QuadPoints_2D*12];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*12;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*sizeof(double));
  vfct = new TFEFunction2D(USpace, VString, v);

// ########################################################################
// loop over all cells
// ########################################################################
  auto Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==4) 
          {
            auto eleCell =  USpace->get_fe(i); // finite element of cell
            auto FEDesc = eleCell.GetFEDesc();   // fe descriptor
            N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
            DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
            auto DOF = USpace->GetGlobalDOF(i); // pointer to global dofs
            for (k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
              v[DOF[DOF_Circ[k]]] = 1;
          }
      }      
    }
  }
  
  cd = 0;
  cl = 0;

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    auto fe_u = USpace->get_fe(i);
    auto fe_p = PSpace->get_fe(i);
    LocalUsedElements[0] = &fe_u;
    LocalUsedElements[1] = &fe_p;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    // calculate all needed values of p 
    N_ = fe_p.GetN_DOF();

    auto DOF = PSpace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = FEDatabase::GetOrigElementValues(*fe_p.GetBaseFunct(),
                                                    MultiIndex2D::D00);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2 
    N_ = fe_u.GetN_DOF();

    DOF = USpace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = v[DOF[l]];
      FEFunctValues4[l] = u1old[DOF[l]];
      FEFunctValues5[l] = u2old[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*fe_u.GetBaseFunct(),
                                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
	} // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+4] = value2;
        Derivatives[j][k+7] = value3;
	if (k==0)
	{
	  value4 = 0;
	  value5 = 0;
	  for(l=0;l<N_;l++)
	    {
	      value4 += FEFunctValues4[l] * Orig[l];
	      value5 += FEFunctValues5[l] * Orig[l];
	    } // endfor l
	  Derivatives[j][k+10] = value4;
	  Derivatives[j][k+11] = value5;	  
	} // end if
      } // endfor j
    } // endfor k

    // calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      // Output::print(Der[1], " ", Der[4], " ", Der[10], " ", Der[11]);
      value1  = (Der[1]-Der[10])*Der[7]/dt+ nu*(Der[2]*Der[8]+Der[3]*Der[9]);
      value1 += (Der[1]*Der[2]+Der[4]*Der[3])*Der[7];
      value1 -= Der[0]*Der[8];

      value2  = (Der[4]-Der[11])*Der[7]/dt+ nu*(Der[5]*Der[8]+Der[6]*Der[9]);
      value2 += (Der[1]*Der[5]+Der[4]*Der[6])*Der[7];
      value2 -= Der[0]*Der[9];

      double weight = qf_orig.get_weight(j);
      cd += weight * value1;
      cl += weight * value2;
    }

  } // endfor i

  cd *= -20;
  cl *= -20;

  delete Derivatives[0];
  delete vfct;
  delete v;
}

void compute_drag_lift_pdiff(TimeNavierStokes<2>& time_nse2d, double&)
{
  double drag, lift, dP1[4], dP2[4], u1_at[4], u2_at[4];

  const TFEVectFunct2D& u(time_nse2d.get_velocity());
  TFEFunction2D& p(time_nse2d.get_pressure());

  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);

  /// @todo change implementation such that previous velocity can be used for computing drag and lift
  GetCdCl(u1.get(), u2.get(), &p, u1.get(), u2.get(), drag, lift);
  p.FindGradient(0.15, 0.2, dP1);
  p.FindGradient(0.25, 0.2, dP2);
  double pdiff = dP1[0]-dP2[0];

  u1->FindGradient(0.5, 0.2, u1_at);
  u2->FindGradient(0.5, 0.2, u2_at);

  double t = time_nse2d.get_time_stepping_scheme().current_time_;

  // print them reference values - f.y.i. some reference values are:
  Output::print(">>>>> Flow Around Cylinder (steady inflow) 2D: Postprocessing Output <<<<<");
  Output::print("time: ", t, " Drag = ",setprecision(16), drag);
  Output::print("time: ", t, " Lift = ", setprecision(16), lift);
  Output::print("time: ", t, " deltaP = ", setprecision(16), pdiff);
  Output::print("time: ", t, " velocity at (0.5 ; 0.2) = ", setprecision(16), u1_at[0], "  ", u2_at[0]);

  ExampleOutput::UpdateOutputVariable("Drag", drag);
  ExampleOutput::UpdateOutputVariable("Lift", lift);
  ExampleOutput::UpdateOutputVariable("deltaP", pdiff);
  ExampleOutput::UpdateOutputVariable("v_0.5_0.2_x", u1_at[0]);
  ExampleOutput::UpdateOutputVariable("v_0.5_0.2_y", u1_at[1]);
}

