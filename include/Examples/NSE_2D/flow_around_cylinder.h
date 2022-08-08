/**
 * The stationary flow around cylinder problem in 2D. A benchmark problem for
 * flow solvers, it is described in detail e.g. in:
 *
 * V. John, G. Matthies, "Higher Order Finite Element Discretizations in a
 * Benchmark Problem for Incompressible Flows", Int. J. Num. Meth. Fluids 37,
 * 885 - 903, 2001
 *
 */

#include<BoundEdge.h>

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;
constexpr double reference_nu = 0.001;
constexpr double reference_drag = 5.57953523384;
constexpr double reference_lift = 0.010618948146;
constexpr double reference_deltap = 0.11752016697;

//side effect: sets the global parameter
void ExampleFile()
{
  Output::print("Example: flow_around_cylinder.h (stationary, 2D)");
  if(DIMENSIONLESS_VISCOSITY == reference_nu)
  {
    Output::print("         This is Example D.5 in Volker John's book "
                  "\"Finite Element Methods for Incompressible Flow Problems\""
                  ", 2016");
    Output::print(setprecision(12),"         Published values are: drag = ",
                  reference_drag);
    Output::print(setprecision(12), "                               lift = ",
                  reference_lift);
    Output::print(setprecision(12), "                               deltaP = ",
                  reference_deltap);
  }
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
  if (i==1)
  {
    cond = NEUMANN;
  }
  else
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value= 0;
            break;
    case 2: value = 0;
            break;
    case 3: value=1.2*Param*(1-Param); // 4*0.3
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
    coeff[3] = 0; // g (divergence)
     
    // additional coefficient (used only in the Brinkman problem)
    coeff[4] = 0.;
  }
}

/** calculate characteristic values */
void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,
             TFEFunction2D *pfct, double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Edges,comp;
  std::vector<const FiniteElement*> LocalUsedElements(2, nullptr);
  double **OrigFEValues, *Orig;
  bool SecondDer[2] = { false, false };
  double *u1, *u2, *p;
  int N_Cells;
  double value, value1, value2, value3;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValues2[MaxN_BaseFunctions2D];
  double FEFunctValues3[MaxN_BaseFunctions2D];
  int N_DerivativesU = 3;
  double *Derivatives[MaxN_BaseFunctions2D];
  MultiIndex2D NeededDerivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                        MultiIndex2D::D01 };
  double *v, nu = DIMENSIONLESS_VISCOSITY;
  double *Der, *aux;
  const TBoundComp *BoundComp;
  int N_DOF_Circ, *DOF_Circ;

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  p = pfct->GetValues();

  auto USpace = u1fct->GetFESpace2D();
  auto PSpace = pfct->GetFESpace2D();


  aux = new double [MaxN_QuadPoints_2D*10];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*10;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*sizeof(double));

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
      auto joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge
      {

        auto boundedge=(const TBoundEdge *)joint;
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id
        if (comp==4)
          {
            auto eleCell = USpace->get_fe(i);
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
    auto& fe_u = USpace->get_fe(i);
    auto& fe_p = PSpace->get_fe(i);
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
      } // endfor j
    } // endfor k

    // calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];

      // nu * (u1_x*v_x, u1_y*v_y), v= (v,0)
      value1  = nu*(Der[2]*Der[8]+Der[3]*Der[9]);
      // (u1 * u1_x + u2* u1_y) * (1,0)
      value1 += (Der[1]*Der[2]+Der[4]*Der[3])*Der[7];
      // pressure times divergence of test function (1,0)
      value1 -= Der[0]*Der[8];

      value2  = nu*(Der[5]*Der[8]+Der[6]*Der[9]);
      value2 += (Der[1]*Der[5]+Der[4]*Der[6])*Der[7];
      value2 -= Der[0]*Der[9];

      double weight = qf_orig.get_weight(j);
      cd += weight * value1;
      cl += weight * value2;
    }

  } // endfor i

  cd *= -500;
  cl *= -500;

  delete Derivatives[0];
  delete v;
}

void compute_drag_lift_pdiff(NavierStokes<2>& nse2d)
{
  double drag, lift, dP1[4], dP2[4];

  const TFEVectFunct2D& u(nse2d.get_velocity());
  TFEFunction2D& p(nse2d.get_pressure());

  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);

  GetCdCl(u1.get(), u2.get(), &p, drag, lift);
  p.FindGradient(0.15, 0.2, dP1);
  p.FindGradient(0.25, 0.2, dP2);

  double pdiff = dP1[0]-dP2[0];

  Output::print(">>>>> Flow Around Cylinder (stat) 2D: Postprocessing Output <<<<<");
  if(DIMENSIONLESS_VISCOSITY == reference_nu)
  {
    Output::print(" Drag = ", setprecision(16), drag, "  \terror: ",
                  std::abs(drag - reference_drag), "  \trelative_error: ",
                  std::abs(drag - reference_drag)/drag);
    Output::print(" Lift = ", setprecision(16), lift, "  \terror: ",
                  std::abs(lift - reference_lift), "  \trelative_error: ",
                  std::abs(lift - reference_lift)/lift);
    Output::print(" deltaP = ", setprecision(16), pdiff, "  \terror: ",
                  std::abs(pdiff - reference_deltap), "  \trelative_error: ",
                  std::abs(pdiff - reference_deltap)/pdiff);
  }
  else
  {
    Output::print(" Drag = ", setprecision(16), drag);
    Output::print(" Lift = ", setprecision(16), lift);
    Output::print(" deltaP = ", setprecision(16), pdiff);
  }
}
