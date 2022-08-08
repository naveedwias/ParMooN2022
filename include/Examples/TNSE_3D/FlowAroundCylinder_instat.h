/**
 * The flow around cylinder benchmark example as described e.g. in
 * John, V. (2002): Higher order finite element methods and multigrid
 * solvers in a benchmark problem for the 3D Navier--Stokes equations. Int. J.
 * Numer. Meth. Fluids 40: 775-798.
 * Boundary conditions are non-time dependent. Dirichlet at inflow, no-slip
 * Dirichlet at the walls and do-nothing conditions at the outflow.
 *
 * @note Make sure to use the correct domain and initial mesh description, as
 * the example won't check that.
 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

//side effect: sets the global parameter
void ExampleFile()
{
  Output::root_info<1>("Example","FlowAroundCylinder_instat.h");
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialU2(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialU3(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialP(double, double, double, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP = unknown_solution_3d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double x, double, double, BoundCond &cond)
{
  if (std::abs(x-2.5)<1e-6)
  {
    cond = NEUMANN;
  }
  else
    cond  = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  //if( (std::abs(x)<1e-6) || (std::abs(x-2.5)<1e-6) )
  if((std::abs(x)<1e-6)) //inflow boundary
  {
    double H = 0.41; //height of the cylinder in m
    double U = 2.25; //peak inflow velocity
    //value = 16*U*std::sin(M_PI*t/8)*y*z*(H-y)*(H-z) / std::pow(H,4);
    //value = 16*U*y*z*(H-y)*(H-z) / std::pow(H,4);
    value = 16*U*(0.5*std::sin(2*M_PI*t/0.5-M_PI/2.0)+1.5) *y*z*(H-y)*(H-z) / std::pow(H,4);
  }
  else
  {
    value = 0;
  }

}

// value of boundary condition
void U2BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{

  double eps = DIMENSIONLESS_VISCOSITY; // the kinematic viscosity (1e-3 in the paper cited above)

  for(int i=0;i<n_points;i++)
  {
    double* coeff = coeffs[i];
      
    coeff[0] = eps; // kinematic viscosity
    coeff[1] = 0;   // f1
    coeff[2] = 0;   // f2
    coeff[3] = 0;   // f3
    coeff[4] = 0;   // g
  }
}

//treat this like a private method
/** calculate characteristic values */
void get_cdrag_clift(TFEFunction3D *u1fct, TFEFunction3D *u2fct,
             TFEFunction3D *u3fct,
             TFEFunction3D *pfct,
             double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Faces,comp;
  double **OrigFEValues, *Orig;
  bool SecondDer[2] = { false, false };
  double *u1, *u2, *u3, *p;
  int N_Cells;
  double value, value1, value2, value3, value4;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValues3[MaxN_BaseFunctions3D];
  double FEFunctValues4[MaxN_BaseFunctions3D];
  int N_DerivativesU = 4;
  double *Derivatives[MaxN_QuadPoints_3D];
  MultiIndex3D NeededDerivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                        MultiIndex3D::D010, MultiIndex3D::D001 };
  TFEFunction3D *vfct;
  double *v;
  double nu = DIMENSIONLESS_VISCOSITY;
  double *Der, *aux;
  const TBoundComp3D *BoundComp;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  u3 = u3fct->GetValues();
  p = pfct->GetValues();

  auto USpace = u1fct->GetFESpace3D();
  auto PSpace = pfct->GetFESpace3D();

  aux = new double [MaxN_QuadPoints_3D*20];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*20;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*sizeof(double));
  vfct = new TFEFunction3D(USpace, VString, v);

// ########################################################################
// loop over all cells
// ########################################################################
  auto Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      auto joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge
      {

        auto boundface = (const TBoundFace *)joint;
        BoundComp = boundface->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id
        if ((comp>=4)&&(comp<=7))
          {
            auto FEEle = USpace->get_fe(i);   // finite element of cell
            auto FEDesc = FEEle.GetFEDesc();   // fe descriptor
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

// ########################################################################
// loop over all cells
// ########################################################################
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector<const FiniteElement*> LocalUsedElements(2, nullptr);
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

#ifdef _MPI
    if(cell->IsHaloCell())
    { //only perform the following calculations on OwnCells
      continue;
    }
#endif

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
    int N_Points = qf_orig.GetN_QuadPoints();

    // calculate all needed values of p
    auto BaseFunct = fe_p.GetBaseFunct();
    N_ = fe_p.GetN_DOF();

    auto DOF = PSpace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = FEDatabase::GetOrigElementValues(*BaseFunct,
                                                    MultiIndex3D::D000);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2, u3
    BaseFunct = fe_u.GetBaseFunct();
    N_ = fe_u.GetN_DOF();
    DOF = USpace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = u3[DOF[l]];
      FEFunctValues4[l] = v[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*BaseFunct,
                                                      NeededDerivatives[k]);

      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        value4 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
          value4 += FEFunctValues4[l] * Orig[l];
        } // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+5] = value2;
        Derivatives[j][k+9] = value3;
        Derivatives[j][k+13] = value4;
      } // endfor j
    } // endfor k

// calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      // Der[0] = p
      // Der[1] = u1
      // Der[2] = u1_x
      // Der[3] = u1_y
      // Der[4] = u1_z

      // nu * (u1_x*v_x+ u1_y*v_y + u1_z*v_z), v= (v,0,0)
      value1  = nu*(Der[2]*Der[14]+Der[3]*Der[15]+Der[4]*Der[16]);
      // (u1 * u1_x + u2* u1_y + u3* u1_z) * (1,0,0)
      value1 += (Der[1]*Der[2]+Der[5]*Der[3]+Der[9]*Der[4])*Der[13];
      // pressure times divergence of test function (1,0,0)
      value1 -= Der[0]*Der[14];

      value2  = nu*(Der[6]*Der[14]+Der[7]*Der[15]+Der[8]*Der[16]);
      value2 += (Der[1]*Der[6]+Der[5]*Der[7]+Der[9]*Der[8])*Der[13];
      value2 -= Der[0]*Der[15];

      double weight = qf_orig.get_weight(j);
      cd += weight * value1;
      cl += weight * value2;
    }

  } // endfor i

#ifdef _MPI
  //communicate the values of cd and cl and sum them up
  MPI_Comm comm = MPI_COMM_WORLD;

  double sendbuf[2] = {cd, cl};
  double recvbuf[2] = {0.0, 0.0};
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, comm);
  cd = recvbuf[0];
  cl = recvbuf[1];
#endif

  cd *= -20/0.41;
  cl *= -20/0.41;

  delete[] aux;
  delete[] v;
  delete vfct;
}

/**
 * Get the difference of pressure function p at point_A and point_B.
 * @param point_A Point A
 * @param point_B Point B
 * @param p The pressure function (or any other TFEFunction3D...)
 * @return THe difference of p's function values between points A and B.
 */
double get_p_diff(const std::array<double,3>& point_A,
                  const std::array<double,3>& point_B,
                  const TFEFunction3D& p)
{
  std::vector<double> dP1(4,0.0);
  std::vector<double> dP2(4,0.0);

  double pdiff = 0;

#ifdef _MPI
  {//this whole scope is only of interest in mpi case
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);

    //find out two processes which contain the points of interest
    int proc_A = p.GetFESpace3D()->GetCollection()->find_process_of_point(0.45, 0.2, 0.205);
    int proc_B = p.GetFESpace3D()->GetCollection()->find_process_of_point(0.55, 0.2, 0.205);

    //calculate value in point A
    if(my_rank==proc_A)
    {
      if(!p.FindGradient(point_A[0],point_A[1],point_A[2], dP1))
      {
        ErrThrow("Hey! I could not find point_A on the promised process!");
      }
    }

    //calculate value in point B
    if (my_rank==proc_B)
    {
      if(!p.FindGradient(point_B[0],point_B[1],point_B[2], dP2))
      {
        ErrThrow("Hey! I could not find point_B on the promised process!");
      }
    }

    // I'm using collective instead of point-to-point communication,
    // which is slower but right now much easier to implement...
    MPI_Bcast(&dP1.at(0), 4, MPI_DOUBLE, proc_A, comm);
    MPI_Bcast(&dP2.at(0), 4, MPI_DOUBLE, proc_B, comm);
    pdiff = dP1[0] - dP2[0];

  }
#else
  p.FindGradient(point_A[0],point_A[1],point_A[2], dP1);
  p.FindGradient(point_B[0],point_B[1],point_B[2], dP2);
  pdiff = dP1[0] - dP2[0];
#endif
  return pdiff;
}

// this is the actual interface
void compute_drag_lift_pdiff(TimeNavierStokes<3>& tnse3d)
{
#ifdef _MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);
#else
  int my_rank = 0;
#endif
  //compute drag and lift and print them
  double drag, lift;

  const TFEVectFunct3D& u(tnse3d.get_velocity());
  TFEFunction3D& p(tnse3d.get_pressure()); //I want this const!!!

  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
  auto u3 = u.GetComponent(2);

  get_cdrag_clift(u1.get(), u2.get(), u3.get(), &p, drag, lift);

  std::array<double,3> point_A = {{0.45, 0.2, 0.205}};
  std::array<double,3> point_B = {{0.55, 0.2, 0.205}};

  double pdiff = get_p_diff(point_A, point_B, p);

  // print them reference values - f.y.i. the reference intervals:
  // drag \in [6.05,6.25], lift \in [0.008,0.01], pdiff \in [0.165,0.175]
  // note: these hold for KINEMATIC_VISCOSITY = 1e-3 and geometry
  // as described in John 2002.
  if(my_rank == 0)
  {
    //Output::print(">>>>> Flow Around Cylinder (stat) 3D: Postprocessing Output <<<<<");
    //Output::print( " Drag = ",setprecision(16), drag);
    //Output::print( " Lift = ", setprecision(16), lift);
    //Output::print( " deltaP = ", setprecision(16), pdiff);
    Output::print(TDatabase::TimeDB->CURRENTTIME, " drag = ",setprecision(16), drag, " lift = ", setprecision(16), lift, " deltaP = ", setprecision(16), pdiff);
  }

}
