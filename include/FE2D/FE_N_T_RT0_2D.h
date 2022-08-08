// ***********************************************************************
// P0 Raviart-Thomas vector element, nonconforming , 2D
// History:  02.09.2010 implementation (Alfonso)
// ***********************************************************************

// number of degrees of freedom
static int N_T_RT0_2D_NDOF = 3;

// number of dofs on the closure of the joints
static int N_T_RT0_2D_JointDOF = 1;

// which local dofs are on the joints
static int N_T_RT0_2D_J0[1] = { 0 };
static int N_T_RT0_2D_J1[1] = { 1 };
static int N_T_RT0_2D_J2[1] = { 2 };
 
static int *N_T_RT0_2D_J[3] = { N_T_RT0_2D_J0, N_T_RT0_2D_J1,
                                 N_T_RT0_2D_J2};

// number of inner dofs
static int N_T_RT0_2D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_T_RT0_2D_Inner = nullptr;

// number of outer dofs
static int N_T_RT0_2D_NOuter = 3;

// array containing the numbers for the outer dofs
static int N_T_RT0_2D_Outer[3] = { 0, 1, 2};
