// ***********************************************************************
// P2 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P2_2D_NDOF = 7;

// number of dofs on the closure of the joints
static int N_T_P2_2D_JointDOF = 2;

// which local dofs are on the joints
static int N_T_P2_2D_J0[2] = { 0, 3 };
static int N_T_P2_2D_J1[2] = { 1, 4 };
static int N_T_P2_2D_J2[2] = { 2, 5 };
 
static int *N_T_P2_2D_J[3] = { N_T_P2_2D_J0, N_T_P2_2D_J1, N_T_P2_2D_J2 };

// number of inner dofs
static int N_T_P2_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int N_T_P2_2D_Inner[1] = { 6 };

// number of outer dofs
static int N_T_P2_2D_NOuter = 6;

// array containing the numbers for the outer dofs
static int N_T_P2_2D_Outer[6] = { 0, 1, 2, 3, 4, 5 };
