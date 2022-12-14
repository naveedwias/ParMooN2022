// ***********************************************************************
// P4 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P4_2D_NDOF = 18;

// number of dofs on the closure of the joints
static int N_T_P4_2D_JointDOF = 4;

// which local dofs are on the joints
static int N_T_P4_2D_J0[4] = { 0, 3, 6,  9 };
static int N_T_P4_2D_J1[4] = { 1, 4, 7, 10 };
static int N_T_P4_2D_J2[4] = { 2, 5, 8, 11 };
 
static int *N_T_P4_2D_J[3] = { N_T_P4_2D_J0, N_T_P4_2D_J1, N_T_P4_2D_J2 };

// number of inner dofs
static int N_T_P4_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int N_T_P4_2D_Inner[6] = { 12, 13, 14, 15, 16, 17 };

// number of outer dofs
static int N_T_P4_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_T_P4_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
