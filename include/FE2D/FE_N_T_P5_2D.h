// ***********************************************************************
// P5 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P5_2D_NDOF = 25;

// number of dofs on the closure of the joints
static int N_T_P5_2D_JointDOF = 5;

// which local dofs are on the joints
static int N_T_P5_2D_J0[5] = { 0, 3, 6,  9, 12 };
static int N_T_P5_2D_J1[5] = { 1, 4, 7, 10, 13 };
static int N_T_P5_2D_J2[5] = { 2, 5, 8, 11, 14 };
 
static int *N_T_P5_2D_J[3] = { N_T_P5_2D_J0, N_T_P5_2D_J1, N_T_P5_2D_J2 };

// number of inner dofs
static int N_T_P5_2D_NInner = 10;

// array containing the numbers for the inner dofs
static int N_T_P5_2D_Inner[10] = { 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };

// number of outer dofs
static int N_T_P5_2D_NOuter = 15;

// array containing the numbers for the outer dofs
static int N_T_P5_2D_Outer[15] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                   12, 13, 14 };
