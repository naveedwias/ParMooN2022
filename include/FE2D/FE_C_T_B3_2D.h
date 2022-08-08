// ***********************************************************************
// conforming P3 element with cell bubbles
// ***********************************************************************

// number of degrees of freedom
static int C_T_B3_2D_NDOF = 12;

// number of dofs on the closure of the joints
static int C_T_B3_2D_JointDOF = 4;

// which local dofs are on the joints
static int C_T_B3_2D_J0[4] = { 0, 1, 2, 3 };
static int C_T_B3_2D_J1[4] = { 3, 5, 7, 8 };
static int C_T_B3_2D_J2[4] = { 8, 6, 4, 0 };

static int *C_T_B3_2D_J[3] = { C_T_B3_2D_J0, C_T_B3_2D_J1,
                              C_T_B3_2D_J2 };

// number of inner dofs
static int C_T_B3_2D_NInner = 3;

// array containing the numbers for the inner dofs
static int C_T_B3_2D_Inner[3] = { 9, 10, 11 };

// number of outer dofs
static int C_T_B3_2D_NOuter = 9;

// array containing the numbers for the outer dofs
static int C_T_B3_2D_Outer[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
