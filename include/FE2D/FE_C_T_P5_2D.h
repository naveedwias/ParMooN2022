// ***********************************************************************
// conforming P5 element
// ***********************************************************************

// number of degrees of freedom
static int C_T_P5_2D_NDOF = 21;

// number of dofs on the closure of the joints
static int C_T_P5_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_T_P5_2D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_T_P5_2D_J1[6] = {  5, 10, 14, 17, 19, 20 };
static int C_T_P5_2D_J2[6] = { 20, 18, 15, 11,  6,  0 };

static int *C_T_P5_2D_J[3] = { C_T_P5_2D_J0, C_T_P5_2D_J1,
                               C_T_P5_2D_J2 };

// number of inner dofs
static int C_T_P5_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int C_T_P5_2D_Inner[6] = { 7, 8, 9, 12, 13, 16 };

// number of outer dofs
static int C_T_P5_2D_NOuter = 15;

// array containing the numbers for the outer dofs
static int C_T_P5_2D_Outer[15] = { 0,  1,  2,  3,  4, 5, 6, 10, 11, 14,
                                  15, 17, 18, 19, 20 };
