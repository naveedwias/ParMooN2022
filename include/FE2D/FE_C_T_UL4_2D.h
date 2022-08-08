// ***********************************************************************
// P4 element with bubbles, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_T_UL4_2D_NDOF = 22;

// number of dofs on the closure of the joints
static int C_T_UL4_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_T_UL4_2D_J0[5] = {  0,  1,  2,  3,  4 };
static int C_T_UL4_2D_J1[5] = {  4,  5,  6,  7,  8 };
static int C_T_UL4_2D_J2[5] = {  8,  9, 10, 11,  0 };

static int *C_T_UL4_2D_J[3] = { C_T_UL4_2D_J0, C_T_UL4_2D_J1, 
                                C_T_UL4_2D_J2 };

// number of inner dofs
static int C_T_UL4_2D_NInner = 10;

// array containing the numbers for the inner dofs
static int C_T_UL4_2D_Inner[10] = { 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };

// number of outer dofs
static int C_T_UL4_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int C_T_UL4_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
