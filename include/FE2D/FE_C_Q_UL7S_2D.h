// ***********************************************************************
// UL7S element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL7S_2D_NDOF = 55;

// number of dofs on the closure of the joints
static int C_Q_UL7S_2D_JointDOF = 8;

// which local dofs are on the joints
static int C_Q_UL7S_2D_J0[8] = {  0,  1,  2,  3,  4,  5,  6,  7 };
static int C_Q_UL7S_2D_J1[8] = {  7,  8,  9, 10, 11, 12, 13, 14 };
static int C_Q_UL7S_2D_J2[8] = { 14, 15, 16, 17, 18, 19, 20, 21 };
static int C_Q_UL7S_2D_J3[8] = { 21, 22, 23, 24, 25, 26, 27,  0 };

static int *C_Q_UL7S_2D_J[4] = { C_Q_UL7S_2D_J0, C_Q_UL7S_2D_J1,
                                 C_Q_UL7S_2D_J2, C_Q_UL7S_2D_J3 };

// number of inner dofs
static int C_Q_UL7S_2D_NInner = 27;

// array containing the numbers for the inner dofs
static int C_Q_UL7S_2D_Inner[27] = { 28, 29, 30, 31, 32, 33, 34, 35,
                                     36, 37, 38, 39, 40, 41, 42, 43,
                                     44, 45, 46, 47, 48, 49, 50, 51,
                                     52, 53, 54 };

// number of outer dofs
static int C_Q_UL7S_2D_NOuter = 28;

// array containing the numbers for the outer dofs
static int C_Q_UL7S_2D_Outer[28] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                    20, 21, 22, 23, 24, 25, 26, 27 };
