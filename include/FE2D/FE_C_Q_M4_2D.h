// ***********************************************************************
// UL4S element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_M4_2D_NDOF = 17;

// number of dofs on the closure of the joints
static int C_Q_M4_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_Q_M4_2D_J0[5] = {  0,  1,  2,  3,  4 };
static int C_Q_M4_2D_J1[5] = {  4,  5,  6,  7,  8 };
static int C_Q_M4_2D_J2[5] = {  8,  9, 10, 11, 12 };
static int C_Q_M4_2D_J3[5] = { 12, 13, 14, 15,  0 };

static int *C_Q_M4_2D_J[4] = { C_Q_M4_2D_J0, C_Q_M4_2D_J1,
                               C_Q_M4_2D_J2, C_Q_M4_2D_J3 };

// number of inner dofs
static int C_Q_M4_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_Q_M4_2D_Inner[1] = { 16 };

// number of outer dofs
static int C_Q_M4_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int C_Q_M4_2D_Outer[16] = { 0,  1,  2,  3,  4,  5, 6, 7, 8, 9,
                                  10, 11, 12, 13, 14, 15 };
