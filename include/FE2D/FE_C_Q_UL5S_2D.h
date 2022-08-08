// ***********************************************************************
// UL5S element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL5S_2D_NDOF = 31;

// number of dofs on the closure of the joints
static int C_Q_UL5S_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_Q_UL5S_2D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_Q_UL5S_2D_J1[6] = {  5,  6,  7,  8,  9, 10 };
static int C_Q_UL5S_2D_J2[6] = { 10, 11, 12, 13, 14, 15 };
static int C_Q_UL5S_2D_J3[6] = { 15, 16, 17, 18, 19,  0 };

static int *C_Q_UL5S_2D_J[4] = { C_Q_UL5S_2D_J0, C_Q_UL5S_2D_J1,
                                 C_Q_UL5S_2D_J2, C_Q_UL5S_2D_J3 };

// number of inner dofs
static int C_Q_UL5S_2D_NInner = 11;

// array containing the numbers for the inner dofs
static int C_Q_UL5S_2D_Inner[11] = { 20, 21, 22, 23, 24, 25, 26, 27,
                                     28, 29, 30 };

// number of outer dofs
static int C_Q_UL5S_2D_NOuter = 20;

// array containing the numbers for the outer dofs
static int C_Q_UL5S_2D_Outer[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
