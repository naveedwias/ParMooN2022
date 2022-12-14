// ***********************************************************************
// UL8SE element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL8SE_2D_NDOF = 68;

// number of dofs on the closure of the joints
static int C_Q_UL8SE_2D_JointDOF = 9;

// which local dofs are on the joints
static int C_Q_UL8SE_2D_J0[9] = {  0,  1,  2,  3,  4,  5,  6,  7,  8 };
static int C_Q_UL8SE_2D_J1[9] = {  8,  9, 10, 11, 12, 13, 14, 15, 16 };
static int C_Q_UL8SE_2D_J2[9] = { 16, 17, 18, 19, 20, 21, 22, 23, 24 };
static int C_Q_UL8SE_2D_J3[9] = { 24, 25, 26, 27, 28, 29, 30, 31,  0 };

static int *C_Q_UL8SE_2D_J[4] = { C_Q_UL8SE_2D_J0, C_Q_UL8SE_2D_J1,
                                  C_Q_UL8SE_2D_J2, C_Q_UL8SE_2D_J3 };

// number of inner dofs
static int C_Q_UL8SE_2D_NInner = 36;

// array containing the numbers for the inner dofs
static int C_Q_UL8SE_2D_Inner[36] = { 32, 33, 34, 35, 36, 37, 38, 39,
                                      40, 41, 42, 43, 44, 45, 46, 47,
                                      48, 49, 50, 51, 52, 53, 54, 55,
                                      56, 57, 58, 59, 60, 61, 62, 63,
                                      64, 65, 66, 67 };

// number of outer dofs
static int C_Q_UL8SE_2D_NOuter = 32;

// array containing the numbers for the outer dofs
static int C_Q_UL8SE_2D_Outer[32] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                     10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                     30, 31 };
