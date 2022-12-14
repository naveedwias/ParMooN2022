// ***********************************************************************
// UL9SE element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL9SE_2D_NDOF = 85;

// number of dofs on the closure of the joints
static int C_Q_UL9SE_2D_JointDOF = 10;

// which local dofs are on the joints
static int C_Q_UL9SE_2D_J0[10] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9 };
static int C_Q_UL9SE_2D_J1[10] = {  9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
static int C_Q_UL9SE_2D_J2[10] = { 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
static int C_Q_UL9SE_2D_J3[10] = { 27, 28, 29, 30, 31, 32, 33, 34, 35,  0 };

static int *C_Q_UL9SE_2D_J[4] = { C_Q_UL9SE_2D_J0, C_Q_UL9SE_2D_J1,
                                  C_Q_UL9SE_2D_J2, C_Q_UL9SE_2D_J3 };

// number of inner dofs
static int C_Q_UL9SE_2D_NInner = 49;

// array containing the numbers for the inner dofs
static int C_Q_UL9SE_2D_Inner[49] = { 36, 37, 38, 39, 40, 41, 42, 43,
                                      44, 45, 46, 47, 48, 49, 50, 51,
                                      52, 53, 54, 55, 56, 57, 58, 59,
                                      60, 61, 62, 63, 64, 65, 66, 67,
                                      68, 69, 70, 71, 72, 73, 74, 75,
                                      76, 77, 78, 79, 80, 81, 82, 83,
                                      84 };

// number of outer dofs
static int C_Q_UL9SE_2D_NOuter = 36;

// array containing the numbers for the outer dofs
static int C_Q_UL9SE_2D_Outer[36] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                     10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                     30, 31, 32, 33, 34, 35 };
