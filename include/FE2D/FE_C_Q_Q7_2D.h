//************************************************************
// Q7 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q7_2D_NDOF = 64;

// number of dofs on the closure of the joints
static int C_Q_Q7_2D_JointDOF = 8;

// which local dofs are on the joints
static int C_Q_Q7_2D_J0[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static int C_Q_Q7_2D_J1[8] = { 7, 15, 23, 31, 39, 47, 55, 63 };
static int C_Q_Q7_2D_J2[8] = { 63, 62, 61, 60, 59, 58, 57, 56 };
static int C_Q_Q7_2D_J3[8] = { 56, 48, 40, 32, 24, 16, 8, 0 };

static int *C_Q_Q7_2D_J[8] = {  C_Q_Q7_2D_J0,  C_Q_Q7_2D_J1,
                              C_Q_Q7_2D_J2,  C_Q_Q7_2D_J3 };

// number of inner dofs
static int C_Q_Q7_2D_NInner = 36;

// array containing the numbers for the inner dofs
static int C_Q_Q7_2D_Inner[36] = { 9, 10, 11, 12, 13, 14, 17, 18, 19, 20,
                                   21, 22, 25, 26, 27, 28, 29, 30, 33, 34,
                                   35, 36, 37, 38, 41, 42, 43, 44, 45, 46,
                                   49, 50, 51, 52, 53, 54 };

// number of outer dofs
static int C_Q_Q7_2D_NOuter = 28;

// array containing the numbers for the outer dofs
static int C_Q_Q7_2D_Outer[28] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 23,
                                   24, 31, 32, 39, 40, 47, 48, 55, 56,
                                   57, 58, 59, 60, 61, 62, 63 };
