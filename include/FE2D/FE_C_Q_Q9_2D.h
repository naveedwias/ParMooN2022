//************************************************************
// Q9 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q9_2D_NDOF = 100;

// number of dofs on the closure of the joints
static int C_Q_Q9_2D_JointDOF = 10;

// which local dofs are on the joints
static int C_Q_Q9_2D_J0[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static int C_Q_Q9_2D_J1[10] = { 9, 19, 29, 39, 49, 59, 69, 79, 89, 99 };
static int C_Q_Q9_2D_J2[10] = { 99, 98, 97, 96, 95, 94, 93, 92, 92, 90 };
static int C_Q_Q9_2D_J3[10] = { 90, 80, 70, 60, 50, 40, 30, 20, 10, 0 };

static int *C_Q_Q9_2D_J[10] = {  C_Q_Q9_2D_J0,  C_Q_Q9_2D_J1,
                               C_Q_Q9_2D_J2,  C_Q_Q9_2D_J3 };

// number of inner dofs
static int C_Q_Q9_2D_NInner = 64;

// array containing the numbers for the inner dofs
static int C_Q_Q9_2D_Inner[64] = { 11, 12, 13, 14, 15, 16, 17, 18, 21, 22,
                                   23, 24, 25, 26, 27, 28, 31, 32, 33, 34,
                                   35, 36, 37, 38, 41, 42, 43, 44, 45, 46,
                                   47, 48, 51, 52, 53, 54, 55, 56, 57, 58,
                                   61, 62, 63, 64, 65, 66, 67, 68, 71, 72,
                                   73, 74, 75, 76 ,77, 78, 81, 82, 83, 84,
                                   85, 86, 87, 88 };

// number of outer dofs
static int C_Q_Q9_2D_NOuter = 36;

// array containing the numbers for the outer dofs
static int C_Q_Q9_2D_Outer[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                   19, 20, 29, 30, 39, 40, 49, 50, 59, 60,
                                   69, 70, 79, 80, 89, 90, 91, 92, 93, 94,
                                   95, 96, 97, 98, 99 };
