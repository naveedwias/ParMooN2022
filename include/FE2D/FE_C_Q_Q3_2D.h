// ***********************************************************************
// Q3 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q3_2D_NDOF = 16;

// number of dofs on the closure of the joints
static int C_Q_Q3_2D_JointDOF = 4;

// which local dofs are on the joints
static int C_Q_Q3_2D_J0[4] = {  0,  1,  2,  3 };
static int C_Q_Q3_2D_J1[4] = {  3,  7, 11, 15 };
static int C_Q_Q3_2D_J2[4] = { 15, 14, 13, 12 };
static int C_Q_Q3_2D_J3[4] = { 12,  8,  4,  0 };

static int *C_Q_Q3_2D_J[4] = { C_Q_Q3_2D_J0, C_Q_Q3_2D_J1, 
                             C_Q_Q3_2D_J2, C_Q_Q3_2D_J3 };

// number of inner dofs
static int C_Q_Q3_2D_NInner = 4;

// array containing the numbers for the inner dofs
static int C_Q_Q3_2D_Inner[4] = { 5, 6, 9, 10 };

// number of outer dofs
static int C_Q_Q3_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int C_Q_Q3_2D_Outer[12] = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 };
