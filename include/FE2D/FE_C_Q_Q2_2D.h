// ***********************************************************************
// Q2 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q2_2D_NDOF = 9;

// number of dofs on the closure of the joints
static int C_Q_Q2_2D_JointDOF = 3;

// which local dofs are on the joints
static int C_Q_Q2_2D_J0[3] = { 0, 1, 2 };
static int C_Q_Q2_2D_J1[3] = { 2, 5, 8 };
static int C_Q_Q2_2D_J2[3] = { 8, 7, 6 };
static int C_Q_Q2_2D_J3[3] = { 6, 3, 0 };

static int *C_Q_Q2_2D_J[4] = { C_Q_Q2_2D_J0, C_Q_Q2_2D_J1,
                             C_Q_Q2_2D_J2, C_Q_Q2_2D_J3 };

// number of inner dofs
static int C_Q_Q2_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_Q_Q2_2D_Inner[1] = { 4 };

// number of outer dofs
static int C_Q_Q2_2D_NOuter = 8;

// array containing the numbers for the outer dofs
static int C_Q_Q2_2D_Outer[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };
