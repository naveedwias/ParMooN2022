// ***********************************************************************
// Q5 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_Q_Q5_2D_NDOF = 30;

// number of dofs on the closure of the joints
static int N_Q_Q5_2D_JointDOF = 5;

// which local dofs are on the joints
static int N_Q_Q5_2D_J0[5] = { 0, 4,  8, 12, 16 };
static int N_Q_Q5_2D_J1[5] = { 1, 5,  9, 13, 17 };
static int N_Q_Q5_2D_J2[5] = { 2, 6, 10, 14, 18 };
static int N_Q_Q5_2D_J3[5] = { 3, 7, 11, 15, 19 };
 
static int *N_Q_Q5_2D_J[4] = { N_Q_Q5_2D_J0, N_Q_Q5_2D_J1,
                                 N_Q_Q5_2D_J2, N_Q_Q5_2D_J3 };

// number of inner dofs
static int N_Q_Q5_2D_NInner = 10;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_Q_Q5_2D_Inner[10] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };

// number of outer dofs
static int N_Q_Q5_2D_NOuter = 20;

// array containing the numbers for the outer dofs
static int N_Q_Q5_2D_Outer[20] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                   12, 13, 14, 15, 16, 17, 18, 19 };
