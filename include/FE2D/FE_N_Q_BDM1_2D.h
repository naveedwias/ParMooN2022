// ***********************************************************************
// Q1 BDM vector element, nonconforming , 2D
// History:  17.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// number of degrees of freedom
static int N_Q_BDM1_2D_NDOF = 8;

// number of dofs on the closure of the joints
static int N_Q_BDM1_2D_JointDOF = 2;

// which local dofs are on the joints
static int N_Q_BDM1_2D_J0[2] = { 0, 1 };
static int N_Q_BDM1_2D_J1[2] = { 2, 3 };
static int N_Q_BDM1_2D_J2[2] = { 4, 5 };
static int N_Q_BDM1_2D_J3[2] = { 6, 7 };
 
static int *N_Q_BDM1_2D_J[4] = { N_Q_BDM1_2D_J0, N_Q_BDM1_2D_J1,
                                 N_Q_BDM1_2D_J2, N_Q_BDM1_2D_J3 };
// number of inner dofs
static int N_Q_BDM1_2D_NInner = 0;

// array containing the numbers for the inner dofs 
static int *N_Q_BDM1_2D_Inner = nullptr;

// number of outer dofs
static int N_Q_BDM1_2D_NOuter = 8;

// array containing the numbers for the outer dofs
static int N_Q_BDM1_2D_Outer[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
