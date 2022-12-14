// ***********************************************************************
// Q2 Raviart-Thomas vector element, nonconforming , 2D
// History:  05.12.2011 implementation (Ulrich)
// ***********************************************************************

// number of degrees of freedom
static int N_Q_RT2_2D_NDOF = 24;

// number of dofs on the closure of the joints
static int N_Q_RT2_2D_JointDOF = 3;

// which local dofs are on the joints
static int N_Q_RT2_2D_J0[3] = { 0, 1, 2};
static int N_Q_RT2_2D_J1[3] = { 3, 4, 5};
static int N_Q_RT2_2D_J2[3] = { 6, 7, 8};
static int N_Q_RT2_2D_J3[3] = { 9, 10, 11};

 
static int *N_Q_RT2_2D_J[4] = { N_Q_RT2_2D_J0,
                                N_Q_RT2_2D_J1,
                                N_Q_RT2_2D_J2,
                                N_Q_RT2_2D_J3
                              };

// number of inner dofs
static int N_Q_RT2_2D_NInner = 12;

// array containing the numbers for the inner dofs
static int N_Q_RT2_2D_Inner[12] = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

// number of outer dofs (dofs on edges)
static int N_Q_RT2_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_Q_RT2_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
