// ***********************************************************************
// BDM vector element, nonconforming , 2D
// History:  23.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// number of degrees of freedom
static int N_Q_BDM3_2D_NDOF = 22;

// number of dofs on the closure of the joints
static int N_Q_BDM3_2D_JointDOF = 4;

// which local dofs are on the joints
static int N_Q_BDM3_2D_J0[4] = { 0, 1, 2, 3 };
static int N_Q_BDM3_2D_J1[4] = { 4, 5, 6, 7 };
static int N_Q_BDM3_2D_J2[4] = { 8, 9, 10,11};
static int N_Q_BDM3_2D_J3[4] = { 12,13,14,15};

 
static int *N_Q_BDM3_2D_J[4] = { N_Q_BDM3_2D_J0,
                                N_Q_BDM3_2D_J1,
                                N_Q_BDM3_2D_J2,
                                N_Q_BDM3_2D_J3
                              };

// number of inner dofs
static int N_Q_BDM3_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int N_Q_BDM3_2D_Inner[6] = {16,17,18,19,20,21};

// number of outer dofs (dofs on edges)
static int N_Q_BDM3_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int N_Q_BDM3_2D_Outer[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
