// ***********************************************************************
// Raviart-Thomas element of first order on tetrahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_RT1_3D_NDOF = 15;

// number of dofs on the closure of each joints
static int N_T_RT1_3D_JointDOF = 3;

// which local dofs are on the joints
static int N_T_RT1_3D_J0[3] = { 0, 1, 2 };
static int N_T_RT1_3D_J1[3] = { 3, 4, 5 };
static int N_T_RT1_3D_J2[3] = { 6, 7, 8 };
static int N_T_RT1_3D_J3[3] = { 9,10,11 };

static int *N_T_RT1_3D_J[4] = { N_T_RT1_3D_J0, N_T_RT1_3D_J1,
                                N_T_RT1_3D_J2, N_T_RT1_3D_J3 };

// number of inner dofs
static int N_T_RT1_3D_NInner = 3;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_T_RT1_3D_Inner[3] = { 12,13,14 };

// number of outer dofs
static int N_T_RT1_3D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_T_RT1_3D_Outer[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
