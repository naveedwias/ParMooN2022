// ***********************************************************************
// Raviart-Thomas element of third order on tetrahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_RT3_3D_NDOF = 70;

// number of dofs on the closure of each joints
static int N_T_RT3_3D_JointDOF = 10;

// which local dofs are on the joints
static int N_T_RT3_3D_J0[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static int N_T_RT3_3D_J1[10] = {10,11,12,13,14,15,16,17,18,19 };
static int N_T_RT3_3D_J2[10] = {20,21,22,23,24,25,26,27,28,29 };
static int N_T_RT3_3D_J3[10] = {30,31,32,33,34,35,36,37,38,39 };

static int *N_T_RT3_3D_J[4] = { N_T_RT3_3D_J0, N_T_RT3_3D_J1,
                                N_T_RT3_3D_J2, N_T_RT3_3D_J3 };

// number of inner dofs
static int N_T_RT3_3D_NInner = 30;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_T_RT3_3D_Inner[30] = { 40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69 };

// number of outer dofs
static int N_T_RT3_3D_NOuter = 40;

// array containing the numbers for the outer dofs
static int N_T_RT3_3D_Outer[40] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 };
