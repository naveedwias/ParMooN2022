// ***********************************************************************
// P2 element, with face and cell bubble conforming, 3D
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

// number of degrees of freedom
static int C_T_B2_3D_NDOF = 15;

// number of dofs on the closure of the joints
static int C_T_B2_3D_JointDOF = 7;

// which local dofs are on the joints
static int C_T_B2_3D_J0[7] = {  0,  1,  2,  3,  4,  5, 10 };
static int C_T_B2_3D_J1[7] = {  0,  6,  9,  1,  7,  2, 11 };
static int C_T_B2_3D_J2[7] = {  5,  4,  2,  8,  7,  9, 12 };
static int C_T_B2_3D_J3[7] = {  0,  3,  5,  6,  8,  9, 13 };

static int *C_T_B2_3D_J[4] = { C_T_B2_3D_J0, C_T_B2_3D_J1,
                             C_T_B2_3D_J2, C_T_B2_3D_J3 };

// number of dofs on the closure of the edges
static int C_T_B2_3D_EdgeDOF = 3;

// which local dofs are on the joints
static int C_T_B2_3D_E0[3] = { 0, 1, 2 };
static int C_T_B2_3D_E1[3] = { 2, 4, 5 };
static int C_T_B2_3D_E2[3] = { 5, 3, 0 };
static int C_T_B2_3D_E3[3] = { 0, 6, 9 };
static int C_T_B2_3D_E4[3] = { 2, 7, 9 };
static int C_T_B2_3D_E5[3] = { 5, 8, 9 };


static int *C_T_B2_3D_E[6] = { C_T_B2_3D_E0, C_T_B2_3D_E1, C_T_B2_3D_E2, C_T_B2_3D_E3,
                               C_T_B2_3D_E4, C_T_B2_3D_E5};

// number of dofs on the closure of the vertices
static int C_T_B2_3D_VertDOF = 1;

// array containing the numbers for the vertices dofs
static int C_T_B2_3D_Vert[4] =  {0, 2, 5, 9};

// number of inner dofs
static int C_T_B2_3D_NInner = 1;

// array containing the numbers for the inner dofs (here is no inner dof)
static int C_T_B2_3D_Inner[1] ={ 14 };
