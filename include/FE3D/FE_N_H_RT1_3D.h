// ***********************************************************************
// Raviart-Thomas element of first order on hexahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_RT1_3D_NDOF = 36;

// number of dofs on the closure of each joints
static int N_H_RT1_3D_JointDOF = 4;

// which local dofs are on the joints
static int N_H_RT1_3D_J0[4] = { 0, 1, 2, 3 };
static int N_H_RT1_3D_J1[4] = { 4, 5, 6, 7 };
static int N_H_RT1_3D_J2[4] = { 8, 9,10,11 };
static int N_H_RT1_3D_J3[4] = {12,13,14,15 };
static int N_H_RT1_3D_J4[4] = {16,17,18,19 };
static int N_H_RT1_3D_J5[4] = {20,21,22,23 };

static int *N_H_RT1_3D_J[6] = { N_H_RT1_3D_J0, N_H_RT1_3D_J1,
                                N_H_RT1_3D_J2, N_H_RT1_3D_J3,
                                N_H_RT1_3D_J4, N_H_RT1_3D_J5 };

// number of inner dofs
static int N_H_RT1_3D_NInner = 12;

// array containing the numbers for the inner dofs
static int N_H_RT1_3D_Inner[12] = { 24,25,26,27,28,29,30,31,32,33,34,35 };

// number of outer dofs
static int N_H_RT1_3D_NOuter = 24;

// array containing the numbers for the outer dofs
static int N_H_RT1_3D_Outer[24] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 };
