// ***********************************************************************
// Q2 element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_Q2_3D_NDOF = 19;

// number of dofs on the closure of the joints
static int N_H_Q2_3D_JointDOF = 3;

// which local dofs are on the joints
static int N_H_Q2_3D_J0[3] = { 0,  6, 12 };
static int N_H_Q2_3D_J1[3] = { 1,  7, 13 };
static int N_H_Q2_3D_J2[3] = { 2,  8, 14 };
static int N_H_Q2_3D_J3[3] = { 3,  9, 15 };
static int N_H_Q2_3D_J4[3] = { 4, 10, 16 };
static int N_H_Q2_3D_J5[3] = { 5, 11, 17 };

static int *N_H_Q2_3D_J[6] = { N_H_Q2_3D_J0, N_H_Q2_3D_J1,
                             N_H_Q2_3D_J2, N_H_Q2_3D_J3,
                             N_H_Q2_3D_J4, N_H_Q2_3D_J5};

// number of inner dofs
static int N_H_Q2_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int N_H_Q2_3D_Inner[1] = { 18 };

// number of outer dofs
static int N_H_Q2_3D_NOuter = 18;

// array containing the numbers for the outer dofs
static int N_H_Q2_3D_Outer[18] = {  0,  1,  2,  3,  4,  5,
                                    6,  7,  8,  9, 10, 11,
                                   12, 13, 14, 15, 16, 17 };
