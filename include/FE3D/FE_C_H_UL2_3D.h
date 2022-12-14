// ***********************************************************************
// Q2 element with bubbles, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_UL2_3D_NDOF = 30;

// number of dofs on the closure of the joints
static int C_H_UL2_3D_JointDOF = 9;

// which local dofs are on the joints
static int C_H_UL2_3D_J0[9] = {  0,  1,  2,  3,  4,  5,  6,  7,  8 };
static int C_H_UL2_3D_J1[9] = {  0,  9, 17,  1, 10, 18,  2, 11, 19 };
static int C_H_UL2_3D_J2[9] = {  2, 11, 19,  5, 12, 22,  8, 13, 25 };
static int C_H_UL2_3D_J3[9] = {  8, 13, 25,  7, 14, 24,  6, 15, 23 };
static int C_H_UL2_3D_J4[9] = {  0,  3,  6,  9, 16, 15, 17, 20, 23 };
static int C_H_UL2_3D_J5[9] = { 17, 20, 23, 18, 21, 24, 19, 22, 25 };

static int *C_H_UL2_3D_J[6] = { C_H_UL2_3D_J0, C_H_UL2_3D_J1,
                             C_H_UL2_3D_J2, C_H_UL2_3D_J3,
                             C_H_UL2_3D_J4, C_H_UL2_3D_J5};

// number of inner dofs
static int C_H_UL2_3D_NInner = 4;

// array containing the numbers for the inner dofs
static int C_H_UL2_3D_Inner[4] = { 26, 27, 28, 29 };

// number of outer dof
static int C_H_UL2_3D_NOuter = 26;

// array containing the numbers for the outer dofs
static int C_H_UL2_3D_Outer[26] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,
                                     9, 10, 11, 12, 13, 14, 15, 16, 17,
                                    18, 19, 20, 21, 22, 23, 24, 25 };
