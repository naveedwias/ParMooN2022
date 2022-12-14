// ***********************************************************************
// Q3 element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_Q3_3D_NDOF = 40;

// number of dofs on the closure of the joints
static int N_H_Q3_3D_JointDOF = 6;

// which local dofs are on the joints
static int N_H_Q3_3D_J0[6] = { 0,  6, 12, 18, 24, 30 };
static int N_H_Q3_3D_J1[6] = { 1,  7, 13, 19, 25, 31 };
static int N_H_Q3_3D_J2[6] = { 2,  8, 14, 20, 26, 32 };
static int N_H_Q3_3D_J3[6] = { 3,  9, 15, 21, 27, 33 };
static int N_H_Q3_3D_J4[6] = { 4, 10, 16, 22, 28, 34 };
static int N_H_Q3_3D_J5[6] = { 5, 11, 17, 23, 29, 35 };

static int *N_H_Q3_3D_J[6] = { N_H_Q3_3D_J0, N_H_Q3_3D_J1,
                               N_H_Q3_3D_J2, N_H_Q3_3D_J3,
                               N_H_Q3_3D_J4, N_H_Q3_3D_J5};

// number of inner dofs
static int N_H_Q3_3D_NInner = 4;

// array containing the numbers for the inner dofs
static int N_H_Q3_3D_Inner[4] = { 36, 37, 38, 39 };
