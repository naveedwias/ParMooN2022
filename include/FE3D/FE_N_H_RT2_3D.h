// ***********************************************************************
// Raviart-Thomas element of second order on hexahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_RT2_3D_NDOF = 108;

// number of dofs on the closure of each joints
static int N_H_RT2_3D_JointDOF = 9;

// which local dofs are on the joints
static int N_H_RT2_3D_J0[9] = {  0,  1,  2,  3,  4,  5,  6,  7,  8 };
static int N_H_RT2_3D_J1[9] = {  9, 10, 11, 12, 13, 14, 15, 16, 17 };
static int N_H_RT2_3D_J2[9] = { 18, 19, 20, 21, 22, 23, 24, 25, 26 };
static int N_H_RT2_3D_J3[9] = { 27, 28, 29, 30, 31, 32, 33, 34, 35 };
static int N_H_RT2_3D_J4[9] = { 36, 37, 38, 39, 40, 41, 42, 43, 44 };
static int N_H_RT2_3D_J5[9] = { 45, 46, 47, 48, 49, 50, 51, 52, 53 };

static int *N_H_RT2_3D_J[6] = { N_H_RT2_3D_J0, N_H_RT2_3D_J1,
                                N_H_RT2_3D_J2, N_H_RT2_3D_J3,
                                N_H_RT2_3D_J4, N_H_RT2_3D_J5 };

// number of inner dofs
static int N_H_RT2_3D_NInner = 54;

// array containing the numbers for the inner dofs
static int N_H_RT2_3D_Inner[54] = {
  54,55,56,57,58,59,60,61,62,63,64,65,66,67,
  68,69,70,71,72,73,74,75,76,77,78,79,80,81,
  82,83,84,85,86,87,88,89,90,91,92,93,94,95,
  96,97,98,99,100,101,102,103,104,105,106,107};

// number of outer dofs
static int N_H_RT2_3D_NOuter = 54;

// array containing the numbers for the outer dofs
static int N_H_RT2_3D_Outer[54] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
  18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
  33,34,35,36,37,38,39,40,41,42,43,44,45,
  46,47,48,49,50,51,52,53};
