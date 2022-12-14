// ***********************************************************************
// Q1Rot element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_Q1_3D_NDOF = 6;

// number of dofs on the closure of the joints
static int N_H_Q1_3D_JointDOF = 1;

// which local dofs are on the joints
static int N_H_Q1_3D_J0[1] = { 0 };
static int N_H_Q1_3D_J1[1] = { 1 };
static int N_H_Q1_3D_J2[1] = { 2 };
static int N_H_Q1_3D_J3[1] = { 3 };
static int N_H_Q1_3D_J4[1] = { 4 };
static int N_H_Q1_3D_J5[1] = { 5 };

static int *N_H_Q1_3D_J[6] = { N_H_Q1_3D_J0, N_H_Q1_3D_J1,
                             N_H_Q1_3D_J2, N_H_Q1_3D_J3,
                             N_H_Q1_3D_J4, N_H_Q1_3D_J5};

// number of inner dofs
static int N_H_Q1_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_H_Q1_3D_Inner = nullptr;

// number of dofs on the closure of the edges
static int C_N_Q1_3D_EdgeDOF = 0;

// which local dofs are on the joints
static int *C_N_Q1_3D_E0 = nullptr;
static int *C_N_Q1_3D_E1 = nullptr;
static int *C_N_Q1_3D_E2 = nullptr;
static int *C_N_Q1_3D_E3 = nullptr;

static int *C_N_Q1_3D_E4 = nullptr;
static int *C_N_Q1_3D_E5 = nullptr;
static int *C_N_Q1_3D_E6 = nullptr;
static int *C_N_Q1_3D_E7 = nullptr;

static int *C_N_Q1_3D_E8 = nullptr;
static int *C_N_Q1_3D_E9 = nullptr;
static int *C_N_Q1_3D_E10 = nullptr;
static int *C_N_Q1_3D_E11 = nullptr;

static int *C_N_Q1_3D_E[12] = { C_N_Q1_3D_E0, C_N_Q1_3D_E1, C_N_Q1_3D_E2, C_N_Q1_3D_E3,
                                C_N_Q1_3D_E4, C_N_Q1_3D_E5, C_N_Q1_3D_E6, C_N_Q1_3D_E7,
                                C_N_Q1_3D_E8, C_N_Q1_3D_E9, C_N_Q1_3D_E10, C_N_Q1_3D_E11};

// number of dofs on the closure of the vertices
static int C_N_Q1_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *C_N_Q1_3D_Vert = nullptr;
