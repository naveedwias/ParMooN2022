// ***********************************************************************
// Q1 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q1_3D_NDOF = 8;

// number of dofs on the closure of the joints
static int C_H_Q1_3D_JointDOF = 4;

// which local dofs are on the joints
static int C_H_Q1_3D_J0[4] = { 0, 1, 2, 3 };
static int C_H_Q1_3D_J1[4] = { 0, 4, 1, 5 };
static int C_H_Q1_3D_J2[4] = { 1, 5, 3, 7 };
static int C_H_Q1_3D_J3[4] = { 3, 7, 2, 6 };
static int C_H_Q1_3D_J4[4] = { 0, 2, 4, 6 };
static int C_H_Q1_3D_J5[4] = { 4, 6, 5, 7 };

static int *C_H_Q1_3D_J[6] = { C_H_Q1_3D_J0, C_H_Q1_3D_J1,
                             C_H_Q1_3D_J2, C_H_Q1_3D_J3,
                             C_H_Q1_3D_J4, C_H_Q1_3D_J5};

// number of dofs on the closure of the edges
static int C_H_Q1_3D_EdgeDOF = 2;

// which local dofs are on the joints
static int C_H_Q1_3D_E0[2] = { 0, 1 };
static int C_H_Q1_3D_E1[2] = { 1, 3 };
static int C_H_Q1_3D_E2[2] = { 3, 2 };
static int C_H_Q1_3D_E3[2] = { 2, 0 };

static int C_H_Q1_3D_E4[2] = { 0, 4 };
static int C_H_Q1_3D_E5[2] = { 1, 5 };
static int C_H_Q1_3D_E6[2] = { 3, 7 };
static int C_H_Q1_3D_E7[2] = { 2, 6 };

static int C_H_Q1_3D_E8[2] = { 4, 5 };
static int C_H_Q1_3D_E9[2] = { 5, 7 };
static int C_H_Q1_3D_E10[2] = { 7, 6 };
static int C_H_Q1_3D_E11[2] = { 6, 4 };

static int *C_H_Q1_3D_E[12] = { C_H_Q1_3D_E0, C_H_Q1_3D_E1, C_H_Q1_3D_E2, C_H_Q1_3D_E3,
                                C_H_Q1_3D_E4, C_H_Q1_3D_E5, C_H_Q1_3D_E6, C_H_Q1_3D_E7,
                                C_H_Q1_3D_E8, C_H_Q1_3D_E9, C_H_Q1_3D_E10, C_H_Q1_3D_E11};

// number of dofs on the closure of the vertices
static int C_H_Q1_3D_VertDOF = 1;

// array containing the numbers for the vertices dofs
static int C_H_Q1_3D_Vert[8] =  {0, 1, 3, 2, 4, 5, 7, 6};

// number of inner dofs
static int C_H_Q1_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *C_H_Q1_3D_Inner = nullptr;
