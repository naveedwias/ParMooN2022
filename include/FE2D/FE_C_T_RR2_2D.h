// red refined triangle, second order, used for LPS

// ***********************************************************************
// P2 element, conforming, 2D, on a subdivision into four triangles
// ***********************************************************************

// number of degrees of freedom
static int C_T_RR2_2D_NDOF = 15;

// number fo dofs on the closure of joints
static int C_T_RR2_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_T_RR2_2D_J0[5]={ 0, 1,  2,  3, 4 };
static int C_T_RR2_2D_J1[5]={ 4, 5,  6,  7, 8 };
static int C_T_RR2_2D_J2[5]={ 8, 9, 10, 11, 0 };

static int *C_T_RR2_2D_J[3]={ C_T_RR2_2D_J0, C_T_RR2_2D_J1,
                              C_T_RR2_2D_J2 };

// number of inner dofs
static int C_T_RR2_2D_NInner = 3;

// array containing the numbers for the inner dofs
static int C_T_RR2_2D_Inner[3] = { 12, 13, 14 };
