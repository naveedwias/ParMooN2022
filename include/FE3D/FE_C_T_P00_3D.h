// ***********************************************************************
// P00 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P00_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_T_P00_3D_JointDOF = 0;

// which local dofs are on the joints
static int *C_T_P00_3D_J0 = nullptr;
static int *C_T_P00_3D_J1 = nullptr;
static int *C_T_P00_3D_J2 = nullptr;
static int *C_T_P00_3D_J3 = nullptr;

static int *C_T_P00_3D_J[4] = { C_T_P00_3D_J0, C_T_P00_3D_J1,
                               C_T_P00_3D_J2, C_T_P00_3D_J3 };

// number of inner dofs
static int C_T_P00_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_T_P00_3D_Inner[1] = { 0 };
