// ***********************************************************************
// P0 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P0_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_T_P0_3D_JointDOF = 0;

// which local dofs are on the joints
static int *C_T_P0_3D_J0 = nullptr;
static int *C_T_P0_3D_J1 = nullptr;
static int *C_T_P0_3D_J2 = nullptr;
static int *C_T_P0_3D_J3 = nullptr;

static int *C_T_P0_3D_J[4] = { C_T_P0_3D_J0, C_T_P0_3D_J1,
                               C_T_P0_3D_J2, C_T_P0_3D_J3 };

// number of inner dofs
static int C_T_P0_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_T_P0_3D_Inner[1] = { 0 };

// number of dofs on the closure of the edges
static int  C_T_P0_3D_EdgeDOF = 0;

// which local dofs are on the joints
static int *C_T_P0_3D_E0 = nullptr;
static int *C_T_P0_3D_E1 = nullptr;
static int *C_T_P0_3D_E2 = nullptr;
static int *C_T_P0_3D_E3 = nullptr;
static int *C_T_P0_3D_E4 = nullptr;
static int *C_T_P0_3D_E5 = nullptr;

static int *C_T_P0_3D_E[6] = { C_T_P0_3D_E0, C_T_P0_3D_E1, C_T_P0_3D_E2, C_T_P0_3D_E3,
                               C_T_P0_3D_E4, C_T_P0_3D_E5};
			       
// number of dofs on the closure of the vertices
static int C_T_P0_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *C_T_P0_3D_Vert =  nullptr;
