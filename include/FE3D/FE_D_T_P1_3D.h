// ***********************************************************************
// P1 element, discontinuous, 3D
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

// number of degrees of freedom
static int D_T_P1_3D_NDOF = 4;

// number of dofs on the closure of the joints
static int D_T_P1_3D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_P1_3D_J0 = nullptr;
static int *D_T_P1_3D_J1 = nullptr;
static int *D_T_P1_3D_J2 = nullptr;
static int *D_T_P1_3D_J3 = nullptr;

static int *D_T_P1_3D_J[4] = { D_T_P1_3D_J0, D_T_P1_3D_J1,
                               D_T_P1_3D_J2, D_T_P1_3D_J3 };

// number of dofs on the closure of the edges
static int D_T_P1_3D_EdgeDOF = 0;

static int *D_T_P1_3D_E0 = nullptr;
static int *D_T_P1_3D_E1 = nullptr;
static int *D_T_P1_3D_E2 = nullptr;
static int *D_T_P1_3D_E3 = nullptr;
static int *D_T_P1_3D_E4 = nullptr;
static int *D_T_P1_3D_E5 = nullptr;

static int *D_T_P1_3D_E[6] = {D_T_P1_3D_E0, D_T_P1_3D_E1, D_T_P1_3D_E2, D_T_P1_3D_E3,
                               D_T_P1_3D_E4, D_T_P1_3D_E5};


// number of dofs on the closure of the vertex
static int D_T_P1_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *D_T_P1_3D_Vert = nullptr;


// number of inner dofs
static int D_T_P1_3D_NInner = 4;

// array containing the numbers for the inner dofs 
static int D_T_P1_3D_Inner[] = { 0, 1, 2, 3 };

