// ***********************************************************************
// P1 element, discontinous, 2D, triangle
// ***********************************************************************

// number of degrees of freedom
static int D_T_P1_2D_NDOF = 3;

// number of dofs on the closure of the joints
static int D_T_P1_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_P1_2D_J0 = nullptr;
static int *D_T_P1_2D_J1 = nullptr;
static int *D_T_P1_2D_J2 = nullptr;

static int *D_T_P1_2D_J[3] = { D_T_P1_2D_J0, D_T_P1_2D_J1, D_T_P1_2D_J2 };

// number of inner dofs
static int D_T_P1_2D_NInner = 3;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_T_P1_2D_Inner[3] = { 0, 1, 2 };
