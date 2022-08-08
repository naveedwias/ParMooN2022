// ***********************************************************************
// P2 discontinuous element, 1D
// ***********************************************************************

// number of degrees of freedom
static int D_L_P2_1D_NDOF = 3;

// number of dofs on the closure of the joints
static int D_L_P2_1D_JointDOF = 0;

// which local dofs are on the joints
static int *D_L_P2_1D_J0 = nullptr;
static int *D_L_P2_1D_J1 = nullptr;

static int *D_L_P2_1D_J[2] = { D_L_P2_1D_J0, D_L_P2_1D_J1 };

// number of inner dofs
static int D_L_P2_1D_NInner = 3;

// array containing the numbers for the inner dofs
static int D_L_P2_1D_Inner[3] = {0, 1, 2};
