// ***********************************************************************
// P6 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// number of degrees of freedom
static int D_Q_P6_2D_NDOF = 28;

// number of dofs on the closure of the joints
static int D_Q_P6_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_Q_P6_2D_J0 = nullptr;
static int *D_Q_P6_2D_J1 = nullptr;
static int *D_Q_P6_2D_J2 = nullptr;
static int *D_Q_P6_2D_J3 = nullptr;

static int *D_Q_P6_2D_J[4] = { D_Q_P6_2D_J0, D_Q_P6_2D_J1,
                               D_Q_P6_2D_J2, D_Q_P6_2D_J3 };

// number of inner dofs
static int D_Q_P6_2D_NInner = 28;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_Q_P6_2D_Inner[28] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                   20, 21, 22, 23, 24, 25, 26, 27 };
