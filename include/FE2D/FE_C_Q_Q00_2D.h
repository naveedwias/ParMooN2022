// ***********************************************************************
// Q00 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q00_2D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_Q_Q00_2D_JointDOF = 0;

// which local dofs are on the joints
static int *C_Q_Q00_2D_J[3] = { nullptr, nullptr, nullptr };

// number of inner dofs
static int C_Q_Q00_2D_NInner = 1;

// array containing the numbers for the inner dofs 
static int C_Q_Q00_2D_Inner[1] = { 0 };
