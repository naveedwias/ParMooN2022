// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// number of degrees of freedom
static int B_Q_IB2_2D_NDOF = 1;

// number of dofs on the closure of the joints
static int B_Q_IB2_2D_JointDOF = 0;

// which local dofs are on the joints
static int *B_Q_IB2_2D_J[3] = { nullptr, nullptr, nullptr };

// number of inner dofs
static int B_Q_IB2_2D_NInner = 1;

// array containing the numbers for the inner dofs 
static int B_Q_IB2_2D_Inner[1] = { 0 };
