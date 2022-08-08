// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// number of degrees of freedom
static int B_H_IB2_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int B_H_IB2_3D_JointDOF = 0;

// which local dofs are on the joints
static int *B_H_IB2_3D_J0 = nullptr;
static int *B_H_IB2_3D_J1 = nullptr;
static int *B_H_IB2_3D_J2 = nullptr;
static int *B_H_IB2_3D_J3 = nullptr;
static int *B_H_IB2_3D_J4 = nullptr;
static int *B_H_IB2_3D_J5 = nullptr;

static int *B_H_IB2_3D_J[6] = { B_H_IB2_3D_J0, B_H_IB2_3D_J1,
                               B_H_IB2_3D_J2, B_H_IB2_3D_J3,
                               B_H_IB2_3D_J4, B_H_IB2_3D_J5};

// number of inner dofs
static int B_H_IB2_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int B_H_IB2_3D_Inner[1] = { 0 };
