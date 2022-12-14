#include "Enumerations_fe.h"

/**************************** NSTYPE == 1 ******************************/
// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one A block, 
//      B1, B2, B3 (divergence blocks)
// ======================================================================

int NSType1N_Terms = 5;
MultiIndex3D NSType1Derivatives[5] = { D100, D010, D001, D000, D000 };
int NSType1SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
int NSType1N_Matrices = 4;
int NSType1RowSpace[4] = { 0, 1, 1, 1 };
int NSType1ColumnSpace[4] = { 0, 0, 0, 0 };
int NSType1N_Rhs = 3;
int NSType1RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one nonlinear A block
//      WITHOUT right hand sides
// ======================================================================

int NSType1NLN_Terms = 4;
MultiIndex3D NSType1NLDerivatives[4] = { D100, D010, D001, D000 };
int NSType1NLSpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType1NLN_Matrices = 1;
int NSType1NLRowSpace[1] = { 0 };
int NSType1NLColumnSpace[1] = { 0 };
int NSType1NLN_Rhs = 0;
int *NSType1NLRhsSpace = nullptr;

int NSType1NLSDFEMN_Rhs = 3;
int NSType1NLSDFEMRhsSpace[3] = {0, 0, 0};

// ======================================================================
// VMSProjection
// ======================================================================
int NSType1VMSProjectionN_Terms = 6;
MultiIndex3D NSType1VMSProjectionDerivatives[6] = { D100, D010, D001, D000, D000, D000 };
int NSType1VMSProjectionSpaceNumbers[6] = { 0, 0, 0, 0, 1, 2 };
int NSType1VMSProjectionN_Matrices = 11;
int NSType1VMSProjectionRowSpace[11] = { 0, 2, 1, 1, 1, 
					 0, 0, 0, 2, 2, 2};
int NSType1VMSProjectionColumnSpace[11] = { 0, 2, 0, 0, 0,
					     2, 2, 2, 0, 0, 0};

int NSType1_2NLVMSProjectionN_Terms = 5;
MultiIndex3D NSType1_2NLVMSProjectionDerivatives[5] = { D100, D010, D001, D000, D000 };
int NSType1_2NLVMSProjectionSpaceNumbers[5] = { 0, 0, 0, 0, 1 };
int NSType1_2NLVMSProjectionN_Matrices = 4;
int NSType1_2NLVMSProjectionRowSpace[4] = { 0, 0, 0, 0};
int NSType1_2NLVMSProjectionColumnSpace[4] = { 0, 1, 1, 1};

/**************************** NSTYPE == 2 ******************************/

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      B1, B2, B3 (divergence blocks), 
//      B1T, B2T, B3T (gradient blocks)
// ======================================================================

int NSType2N_Terms = 5;
MultiIndex3D NSType2Derivatives[5] = { D100, D010, D001, D000, D000 };
int NSType2SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
int NSType2N_Matrices = 7;
int NSType2RowSpace[7] = { 0, 1, 1, 1, 0, 0, 0 };
int NSType2ColumnSpace[7] = { 0, 0, 0, 0, 1, 1, 1};
int NSType2N_Rhs = 3;
int NSType2RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITH B1, B2, B3,  B1T, B2T, B3T
//      WITH right hand sides
// ======================================================================

int NSType2SDN_Terms = 8;
MultiIndex3D NSType2SDDerivatives[8] = { D100, D010, D001, D000, 
                                         D100, D010, D001, D000};
int NSType2SDSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1  };
int NSType2SDN_Matrices = 7;
int NSType2SDRowSpace[7] = { 0, 1, 1, 1, 0, 0, 0 };
int NSType2SDColumnSpace[7] = { 0, 0, 0, 0, 1, 1, 1 };
int NSType2SDN_Rhs = 3;
int NSType2SDRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int NSType2NLN_Terms = 4;
MultiIndex3D NSType2NLDerivatives[4] = { D100, D010, D001, D000 };
int NSType2NLSpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType2NLN_Matrices = 1;
int NSType2NLRowSpace[1] = { 0 };
int NSType2NLColumnSpace[1] = { 0 };
int NSType2NLN_Rhs = 0;
int *NSType2NLRhsSpace = nullptr;


// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITH B1T, B2T, B3T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType2NLSDN_Terms = 8;
MultiIndex3D NSType2NLSDDerivatives[8] = { D100, D010, D001, D000, 
                                         D100, D010, D001, D000  };
int NSType2NLSDSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1  };
int NSType2NLSDN_Matrices = 4;
int NSType2NLSDRowSpace[4] = { 0, 0, 0, 0 };
int NSType2NLSDColumnSpace[4] = { 0, 1, 1, 1 };
int NSType2NLSDN_Rhs = 3;
int NSType2NLSDRhsSpace[3] = { 0, 0, 0 };

/**************************** NSTYPE == 3 ******************************/

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

int NSType3N_Terms = 5;
MultiIndex3D NSType3Derivatives[5] = { D100, D010, D001, D000, D000 };
int NSType3SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
int NSType3N_Matrices = 12;
int NSType3RowSpace[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
int NSType3ColumnSpace[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType3N_Rhs = 3;
int NSType3RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A22, A33
//      WITHOUT right hand sides
// ======================================================================

int NSType3NLN_Terms = 4;
MultiIndex3D NSType3NLDerivatives[4] = { D100, D010, D001, D000 };
int NSType3NLSpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType3NLN_Matrices = 3;
int NSType3NLRowSpace[3] = { 0, 0, 0 };
int NSType3NLColumnSpace[3] = { 0, 0, 0 };
int NSType3NLN_Rhs = 0;
int *NSType3NLRhsSpace = nullptr;

int NSType3_4NLNewtonN_Matrices = 9;
int NSType3_4NLNewtonRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType3_4NLNewtonColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
int NSType3_4NLNewtonN_Rhs = 3;
int NSType3_4NLNewtonRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all blocks Aij
//      WITHOUT right hand sides
// ======================================================================

int NSType3NLSmagorinskyN_Terms = 4;
MultiIndex3D NSType3NLSmagorinskyDerivatives[4] = { D100, D010, D001, D000 };
int NSType3NLSmagorinskySpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType3NLSmagorinskyN_Matrices = 9;
int NSType3NLSmagorinskyRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType3NLSmagorinskyColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType3NLSmagorinskyN_Rhs = 0;
int *NSType3NLSmagorinskyRhsSpace = nullptr;

/**************************** NSTYPE == 4 ******************************/

/** assemble all matrices and rhs */

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all nine A blocks,
//      B1, B2, B3 (divergence blocks), 
//      B1T, B2T, B3T (gradient blocks)
// ======================================================================

int NSType4N_Terms = 5;
MultiIndex3D NSType4Derivatives[5] = { D100, D010, D001, D000, D000 };
int NSType4SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
int NSType4N_Matrices = 15;
int NSType4RowSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
int NSType4ColumnSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
int NSType4N_Rhs = 3;
int NSType4RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int NSType4SDN_Terms = 8;
MultiIndex3D NSType4SDDerivatives[8] = { D100, D010, D001, D000, 
                                         D100, D010, D001, D000 };
int NSType4SDSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1  };
int NSType4SDN_Matrices = 15;
int NSType4SDRowSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
int NSType4SDColumnSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
int NSType4SDN_Rhs = 3;
int NSType4SDRhsSpace[3] = { 0, 0, 0 };

/** assemble some matrices and no rhs */

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A22
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int NSType4NLN_Terms = 4;
MultiIndex3D NSType4NLDerivatives[4] = { D100, D010, D001, D000 };
int NSType4NLSpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType4NLN_Matrices = 3;
int NSType4NLRowSpace[3] = { 0, 0, 0 };
int NSType4NLColumnSpace[3] = { 0, 0, 0 };
int NSType4NLN_Rhs = 0;
int *NSType4NLRhsSpace = nullptr;

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType4NLSDN_Terms = 8;
MultiIndex3D NSType4NLSDDerivatives[8] = { D100, D010, D001, D000, 
                                           D100, D010, D001, D000  };
int NSType4NLSDSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1  };
int NSType4NLSDN_Matrices = 6;
int NSType4NLSDRowSpace[6] = { 0, 0, 0, 0, 0, 0 };
int NSType4NLSDColumnSpace[6] = { 0, 0, 0, 1, 1, 1 };
int NSType4NLSDN_Rhs = 3;
int NSType4NLSDRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType4NLSD_DivDiv_N_Matrices = 12;
int NSType4NLSD_DivDiv_RowSpace[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType4NLSD_DivDiv_ColumnSpace[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      without right hand sides
// ======================================================================

int NSType4NLSmagorinskyN_Terms = 4;
MultiIndex3D NSType4NLSmagorinskyDerivatives[4] = { D100, D010, D001, D000 };
int NSType4NLSmagorinskySpaceNumbers[4] = { 0, 0, 0, 0 };
int NSType4NLSmagorinskyN_Matrices = 9;
int NSType4NLSmagorinskyRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType4NLSmagorinskyColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int NSType4NLSmagorinskyN_Rhs = 0;
int *NSType4NLSmagorinskyRhsSpace = nullptr;

// ======================================================================
// declaration for pressure separation
//    only rhs
// ======================================================================

int NSPressSepN_Terms = 1;
MultiIndex3D NSPressSepDerivatives[1] = { D000 };
int NSPressSepSpaceNumbers[1] = { 0 };
int NSPressSepN_Matrices = 0;
int *NSPressSepRowSpace = nullptr;
int *NSPressSepColumnSpace = nullptr;
int NSPressSepN_Rhs = 3;
int NSPressSepRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// declaration for pressure separation
// with auxiliary problem
// information for ansatz and test functions
// ======================================================================
int NSPressSepAuxProbN_Terms = 3;
MultiIndex3D NSPressSepAuxProbDerivatives[3] = { D100, D010, D001 };
int NSPressSepAuxProbSpaceNumbers[3] = { 0, 0, 0 };
int NSPressSepAuxProbN_Matrices = 1;
int NSPressSepAuxProbRowSpace[1] = { 0 };
int NSPressSepAuxProbColumnSpace[1] = { 0 };
int NSPressSepAuxProbN_Rhs = 1;
int NSPressSepAuxProbRhsSpace[1] = { 0 };

// ======================================================================
// declaration for computation of rhs for RFB
//    only rhs
// ======================================================================

int NSRFBRhsN_Terms = 1;
MultiIndex3D NSRFBRhsDerivatives[1] = { D000 };
int NSRFBRhsSpaceNumbers[1] = { 0 };
int NSRFBRhsN_Matrices = 0;
int *NSRFBRhsRowSpace = nullptr;
int *NSRFBRhsColumnSpace = nullptr;
int NSRFBRhsN_Rhs = 3;
int NSRFBRhsRhsSpace[3] = { 0, 0, 0 };
