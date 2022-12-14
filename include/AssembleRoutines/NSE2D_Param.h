#include "Enumerations_fe.h"

/**************************** NSTYPE == 1 ******************************/
// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one A block, 
//      B1, B2 (divergence blocks)
// ======================================================================

int NSType1N_Terms = 4;
MultiIndex2D NSType1Derivatives[4] = { D10, D01, D00, D00 };
int NSType1SpaceNumbers[4] = { 0, 0, 0, 1 };
int NSType1N_Matrices = 3;
int NSType1RowSpace[3] = { 0, 1, 1 };
int NSType1ColumnSpace[3] = { 0, 0, 0 };
int NSType1N_Rhs = 2;
int NSType1RhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one nonlinear A block
//      WITHOUT right hand sides
// ======================================================================

int NSType1NLN_Terms = 3;
MultiIndex2D NSType1NLDerivatives[3] = { D10, D01, D00 };
int NSType1NLSpaceNumbers[3] = { 0, 0, 0 };
int NSType1NLN_Matrices = 1;
int NSType1NLRowSpace[1] = { 0 };
int NSType1NLColumnSpace[1] = { 0 };
int NSType1NLN_Rhs = 0;
int *NSType1NLRhsSpace = nullptr;

int NSType1NLSDFEMN_Rhs = 2;
int NSType1NLSDFEMRhsSpace[2] = {0, 0};

// ======================================================================
// VMSProjection
// ======================================================================
int NSType1VMSProjectionN_Terms = 5;
MultiIndex2D NSType1VMSProjectionDerivatives[5] = { D10, D01, D00, D00, D00 };
int NSType1VMSProjectionSpaceNumbers[5] = { 0, 0, 0, 1, 2 };
int NSType1VMSProjectionN_Matrices = 8;
int NSType1VMSProjectionRowSpace[8] = { 0, 2, 1, 1, 
					 0, 0, 2, 2};
int NSType1VMSProjectionColumnSpace[8] = { 0, 2, 0, 0, 
					     2, 2, 0, 0};

int NSType1_2NLVMSProjectionN_Terms = 4;
MultiIndex2D NSType1_2NLVMSProjectionDerivatives[4] = { D10, D01, D00, D00 };
int NSType1_2NLVMSProjectionSpaceNumbers[4] = { 0, 0, 0, 1 };
int NSType1_2NLVMSProjectionN_Matrices = 3;
int NSType1_2NLVMSProjectionRowSpace[3] = { 0, 0, 0};
int NSType1_2NLVMSProjectionColumnSpace[3] = { 0, 1, 1};

/**************************** NSTYPE == 2 ******************************/

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int NSType2N_Terms = 4;
MultiIndex2D NSType2Derivatives[4] = { D10, D01, D00, D00 };
int NSType2SpaceNumbers[4] = { 0, 0, 0, 1 };
int NSType2N_Matrices = 5;
int NSType2RowSpace[5] = { 0, 1, 1, 0, 0 };
int NSType2ColumnSpace[5] = { 0, 0, 0, 1, 1 };
int NSType2N_Rhs = 2;
int NSType2RhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITH B1, B2 B1T, B2T
//      WITH right hand sides
// ======================================================================

int NSType2SDN_Terms = 8;
MultiIndex2D NSType2SDDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType2SDSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType2SDN_Matrices = 5;
int NSType2SDRowSpace[5] = { 0, 1, 1, 0, 0 };
int NSType2SDColumnSpace[5] = { 0, 0, 0, 1, 1 };
int NSType2SDN_Rhs = 2;
int NSType2SDRhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int NSType2NLN_Terms = 3;
MultiIndex2D NSType2NLDerivatives[3] = { D10, D01, D00 };
int NSType2NLSpaceNumbers[3] = { 0, 0, 0 };
int NSType2NLN_Matrices = 1;
int NSType2NLRowSpace[1] = { 0 };
int NSType2NLColumnSpace[1] = { 0 };
int NSType2NLN_Rhs = 0;
int *NSType2NLRhsSpace = nullptr;


// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType2NLSDN_Terms = 8;
MultiIndex2D NSType2NLSDDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType2NLSDSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType2NLSDN_Matrices = 3;
int NSType2NLSDRowSpace[3] = { 0, 0, 0 };
int NSType2NLSDColumnSpace[3] = { 0, 1, 1 };
int NSType2NLSDN_Rhs = 2;
int NSType2NLSDRhsSpace[2] = { 0, 0 };

/**************************** NSTYPE == 3 ******************************/

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

int NSType3N_Terms = 4;
MultiIndex2D NSType3Derivatives[4] = { D10, D01, D00, D00 };
int NSType3SpaceNumbers[4] = { 0, 0, 0, 1 };
int NSType3N_Matrices = 6;
int NSType3RowSpace[6] = { 0, 0, 0, 0, 1, 1 };
int NSType3ColumnSpace[6] = { 0, 0, 0, 0, 0, 0 };
int NSType3N_Rhs = 2;
int NSType3RhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A22
//      WITHOUT right hand sides
// ======================================================================

int NSType3NLN_Terms = 3;
MultiIndex2D NSType3NLDerivatives[3] = { D10, D01, D00 };
int NSType3NLSpaceNumbers[3] = { 0, 0, 0 };
int NSType3NLN_Matrices = 2;
int NSType3NLRowSpace[2] = { 0, 0 };
int NSType3NLColumnSpace[2] = { 0, 0 };
int NSType3NLN_Rhs = 0;
int *NSType3NLRhsSpace = nullptr;

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all blocks A11, A22
//      WITHOUT right hand sides
// ======================================================================

int NSType3NLSmagorinskyN_Terms = 3;
MultiIndex2D NSType3NLSmagorinskyDerivatives[3] = { D10, D01, D00 };
int NSType3NLSmagorinskySpaceNumbers[3] = { 0, 0, 0 };
int NSType3NLSmagorinskyN_Matrices = 4;
int NSType3NLSmagorinskyRowSpace[4] = { 0, 0, 0, 0 };
int NSType3NLSmagorinskyColumnSpace[4] = { 0, 0, 0, 0 };
int NSType3NLSmagorinskyN_Rhs = 0;
int *NSType3NLSmagorinskyRhsSpace = nullptr;

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A12, A21, A22
//      WITH right hand sides
//  for Newton method
// ======================================================================

int NSType3NLNewtonN_Matrices = 4;
int NSType3NLNewtonRowSpace[4] = { 0, 0, 0, 0 };
int NSType3NLNewtonColumnSpace[4] = { 0, 0, 0, 0 };
int NSType3NLNewtonN_Rhs = 2;
int NSType3NLNewtonRhsSpace[2] = {0, 0};

/**************************** NSTYPE == 4 ******************************/

/** assemble all matrices and rhs */

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int NSType4N_Terms = 4;
MultiIndex2D NSType4Derivatives[4] = { D10, D01, D00, D00 };
int NSType4SpaceNumbers[4] = { 0, 0, 0, 1 };
int NSType4N_Matrices = 8;
int NSType4RowSpace[8] = { 0, 0, 0, 0, 1, 1, 0, 0 };
int NSType4ColumnSpace[8] = { 0, 0, 0, 0, 0, 0, 1, 1 };
int NSType4N_Rhs = 2;
int NSType4RhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int NSType4SDN_Terms = 8;
MultiIndex2D NSType4SDDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType4SDSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType4SDN_Matrices = 8;
int NSType4SDRowSpace[8] = { 0, 0, 0, 0, 1, 1, 0, 0 };
int NSType4SDColumnSpace[8] = { 0, 0, 0, 0, 0, 0, 1, 1 };
int NSType4SDN_Rhs = 2;
int NSType4SDRhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int NSType4EquOrdN_Terms = 8;
MultiIndex2D NSType4EquOrdDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType4EquOrdSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType4EquOrdN_Matrices = 9;
int NSType4EquOrdRowSpace[9] = { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
int NSType4EquOrdColumnSpace[9] = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
int NSType4EquOrdN_Rhs = 3;
int NSType4EquOrdRhsSpace[3] = { 0, 0, 1 };

/** assemble some matrices and no rhs */

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A22
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int NSType4NLN_Terms = 3;
MultiIndex2D NSType4NLDerivatives[3] = { D10, D01, D00 };
int NSType4NLSpaceNumbers[3] = { 0, 0, 0 };
int NSType4NLN_Matrices = 2;
int NSType4NLRowSpace[2] = { 0, 0 };
int NSType4NLColumnSpace[2] = { 0, 0 };
int NSType4NLN_Rhs = 0;
int *NSType4NLRhsSpace = nullptr;

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType4NLSDN_Terms = 8;
MultiIndex2D NSType4NLSDDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType4NLSDSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType4NLSDN_Matrices = 4;
int NSType4NLSDRowSpace[4] = { 0, 0, 0, 0 };
int NSType4NLSDColumnSpace[4] = { 0, 0, 1, 1 };
int NSType4NLSDN_Rhs = 2;
int NSType4NLSDRhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType4NLEquOrdN_Terms = 8;
MultiIndex2D NSType4NLEquOrdDerivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
int NSType4NLEquOrdSpaceNumbers[8] = { 0, 0, 0, 0, 0, 1, 1, 1  };
int NSType4NLEquOrdN_Matrices = 8;
int NSType4NLEquOrdRowSpace[8] = { 0, 0, 0, 0, 0, 0, 1, 1 };
int NSType4NLEquOrdColumnSpace[8] = { 0, 0, 0, 0, 1, 1, 0, 0 };
int NSType4NLEquOrdN_Rhs = 2;
int NSType4NLEquOrdRhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all blocks Aij
//      WITH B1T, B2T (gradient blocks)
//      WITH right hand sides
// ======================================================================

int NSType4NLSmagorinskyN_Terms = 3;
MultiIndex2D NSType4NLSmagorinskyDerivatives[3] = { D10, D01, D00 };
int NSType4NLSmagorinskySpaceNumbers[3] = { 0, 0, 0 };
int NSType4NLSmagorinskyN_Matrices = 4;
int NSType4NLSmagorinskyRowSpace[4] = { 0, 0, 0, 0 };
int NSType4NLSmagorinskyColumnSpace[4] = { 0, 0, 0, 0 };
int NSType4NLSmagorinskyN_Rhs = 0;
int *NSType4NLSmagorinskyRhsSpace = nullptr;


// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A12, A21, A22
//      WITH right hand sides
//  for Newton method
// ======================================================================

int NSType4NLNewtonN_Matrices = 4;
int NSType4NLNewtonRowSpace[4] = { 0, 0, 0, 0 };
int NSType4NLNewtonColumnSpace[4] = { 0, 0, 0, 0 };
int NSType4NLNewtonN_Rhs = 2;
int NSType4NLNewtonRhsSpace[2] = {0, 0};
// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A12, A21, A22
//      WITH right hand sides
//  for Newton method and sdfem
// ======================================================================

int NSType4NLSDNewtonN_Matrices = 6;
int NSType4NLSDNewtonRowSpace[6] = { 0, 0, 0, 0 , 0, 0};
int NSType4NLSDNewtonColumnSpace[6] = { 0, 0, 0, 0, 1, 1};
/*
// ======================================================================
//  declarations for auxiliary problem for differential filter
//      one matrix 
//      two rhs
// ======================================================================

int Filter_N_Terms = 3;
MultiIndex2D Filter_Derivatives[3] = { D10, D01, D00};
int Filter_SpaceNumbers[3] = { 0, 0, 0};
int Filter_N_Matrices = 1;
int Filter_RowSpace[1] = { 0 };
int Filter_ColumnSpace[1] = { 0 };
int Filter_N_Rhs = 2;
int Filter_RhsSpace[2] = { 0, 0 };
*/
// ======================================================================
// declaration for pressure separation
//    only rhs
// ======================================================================


int NSPressSepN_Terms = 1;
MultiIndex2D NSPressSepDerivatives[1] = { D00 };
int NSPressSepSpaceNumbers[1] = { 0 };
int NSPressSepN_Matrices = 0;
int *NSPressSepRowSpace = nullptr;
int *NSPressSepColumnSpace = nullptr;
int NSPressSepN_Rhs = 2;
int NSPressSepRhsSpace[2] = { 0, 0 };

// ======================================================================
// declaration for pressure separation
// with auxiliary problem
// ======================================================================

int NSPressSepAuxProbN_Terms = 2;
MultiIndex2D NSPressSepAuxProbDerivatives[2] = { D10, D01 };
int NSPressSepAuxProbSpaceNumbers[2] = { 0, 0 };
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
MultiIndex2D NSRFBRhsDerivatives[1] = { D00 };
int NSRFBRhsSpaceNumbers[1] = { 0 };
int NSRFBRhsN_Matrices = 0;
int *NSRFBRhsRowSpace = nullptr;
int *NSRFBRhsColumnSpace = nullptr;
int NSRFBRhsN_Rhs = 2;
int NSRFBRhsRhsSpace[2] = { 0, 0 };
