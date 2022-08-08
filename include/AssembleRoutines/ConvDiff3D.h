// ======================================================================
// %W% %G%
//
// common declaration for all 3D convection diffusion problems
// TODO CB Refactoring and restructuring of the assembling code.
// Hack: Lots of global variables put static to avoid linker errors when
// including this file more than once.
// ======================================================================

#ifndef __CONVDIFF3D__
#define __CONVDIFF3D__

// part for standard Galerkin - COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int N_Terms = 4;
//static MultiIndex3D Derivatives[4] = { D100, D010, D001, D000 };
//static int SpacesNumbers[4] = { 0, 0, 0, 0 };

// part for SDFEM (without 2nd derivatives) - COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int N_Terms_SD = 4;
//static MultiIndex3D Derivatives_SD[4] = { D100, D010, D001, D000};
//static int SpacesNumbers_SD[4] = { 0, 0, 0, 0 };

/*
// part for UPWIND with lumping of reaction term and rhs
int N_Terms_UPW1 = 2;
MultiIndex3D Derivatives_UPW1[2] = { D10, D01 };
int SpacesNumbers_UPW1[2] = { 0, 0 };
*/

// part for UPWIND without lumping of reaction term and rhs - COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int N_Terms_UPW2 = 4;
//static MultiIndex3D Derivatives_UPW2[4] = { D100, D010, D001, D000 };
//static int SpacesNumbers_UPW2[4] = { 0, 0, 0, 0 };

// part for all - COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int N_Matrices = 1;
//static int RowSpace[1] = { 0 };
//static int ColumnSpace[1] = { 0 };
//static int N_Rhs = 1;
//static int RhsSpace[1] = { 0 };
//
//static MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };


/*
void BilinearAssemble_UPW1(double Mult, double *coeff, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
*/

void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// parameters for DC/CD shock capturing scheme
void DC_CD_Params(double *in, double *out);

// COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int DC_CD_N_FESpaces = 1;
//static int DC_CD_N_Fct = 2;
//static int DC_CD_N_ParamFct = 1;
//static int DC_CD_N_FEValues = 2;
//static int DC_CD_N_Params = 2;
//static int DC_CD_FEFctIndex[2] = { 0, 1 };
//static MultiIndex3D DC_CD_FEMultiIndex[2] = { D000, D000 };
//static ParamFct *DC_CD_Fct[1] = { DC_CD_Params };
//static int DC_CD_BeginParam[1] = { 0 };

// parameters for MBE shock capturing scheme
void MBE_Params(double *in, double *out);

// COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int MBE_N_FESpaces = 1;
//static int MBE_N_Fct = 1;
//static int MBE_N_ParamFct = 1;
//static int MBE_N_FEValues = 4;
//static int MBE_N_Params = 4;
//static int MBE_FEFctIndex[4] = { 0, 0, 0, 0  };
//static MultiIndex3D MBE_FEMultiIndex[4] = { D000, D100, D010, D001 };
//static ParamFct *MBE_Fct[1] = { MBE_Params };
//static int MBE_BeginParam[1] = { 0 };

// parameters for SC_2 shock capturing scheme
void SC_2_Params(double *in, double *out);

// COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int SC_2_N_FESpaces = 2;
//static int SC_2_N_Fct = 3;
//static int SC_2_N_ParamFct = 1;
//static int SC_2_N_FEValues = 5;
//static int SC_2_N_Params = 5;
//static int SC_2_FEFctIndex[5] = { 0, 1, 2, 2, 2 };
//static MultiIndex3D SC_2_FEMultiIndex[5] = { D000, D000, D100, D010, D001 };
//static ParamFct *SC_2_Fct[1] = { SC_2_Params };
//static int SC_2_BeginParam[1] = { 0 };

// parameters for SOLD schemes
void SOLD_Params(double *in, double *out);

// COMMENTED IN ORDER TO AVOID IRRELEVANT COMPILER WARNINGS
//static int SOLD_N_FESpaces = 2;
//static int SOLD_N_Fct = 3;
//static int SOLD_N_ParamFct = 1;
//static int SOLD_N_FEValues = 6;
//static int SOLD_N_Params = 6;
//static int SOLD_FEFctIndex[6] = { 0, 0, 0, 0, 1, 2 };
//static MultiIndex3D SOLD_FEMultiIndex[6] = { D000, D100, D010, D001, D000, D000 };
//static ParamFct *SOLD_Fct[1] = { SOLD_Params };
//static int SOLD_BeginParam[1] = { 0 };


#endif // __CONVDIFF3D__
