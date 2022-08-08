// =======================================================================
// @(#)Assemble2D.h        1.5 04/13/00
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE2D__
#define __ASSEMBLE2D__

class TAuxParam2D;
#include <Constants.h>

#include "LocalAssembling.h"

#include "Joint.h"

/** a function from a finite element space */
void Assemble2D(int n_fespaces,
                const TFESpace2D **fespaces,
                int n_sqmatrices,
                TSquareMatrix2D **sqmatrices,
                int n_matrices,
                TMatrix2D **matrices,
                int n_rhs,
                double **rhs,
                const TFESpace2D **ferhs,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                LocalAssembling2D& la,
                bool assemble_dirichlet_rows = false);

/*
std::vector<parmoon::Point> Transform_Quad_Points( const TBaseCell* cell, const
    BFRefElements& ref_element, const ReferenceTransformation_type&
    ref_trans_id, const int& joint_index, const TQuadFormula* quad_form );*/
/** a function from a finite element space */
void Assemble2D_JumpStab(int n_fespaces,
                const TFESpace2D **fespaces,
                int n_sqmatrices,
                TSquareMatrix2D **sqmatrices,
                int n_matrices,
                TMatrix2D **matrices,
                int n_rhs,
                double **rhs,
                const TFESpace2D **ferhs,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                LocalAssembling2D& la,
                bool assemble_dirichlet_rows = false);

/** assembling of slip type bc */
void Assemble2DSlipBC(int n_fespaces, const TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, const TFESpace2D **ferhs,
                      BoundCondFunct2D **BoundaryConditions,
                      BoundValueFunct2D **BoundaryValues,
                      // TAuxParam2D *parameters,
                      TFEFunction2D *u1, TFEFunction2D *u2);


/** assembling for continuous interior penalty discretization */
void Assemble2D_CIP(const CoeffFct2D& Coeff,int n_fespaces,
                    const TFESpace2D **fespaces,
                    int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                    int n_matrices, TMatrix2D **matrices,
                    int n_rhs, double **rhs, TFESpace2D **ferhs,
                    BoundCondFunct2D **BoundaryConditions,
                    BoundValueFunct2D **BoundaryValues,
                    TAuxParam2D *Parameters);

/** assembling for vector finite elements (Raviart-Thomas (RT) and
 * Brezzi-Douglas-Marini (BDM)) */
#ifdef __2D__
void Assemble2D_VectFE(int n_fespaces, const TFESpace2D **fespaces,
                       int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                       int n_matrices, TMatrix2D **matrices,
                       int n_rhs, double **rhs, const TFESpace2D **ferhs,
                       LocalAssembling2D& la,
                       BoundCondFunct2D **BoundaryConditions,
                       BoundValueFunct2D * const * const BoundaryValues
                       );
#endif // 2D

#endif // __ASSEMBLE2D__

