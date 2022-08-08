// =======================================================================
// %W% %G%
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE3D__
#define __ASSEMBLE3D__

class TFESpace3D;
class TSquareMatrix3D;
class TMatrix3D;
#include <Constants.h>

#include "LocalAssembling.h"

/**
 * @brief Assembling method which takes a LocalAssembling3D object
 *
 * The method is adapted for use in MPI setting.
 *
 */
void Assemble3D(int n_fespaces, const TFESpace3D** fespaces,
				int n_sqmatrices, TSquareMatrix3D** sqmatrices,
				int n_matrices, TMatrix3D** matrices,
				int n_rhs, double** rhs,  const TFESpace3D** ferhs,
				BoundCondFunct3D** BoundaryConditions,
                BoundValueFunct3D** BoundaryValues,
				LocalAssembling3D& la, bool assemble_dirichlet_rows = false);


/** a function from a finite element space */
void Assemble3DSlipBC(int n_fespaces, const TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, const TFESpace3D **ferhs,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues);

void ModifyMatrixSlipBC(TSquareMatrix3D **sqmatrices, TMatrix3D **matrices,
			int N_U, double *rhs);

/** assemble mixed finite elements such as Raviart-Thomas or
 * Brezzi-Douglas-Marini.
 */
void Assemble3D_mixed(int n_fespaces, const TFESpace3D** fespaces,
                      int n_sqmatrices, TSquareMatrix3D** sqmatrices,
                      int n_matrices, TMatrix3D** matrices, int n_rhs,
                      double** rhs, const TFESpace3D** ferhs,
                      LocalAssembling3D& la,
                      BoundCondFunct3D** BoundaryConditions,
                      BoundValueFunct3D** BoundaryValues);



#endif // __ASSEMBLE3D__
