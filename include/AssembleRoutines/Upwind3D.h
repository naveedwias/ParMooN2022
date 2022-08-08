// =======================================================================
// @(#)Upwind.h        1.4 10/18/99
//
// Purpose:     upwind stabilization
//              do upwind assembling for first order nonconforming elements
//
// Authors:     Volker Behns / Gunar Matthies  18.10.99
//
// =======================================================================


#ifndef __UPWIND3D__
#define __UPWIND3D__

#include <SquareMatrix3D.h>
#include <FEFunction3D.h>

template <int d> class LocalAssembling;

/// Do upwinding for Navier--Stokes. TODO Understand and comment this.
/// The method was freed form its global database dependency by setting the
/// four needed values as input values, three of them with defaults (same as
/// in old DB). Trouble is: we do not yet know what these parameters do.
void UpwindForNavierStokes3D(TSquareMatrix3D *sqmatrix, const TFEFunction3D *u1,
                             const TFEFunction3D *u2, const TFEFunction3D *u3,
                             double one_over_nu, // that is former TDatabase::ParamDB->RE_NR,
                                                  // but I found naming this "Reynolds number"
                                                  // somewhat misleading (and the database dependency not good)
                             int upwind_order = 1,
                             double upwind_flux_damp = 1,
                             int upwind_application = 0);

void UpwindForConvDiff(TSquareMatrix3D *sqmatrix, double *RHS,
                       const TFESpace3D *fespace, const LocalAssembling<3>& la);

#endif
