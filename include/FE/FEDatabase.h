#ifndef INCLUDE_FE_FEDATABASE_H_
#define INCLUDE_FE_FEDATABASE_H_
#include "RefTrans2D.h"
#include "RefTrans3D.h"
#include "RefDesc.h"
#include "BaseFunctions.h"
#include "QuadFormula.h"
class BaseFunctions;
class FiniteElement;
class TGridCell;

namespace FEDatabase
{
  /** 
   * @brief free the memory of all fe databases in this namespace.
   * 
   * @note you also have to call this in all parmoon programs.
   */
  void destroy();
  //======================================================================
  //      prolongation matrices in 2D
  //======================================================================
  /** return prolongation matrix for given situation */
  double *GetProlongationMatrix2D(FE_type parent, Refinements refine,
                                  FE_type child, int childnumber);
  double *GetProlongationMatrix3D(FE_type parent, Refinements refine,
                                  FE_type child, int childnumber);

  //======================================================================
  //      function restriction matrices in 2D
  //======================================================================
  /** return function restriction matrix for given situation */
  double *GetRestrictionMatrix2D(FE_type parent, Refinements refine,
                                 FE_type child, int childnumber);
  double *GetRestrictionMatrix3D(FE_type parent, Refinements refine,
                                 FE_type child, int childnumber);

  //======================================================================
  //      FE function values and derivatives 2D
  //======================================================================
  /** return requested FE function values or derivatives */
  double **GetRefElementValues(const BaseFunctions& bf, const TQuadFormula& qf,
                               MultiIndex2D MultiIndex);
  double **GetRefElementValues(const BaseFunctions& bf, const TQuadFormula& qf,
                               MultiIndex3D MultiIndex);

  /** return requested joint values of FE base functions */
  double **GetJointDerivatives2D(const BaseFunctions& bf,
                                 const TQuadFormula& qf, int joint,
                                 MultiIndex2D MultiIndex);
  double **GetJointDerivatives3D(const BaseFunctions& bf,
                                 const TQuadFormula& qf, int joint,
                                 MultiIndex3D MultiIndex);

  /** return requested FE function values or derivatives */
  double **GetOrigElementValues(const BaseFunctions& bf,
                                MultiIndex2D MultiIndex);
  double **GetOrigElementValues(const BaseFunctions& bf,
                                MultiIndex3D MultiIndex);

  /** calculate functions and derivatives from reference element
   * to original element
   */
  void GetOrigValues(ReferenceTransformation_type RefTrans,
                     double xi, double eta, const BaseFunctions *bf,
                     const TCollection *Coll, const TGridCell *cell,
                     double *uref, double *uxiref, double *uetaref,
                     double *uorig, double *uxorig, double *uyorig);
  void GetOrigValues(ReferenceTransformation_type RefTrans,
                     double xi, double eta, double zeta,
                     const BaseFunctions *bf, const TCollection *Coll,
                     const TBaseCell *cell,
                     double *uref, double *uxiref, double *uetaref,
                     double *uzetaref,
                     double *uorig, double *uxorig, double *uyorig,
                     double *uzorig);
    
  void GetOrigValues(ReferenceTransformation_type RefTrans, double zeta,
                     const BaseFunctions *bf, int edgeNumber,
                     const TCollection *Coll, const TBaseCell *cell,
                     double *uref, double *uxiref, double *uetaref,
                     double *uorig, double *uxorig, double *uyorig);

  void GetOrigValues(ReferenceTransformation_type RefTrans,
                     const std::vector<const BaseFunctions*>& BaseFuncts,
                     const TCollection *Coll, const TBaseCell *cell,
                     const TQuadFormula& qf, bool *Needs2ndDer);
  //======================================================================
  //      reference transformation 2D
  //======================================================================
  /** return reference transformation */
  TRefTrans2D *GetRefTrans2D(ReferenceTransformation_type reftrans);
  TRefTrans3D *GetRefTrans3D(ReferenceTransformation_type reftrans);

  /** calculate points on original element */
  void GetOrigFromRef(ReferenceTransformation_type RefTrans, int n_points, 
                      const double *xi, const double *eta,
                      double *X, double *Y);
  void GetOrigFromRef(ReferenceTransformation_type RefTrans, int n_points,
                      const double *xi, const double *eta, const double *zeta,
                      double *X, double *Y, double *Z);
  void GetOrigFromRef(ReferenceTransformation_type RefTrans,
                      const TQuadFormula& qf_ref, TQuadFormula& qf_orig);

  /** calculate base functions with derivatives and coordinates
   * from reference to original element
   */
  ReferenceTransformation_type GetOrig(
    std::vector<const FiniteElement*> used_fe, const TCollection *Coll,
    const TBaseCell *cell, bool *Needs2ndDer,
    TQuadFormula& qf_ref, TQuadFormula& qf_orig);
     

  /** calculate points on reference element */
  void GetRefFromOrig(ReferenceTransformation_type RefTrans, double X, double Y,
                      double &xi, double &eta);
  void GetRefFromOrig(ReferenceTransformation_type RefTrans,
                      double X, double Y, double Z,
                      double &xi, double &eta, double &zeta);

  /** set cell for reference transformation */
  void SetCellForRefTrans(const TBaseCell *cell,
                          ReferenceTransformation_type reftrans);
}

#endif // INCLUDE_FE_FEDATABASE_H_
