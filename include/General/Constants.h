// =======================================================================
// @(#)Constants.h        1.15 04/13/00
// 
// Purpose:     constains important constants which are used by several 
//              other classes
//
// Author:      Volker Behns  09.07.97
//
// =======================================================================

#ifndef __CONSTANTS__
#define __CONSTANTS__

#include <functional>

constexpr int N_BOUNDCOND = 7;
enum BoundCond {DIRICHLET, NEUMANN, ROBIN, SLIP,
                SLIP_FRICTION_PENETRATION_RESISTANCE, DIRICHLET_WEAK, PERIODIC};

enum JointType {Joint, JointEqN, BoundaryPoint, BoundaryEdge, BoundaryFace,
                InterfaceJoint, PeriodicJoint, IsoInterfaceJoint,
                IsoJointEqN, IsoBoundEdge, IsoBoundFace,
                Joint_2to1,
                InterfaceJoint3D, IsoInterfaceJoint3D,
                SubDomainJoint, SubDomainHaloJoint, InnerInterfaceJoint,
                InnerEdge, IsoEdge3D, BDEdge3D};

typedef void DoubleFunct1D(double, double *);
typedef void DoubleFunct2D(double, double, double *);
typedef void DoubleFunct3D(double, double, double, double *);
typedef void DoubleFunctND(int, double *, double *);
typedef void DoubleFunctVect(const double *, double *);
typedef int IntFunct2D(double, double);
typedef double DoubleFunct2Param(double, double);

typedef void BoundCondFunct3D(int, double, double, double, BoundCond &);
typedef void BoundValueFunct3D(int, double, double, double, double &);
typedef void BoundCondFunct2D(int, double, BoundCond &);
typedef void BoundValueFunct2D(int, double, double &);

typedef DoubleFunctVect ParamFct;

typedef std::function<void(int, const double*, const double*,
                           const double*const*, double**)> CoeffFct2D;

typedef std::function<void(int, const double*, const double*, const double*,
                           const double*const*, double**)> CoeffFct3D;

typedef void AssembleFct2D(double, const double *, double, const double **, 
                           const int *, double ***, double **);

typedef std::function<void(double, const double *, const double *, double,
                           const double **, const int *, double ***, double **)>
   AssembleFctParam;

typedef void AssembleFct3D(double, const double *, double, const double **, 
                           const int *, double ***, double **);

class TBaseCell;
typedef void ManipulateFct(int, double **, const double **, const TBaseCell *);

class TCollection;
typedef void EvalAllNF(const TCollection *, const TBaseCell *, const double *, double *);
typedef void EvalJointNF(const TCollection *, const TBaseCell *, int, const double *, double *);

class TFESpace2D;
typedef void CheckWrongNeumannNodesFunct2D(TCollection *, TFESpace2D *,
					   int &, int* &,
					   int* &, 
					   double* &);

class TSquareMatrix;
class TMatrix;
typedef void MatVecProc(TSquareMatrix **, TMatrix **, double *, double *);
typedef void DefectProc(TSquareMatrix **, TMatrix **, double *, double *,
                        double *);

class TVertex;
class TIsoBoundEdge;
typedef void ModifyMeshCoords(double , double , double &, double &, double );
typedef void ModifyBoundCoords(int, TVertex **, TIsoBoundEdge **,  double *, double);

#ifdef _MPI
class TParVectorNSE3D;
typedef void ParDefectProc(TSquareMatrix **, TMatrix **, TParVectorNSE3D  *, TParVectorNSE3D *,
                        TParVectorNSE3D *);
#endif

typedef int TypeBoundSwitchFunct2D(int, double );

class TFEFunction2D;
class TFEVectFunct2D;
typedef void EvaluateSolutionFunct2D(TFEFunction2D **, TFEVectFunct2D **, 
                                     double *, int *);

// https://stackoverflow.com/a/4609795
template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

// maybe we should define more types here (like CD, TCD, TNSE, TSTOKES,...)
#define NSE    0
#define STOKES 1
#define OSEEN  2

#endif
