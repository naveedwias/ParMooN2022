#ifndef __NODALFUNCTIONAL2D__
#define __NODALFUNCTIONAL2D__

#include "Enumerations_fe.h"
#include "Constants.h"

/// @brief maximum number of evaluation points needed to evaluate a 2D nodal
/// functional
constexpr int MaxN_PointsForNodal2D = 100;
/// @brief maximum number of evaluation points needed to evaluate a 3D nodal
/// functional
constexpr int MaxN_PointsForNodal3D = 275;

/** @brief realize nodal functionals in 2D */
class NodalFunctional
{
  protected:
    /** @brief ID for this set of nodal functionals */
    NodalFunctional_type type;

    /** @brief number of all functionals */
    int N_AllFunctionals;

    /** @brief number of functionals on one edge */
    int N_EdgeFunctionals;

    /** @brief number of points needed for all nodal functionals */
    int N_PointsAll;

    /** @brief values at points with following xi-coordinates are needed
        to evaluate all nodal functionals */
    double *Xi;

    /** @brief values at points with following eta-coordinates are needed
        to evaluate all nodal functionals */
    double *Eta;

    /** @brief values at points with following zeta-coordinates are needed
        to evaluate all nodal functionals */
    double *Zeta;

    /** @brief routine for evaluating all functionals */
    EvalAllNF *EvalAll;

    /** @brief Number of points needed for edge nodal functionals */
    int N_PointsEdge;

    /** @brief values at edge points with following edge parameters in [-1,1]
        are needed to evaluate functional on edge. In 3D this is the 
        t-parameter for points for evaluating functional on face */
    double *T;

    /** @brief routine for evaluating the edge functionals */
    EvalJointNF *EvalJoint;

    /** @brief array of number of functionals on one face */
    int *N_FaceFunctionals;
    
    /** @brief array of numbers of points needed for face nodal functionals */
    int *N_PointsFace;

    /** @brief xi-coordinate of points for evaluating functional on face */
    double **XiArray;

    /** @brief eta-coordinate of points for evaluating functional on face */
    double **EtaArray;

    /** @brief zeta-coordinate of points for evaluating functional on face */
    double **ZetaArray;

    /** @brief s-parameter for points for evaluating functional on face */
    double *S;

  public:
    
    explicit NodalFunctional ( NodalFunctional_type id);

    /** @brief return information for points for all functionals */
    void GetPointsForAll(int &n_points, const double* &xi, const double* &eta) const
    { n_points = N_PointsAll; xi = Xi; eta = Eta; }
    
    void GetPointsForAll(int &n_points, const double* &xi, const double* &eta,
                         const double* &zeta) const;

    /** @brief return information for points for edge functionals */
    void GetPointsForEdge(int &n_points, const double* &t) const
    { n_points = N_PointsEdge; t = T; }

    /** @brief return information for points for face functionals 
        on joint j */
    void GetPointsForFace(int j, int &n_points, const double* &xi,
                          const double* &eta, const double* &zeta) const;

    /** @brief return information for points for face functionals */
    void GetPointsForFace(int &n_points, const double* &t, const double* &s)
      const;

    /** @brief return values for all nodal functionals */
    void GetAllFunctionals(const TCollection *Coll, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
      const
    { EvalAll(Coll, Cell, PointValues, Functionals); }

    /** @brief return values for edge nodal functional */
    void GetEdgeFunctionals(const TCollection *Coll, const TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
      const
    { EvalJoint (Coll, Cell, Joint, PointValues, Functionals); }
    
    /** @brief return values for face nodal functional */
    void GetFaceFunctionals(const TCollection *Coll, TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
      const
    { EvalJoint (Coll, Cell, Joint, PointValues, Functionals); }

    /** @brief return ID for this set */
    NodalFunctional_type GetID() const
    { return type; }
    
    int n_functionals() const
    {return N_AllFunctionals; }

    int n_face_functionals(int face) const
    {return N_FaceFunctionals[face]; }

    int n_points() const
    { return N_PointsAll; }
};

#endif
