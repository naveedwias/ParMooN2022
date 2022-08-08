// =======================================================================
// @(#)BdPlane.h        1.1 07/16/99
//
// Class:       TBdPlane
// Superclass:  TBoundComp
// Purpose:     a plane as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#ifndef __BDPLANE__
#define __BDPLANE__

#include <BoundComp3D.h>

/** @brief a plane is described by:
    a point P=(P_x,P_y,P_z)
    a (tangential) vector A=(A_x,A_y,A_z)
    a (normal) vector N=(N_x,N_y,N_z)
    the vector B is computed as B = A x N

 **/
/** a plane as a component of a boundary part */
class TBdPlane : public TBoundComp3D
{
  protected:
      /** coordinates of one point in the plane */
      double P_x, P_y, P_z;
      /** orthogonal vectors in the plain (used for parametization) */
      double A_x, A_y, A_z, B_x, B_y, B_z;
      /** outer normal Vector */
      double N_x, N_y, N_z;

  public:
    // Constructor
    explicit TBdPlane(int id);
    TBdPlane(int id, int ref);
    
    virtual ~TBdPlane () {};

    // Methods
    /** set all parameters to the given values */
    void SetParams (double p_x, double p_y, double p_z,
                    double a_x, double a_y, double a_z,
                    double n_x, double n_y, double n_z);

    void GetParams (double &p_x, double &p_y, double &p_z,
                    double &a_x, double &a_y, double &a_z,
                    double &n_x, double &n_y, double &n_z) const;
    
    /** return the coordinates of parameter value T, S */
    virtual int GetXYZofTS(double T, double S,
                           double &X, double &Y, double &Z) const override;

    /** return the parameter value T, S of coordinates */
    virtual int GetTSofXYZ(double X, double Y, double Z,
                           double &T, double &S) const override;

/** 
    @brief return local parameters (t,s) and coordinates (x,y,z) 
    of a given linear combination of vertices 

    @param[in] N_points: the number of input points
    @param[in] LinComb: the coefficients of the linear combination
    @param[in] tp: the array of local t-coordinates
    @param[in] sp: the array of local s-coordinates
    @param[out] T = sum_i LinComb[i] tp[i]
    @param[out] S = sum_i LinComb[i] sp[i]
    @param[out] X = sum_i LinComb[i] xp[i]
    @param[out] Y = sum_i LinComb[i] yp[i]
    @param[out] Z = sum_i LinComb[i] zp[i]
**/

    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S) const override;

    /** read parameter from input file */
    virtual int ReadIn(std::istream &dat) override;

    virtual void get_normal_vector(double, double, double, double& nx,
                                   double& ny, double &nz) const override
    {
      nx = N_x;
      ny = N_y;
      nz = N_z;
    };
    
};

#endif
