#ifndef __BOUNDFACE__
#define __BOUNDFACE__

#include "BoundComp3D.h"
#include "BoundaryJoint.h"

/** face on a boundary component */
class TBoundFace : public BoundaryJoint
{
  protected:
    // boundary component to which this face belongs
    const TBoundComp3D *BoundComp;

    // parameter values for the vertices
    double Param1[4];
    double Param2[4];

  public:
    // Constructors
    TBoundFace(const TBoundComp3D *bdcomp, double *param1, double *param2);

    explicit TBoundFace(const TBoundComp3D *bdcomp);

    // Methods
    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();
    
    virtual bool is_isoparametric() const
    { return BoundComp->GetType() != Plane; }

    /** return boundary component */
    const TBoundComp3D *GetBoundComp() const
    { return BoundComp; }
    
    /// @brief return the coordinates of a point from its parametrization
    void GetXYZofTS(double T, double S, double &X, double &Y, double &Z) const
    { this->BoundComp->GetXYZofTS(T, S, X, Y, Z); }
    
    /// @brief return the parametrization of a point from its coordinates
    void GetTSofXYZ(double X, double Y, double Z, double &T, double &S) const
    { this->BoundComp->GetTSofXYZ(X, Y, Z, T, S); }
    
    void get_normal_vector(double x, double y, double z,
                           double& nx, double& ny, double &nz) const
    { this->BoundComp->get_normal_vector(x,y,z,nx,ny,nz);}
    
    /** return both parameters arrays */
    void GetParameters(double *param1, double *param2) const;

    void SetParameters(double *param1, double *param2);
};

#endif
