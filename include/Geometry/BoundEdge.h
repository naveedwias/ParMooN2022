#ifndef __BOUNDEDGE__
#define __BOUNDEDGE__

#include "BoundaryJoint.h"
#include "BoundComp2D.h"
#include <cmath>

/** edge on a boundary component */
class TBoundEdge : public BoundaryJoint
{
  protected:
    /** boundary component to which this edge belongs to */
    const TBoundComp2D *BoundComp;

    /** parameter of starting point */
    double T_0;
    /** parameter of end point */
    double T_1;

  public:
    // Constructors
    /** initialize the edge with the boundary component bdcomp and the
        paramter of starting and end point t\_0, t\_1 */
    TBoundEdge(const TBoundComp2D *bdcomp, double t_0, double t_1);

    // Methods
    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();
    
    virtual bool is_isoparametric() const
    { return BoundComp->GetType() != Line; }

    /** return start parameter T0 */
    double GetStartParameter() const
    { return T_0; }

    /** return end paramter T1 */
    double GetEndParameter() const
    { return T_1; }

    /** return parameters */
    void GetParameters(double& t0, double& t1) const
    {
      t0 = T_0;
      t1 = T_1;
    }
    
    void get_vertices(double &x0, double &y0, double &x1, double &y1) const
    {
        GetXYofT( T_0, x0, y0);
        GetXYofT( T_1, x1, y1);
    }

    double get_length() const
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      return std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    }

    void get_normal(double &nx, double &ny) const
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      double length =  std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
      nx =  (y1-y0)/length;
      ny = (x0-x1)/length; //(x1-x0)/length; //
    }

    void get_tangent(double &tx, double &ty) const
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      double length = std::sqrt( (x1-x0) * (x1-x0) + (y1-y0) * (y1-y0) );
      tx = (x1-x0)/length;
      ty = (y1-y0)/length;
    }

    /** @brief return the coordinates {X,Y} of parameter value T */
    void GetXYofT(double T, double &X, double &Y) const
    { BoundComp->GetXYofT(T, X, Y); }
    /** @brief return parameter value T of the coordinates {X,Y} */
    void GetTofXY(double X, double Y, double& T) const
    { BoundComp->GetTofXY(X, Y, T); }

    /** return boundary component */
    const TBoundComp2D *GetBoundComp() const
    { return BoundComp; }
};

#endif
