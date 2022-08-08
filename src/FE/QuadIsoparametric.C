#include "QuadIsoparametric.h"
#include "BaseCell.h"
#include "BoundComp.h"
#include "BoundEdge.h"
#include "FEDescriptor.h"
#include "InterfaceJoint.h"
#include "IsoBoundEdge.h"
#include "IsoInterfaceJoint.h"
#include "IsoJointEqN.h"
#include "LinAlg.h"
#include "QuadFormula.h"

BaseFunction_type BaseFunctFromOrder[] = { 
                BF_C_Q_Q0_2D, BF_C_Q_Q1_2D, BF_C_Q_Q2_2D, BF_C_Q_Q3_2D,
                BF_C_Q_Q4_2D, BF_C_Q_Q5_2D, BF_C_Q_Q6_2D, BF_C_Q_Q7_2D,
                BF_C_Q_Q8_2D, BF_C_Q_Q9_2D };

FEDescriptor_type FEDescFromOrder[] = { 
                FE_C_Q_Q0_2D, FE_C_Q_Q1_2D, FE_C_Q_Q2_2D, FE_C_Q_Q3_2D,
                FE_C_Q_Q4_2D, FE_C_Q_Q5_2D, FE_C_Q_Q6_2D, FE_C_Q_Q7_2D,
                FE_C_Q_Q8_2D, FE_C_Q_Q9_2D };

/** constuctor */
TQuadIsoparametric::TQuadIsoparametric()
{
}

/** transfer from reference element to original element */
void TQuadIsoparametric::GetOrigFromRef(double xi, double eta,
                                        double &X, double &Y) const
{
  GetOrigFromRef(1, &xi, &eta, &X, &Y);
}

/** transfer a set of points from reference to original element */
void TQuadIsoparametric::GetOrigFromRef(int N_Points, const double *xi,
                                        const double *eta, double *X, double *Y)
  const
{
  int i, j, k;
  double Xi, Eta;
  double dx1, dx2;
  double dy1, dy2;
  double AuxVector[3*MaxN_BaseFunctions2D];
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];

    dx1 = xc1 + xc3*Eta;
    dx2 = xc2 + xc3*Xi;
    
    dy1 = yc1 + yc3*Eta;
    dy2 = yc2 + yc3*Xi;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Xi*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Xi*Eta;

    bf.GetDerivatives(MultiIndex2D::D00, Xi, Eta, AuxVector);
    bf.GetDerivatives(MultiIndex2D::D10, Xi, Eta, AuxVector+MaxN_BaseFunctions2D);
    bf.GetDerivatives(MultiIndex2D::D01, Xi, Eta, AuxVector+2*MaxN_BaseFunctions2D);

    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
    }
  }
}

void TQuadIsoparametric::GetOrigFromRef(const TQuadFormula& qf_ref,
                                        TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  double AuxVector[3*MaxN_BaseFunctions2D];
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  for(unsigned int i = 0; i < n_points; ++i)
  {
    double Xi = qf_ref.get_point(i).x;
    double Eta = qf_ref.get_point(i).y;
    
    double dx1 = xc1 + xc3*Eta;
    double dx2 = xc2 + xc3*Xi;
    
    double dy1 = yc1 + yc3*Eta;
    double dy2 = yc2 + yc3*Xi;
    
    double x = xc0 + xc1*Xi + xc2*Eta + xc3*Xi*Eta;
    double y = yc0 + yc1*Xi + yc2*Eta + yc3*Xi*Eta;
    
    bf.GetDerivatives(MultiIndex2D::D00, Xi, Eta, AuxVector);
    bf.GetDerivatives(MultiIndex2D::D10, Xi, Eta, AuxVector+MaxN_BaseFunctions2D);
    bf.GetDerivatives(MultiIndex2D::D01, Xi, Eta, AuxVector+2*MaxN_BaseFunctions2D);
    
    for(int k=0;k<N_AuxPoints;k++)
    {
      int j = IntAux[k];
      x += XDistance[k] * AuxVector[j];
      y += YDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
    }
    double absdet = std::abs(dx1*dy2 - dx2*dy1);
    double weight = qf_ref.get_weight(i) * absdet;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y}});
  }
}

/** transfer from original element to reference element */
void TQuadIsoparametric::GetRefFromOrig(double, double, double &xi, double &eta)
  const
{
  cout << "not implemented yet" << endl;
  xi  = 0.0;
  eta = 0.0;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadIsoparametric::GetOrigValues(double xi, double eta, int N_BaseFunct,
                                       const double *uref,
                                       const double *uxiref,
                                       const double *uetaref,
                                       double *uorig, double *uxorig,
                                       double *uyorig, int) const
{
  int i, j, k;
  double dx1, dx2;
  double dy1, dy2;
  double rec_detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];

  dx1 = xc1 + xc3*eta;
  dx2 = xc2 + xc3*xi;
    
  dy1 = yc1 + yc3*eta;
  dy2 = yc2 + yc3*xi;
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  bf.GetDerivatives(MultiIndex2D::D10, xi, eta, AuxVector+MaxN_BaseFunctions2D);
  bf.GetDerivatives(MultiIndex2D::D01, xi, eta, AuxVector+2*MaxN_BaseFunctions2D);

  for(k=0;k<N_AuxPoints;k++)
  {
    j = IntAux[k];

    dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

    dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
  }

  rec_detjk = 1/(dx1*dy2 - dx2*dy1);

  // D00
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( dy2*uxiref[i]-dy1*uetaref[i]) * rec_detjk;
    uyorig[i]=(-dx2*uxiref[i]+dx1*uetaref[i]) * rec_detjk;
  } // endfor i
}

void TQuadIsoparametric::GetOrigAllDerivatives(
  double xi, double eta, int N_BaseFunct,
  const double *r00, const double *r10, const double *r01,
  const double *r20, const double *r11, const double *r02,
  double *o00, double *o10, double *o01, double *o20, double *o11, double *o02,
  int BaseVectDim) const
{
  GetOrigValues(xi, eta, N_BaseFunct, r00, r10, r01, o00, o10, o01,
                BaseVectDim);
  bool all_zero_2nd_derivatives = true;
  for(int i=0;i<N_BaseFunct;i++)
  {
    if(r20[i] != 0. || r11[i] != 0. || r02[i] != 0.)
    {
      all_zero_2nd_derivatives = false;
      break;
    }
  } // endfor i
  if(!all_zero_2nd_derivatives)
    Output::warn("QuadIsoparametric", "transformation of second order "
                 "derivatives is not yet implementd.");
  // avoid compiler warnings
  (void)o20;
  (void)o11;
  (void)o02;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadIsoparametric::GetOrigValues(int joint, double zeta, int N_BaseFunct,
                                       const double *uref,const double *uxiref,
                                       const double *uetaref,
                                       double *uorig, double *uxorig,
                                       double *uyorig, int) const
{
  int i, k;
  double a11, a12, a21, a22, rec_detjk;
  double xi=0, eta=0;
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  double valxi[MaxN_BaseFunctions2D];
  double valeta[MaxN_BaseFunctions2D];
  double *values[1];

  switch(joint)
  {
    case 0:
      xi  = zeta; eta = -1;
    break;

    case 1:
      xi = 1; eta = zeta;
    break;

    case 2:
      xi  = -zeta; eta = 1;
    break;

    case 3:
      xi = -1 ; eta = -zeta;
    break;
  }

  a11 = xc1 + xc3*eta;
  a21 = xc2 + xc3*xi;
  a12 = yc1 + yc3*eta;
  a22 = yc2 + yc3*xi;

  values[0] = valxi;
  bf.GetDerivatives(MultiIndex2D::D10, 1, &zeta, joint, values);
  values[0] = valeta;
  bf.GetDerivatives(MultiIndex2D::D01, 1, &zeta, joint, values);

  // add correction due to isoparamatric boundary approximation
  for(i=0;i<N_AuxPoints;i++)
  {
    k = IntAux[i];
    a11 += XDistance[i] * valxi[k];
    a21 += XDistance[i] * valeta[k];

    a12 += YDistance[i] * valxi[k];
    a22 += YDistance[i] * valeta[k];
  } // endfor i

  rec_detjk = 1 / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
    uorig[i] = uref[i];
  } // endfor i
}

void TQuadIsoparametric::SetCell(const TBaseCell *cell)
{
  int i, j;
  const TBoundComp2D *comp;
  JointType type;
  double t0, t1, t, dt;
  double xa, ya, xe, ye, xm, ym, xp, yp, dx, dy;
//  int compid;
//  BoundTypes bdtype;
  int *JointDOF;
  int N_Vertices;
  TVertex **Vertices;
#ifdef __3D__
  double z[4], zp;
#endif
  int CurvedJoint;

  N_AuxPoints = 0;
  CurvedJoint = -1;

  Cell = cell;

  // cout << endl;
  for(i=0;i<4;i++)
  {
#ifdef __3D__
    Cell->GetVertex(i)->GetCoords(x[i], y[i], z[i]);
#else
    Cell->GetVertex(i)->GetCoords(x[i], y[i]);
#endif
    // cout << setw(20) << x[i] << setw(20) << y[i] << endl;
  }

  xc0=( x[0] + x[1] + x[2] + x[3]) * 0.25;
  xc1=(-x[0] + x[1] + x[2] - x[3]) * 0.25;
  xc2=(-x[0] - x[1] + x[2] + x[3]) * 0.25;
  xc3=( x[0] - x[1] + x[2] - x[3]) * 0.25;

  yc0=( y[0] + y[1] + y[2] + y[3]) * 0.25;
  yc1=(-y[0] + y[1] + y[2] - y[3]) * 0.25;
  yc2=(-y[0] - y[1] + y[2] + y[3]) * 0.25;
  yc3=( y[0] - y[1] + y[2] - y[3]) * 0.25;

  for(i=0;i<4;i++)
  {
    // check whether the joint i is curved
    auto joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryEdge || type == InterfaceJoint)
    {
      // cout << "joint: " << i << endl;
      xa = x[i];
      ya = y[i];
      xe = x[(i+1) % 4];
      ye = y[(i+1) % 4];

      if(type == BoundaryEdge)
      {
        auto boundedge = (const TBoundEdge *)(joint);
        comp = boundedge->GetBoundComp();
        boundedge->GetParameters(t0, t1);
      }
      else
      {
        auto interface = (const TInterfaceJoint *)(joint);
        comp = interface->GetBoundComp();
        if(Cell == interface->GetNeighbour(0))
          interface->GetParameters(t0, t1); // correct order
        else
          interface->GetParameters(t1, t0); // backward order

      }
      FEDescriptor fedesc(FEDescFromOrder[ApproximationOrder]);
      JointDOF = fedesc.GetJointDOF(i);

      dt = (t1-t0)/ApproximationOrder;
      dx = (xe-xa)/ApproximationOrder;
      dy = (ye-ya)/ApproximationOrder;
      for(j=1;j<ApproximationOrder;j++)
      {
        t = t0 + j * dt;
        comp->GetXYofT(t, xp, yp);
        xm = xa + j * dx;
        ym = ya + j * dy;
 
        // cout << "(" << xm << ", " << ym << ")" << endl;
        // cout << "(" << xp << ", " << yp << ")" << endl;

        if(std::abs(xp-xm) > 1e-8 || std::abs(yp-ym) > 1e-8)
        {
          // curved boundary
          XDistance[N_AuxPoints] = xp - xm;
          YDistance[N_AuxPoints] = yp - ym;
          IntAux[N_AuxPoints] = JointDOF[j];
          // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
          if(N_AuxPoints == 0)
            CurvedJoint = i;

          if(CurvedJoint != i)
          {
            ErrThrow("There is only one curved joint allowed");
          }

          N_AuxPoints++;
        } // endif 
      } // endfor
    } // endif
    else
    {
      if(type == IsoInterfaceJoint || type == IsoBoundEdge ||
         type == IsoJointEqN)
      {
        switch(type)
        {
          case IsoInterfaceJoint:
            N_Vertices = ((TIsoInterfaceJoint *)joint)->GetN_Vertices();
            Vertices = ((TIsoInterfaceJoint *)joint)->GetVertices();
          break;

          case IsoBoundEdge:
            N_Vertices = ((TIsoBoundEdge *)joint)->GetN_Vertices();
            Vertices = ((TIsoBoundEdge *)joint)->GetVertices();
          break;

          case IsoJointEqN:
            N_Vertices = ((TIsoJointEqN *)joint)->GetN_Vertices();
            Vertices = ((TIsoJointEqN *)joint)->GetVertices();
          break;
	  default:
	    cout<< "Not a 2D BD type "<<endl;
	   break;
        } // endswitch

        // cout << "N_AuxVertices: " << N_Vertices << endl;
        ApproximationOrder = N_Vertices+1;
        FEDescriptor fedesc(FEDescFromOrder[ApproximationOrder]);
        JointDOF = fedesc.GetJointDOF(i);

        xa = x[i];
        ya = y[i];
        xe = x[(i+1) % 4];
        ye = y[(i+1) % 4];

        dt = (t1-t0)/ApproximationOrder;
        dx = (xe-xa)/ApproximationOrder;
        dy = (ye-ya)/ApproximationOrder;
        for(j=1;j<ApproximationOrder;j++)
        {
          t = t0 + j * dt;
#ifdef __3D__
          Vertices[j-1]->GetCoords(xp, yp, zp);
#else
          Vertices[j-1]->GetCoords(xp, yp);
#endif
          xm = xa + j * dx;
          ym = ya + j * dy;
  
          // cout << "(" << xm << ", " << ym << ")" << endl;
          // cout << "(" << xp << ", " << yp << ")" << endl;
  
          if(std::abs(xp-xm) > 1e-8 || std::abs(yp-ym) > 1e-8)
          {
            // curved boundary
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            IntAux[N_AuxPoints] = JointDOF[j];
            // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
            if(N_AuxPoints == 0)
              CurvedJoint = i;
  
            if(CurvedJoint != i)
            {
              ErrThrow("Only one curved joint per quadrilateral is allowed");
            }
  
            N_AuxPoints++;
          } // endif 
        } // endfor
      } // endif
    } // endelse
  } // endfor

  if(N_AuxPoints)
  {
    // Output::print("ApproximationOrder: ", ApproximationOrder);
    if(ApproximationOrder == 2)
    {
      XDistance[N_AuxPoints] = XDistance[0]*0.5;
      YDistance[N_AuxPoints] = YDistance[0]*0.5;
      IntAux[N_AuxPoints] = 4;
      N_AuxPoints++;
    } // ApproximationOrder == 2

    if(ApproximationOrder == 3)
    {
      XDistance[N_AuxPoints] = XDistance[0]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/3;
      N_AuxPoints++;

      N_AuxPoints -= 4;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
        break;
      } // endswitch
    } // ApproximationOrder == 3

    if(ApproximationOrder == 4)
    {
      XDistance[N_AuxPoints] = XDistance[0]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/4;
      N_AuxPoints++;

      N_AuxPoints -= 9;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
        break;
      } // endswitch
    } // ApproximationOrder == 4

    if(ApproximationOrder == 5)
    {
      XDistance[N_AuxPoints] = XDistance[0]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/5;
      N_AuxPoints++;

      N_AuxPoints -= 16;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
        break;
      } // switch
    } // ApproximationOrder == 5
  }

  if(N_AuxPoints)
  {
    BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
    unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
    for(unsigned int i=0;i<N_QuadPoints;i++)
    {
      auto p = quadrature_formula->get_point(i);
      bf.GetDerivatives(MultiIndex2D::D00, p.x, p.y, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]]; 

      bf.GetDerivatives(MultiIndex2D::D10, p.x, p.y, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf.GetDerivatives(MultiIndex2D::D01, p.x, p.y, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]]; 

    } // endfor i
  } // endif
} // TQuadIsoparametric::SetCell

void TQuadIsoparametric::SetQuadFormula(QuadratureFormula_type formula)
{
  quadrature_formula = std::unique_ptr<TQuadFormula>(new TQuadFormula(formula));
}

void TQuadIsoparametric::PiolaMapOrigFromRef(double, double, int,
                                             const double *, double *) const
{
  ErrThrow("TQuadIsoparametric::PiolaMapOrigFromRef not implemented");
}
    
void TQuadIsoparametric::PiolaMapOrigFromRef(double, double, int,
                                             const double *, const double *,
                                             const double *, double *, double *)
  const
{
  ErrThrow("TQuadIsoparametric::PiolaMapOrigFromRef not implemented");
}
