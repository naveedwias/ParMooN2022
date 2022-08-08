#include "TriaIsoparametric.h"
#include "BaseCell.h"
#include "FEDescriptor.h"
#include "LinAlg.h"
#include "IsoBoundEdge.h"
#include "IsoInterfaceJoint.h"
#include "IsoJointEqN.h"
#include "BoundComp.h"
#include "QuadFormula.h"

BaseFunction_type TTriaIsoparametric::BaseFunctFromOrder[] = { 
                BF_C_T_P0_2D, BF_C_T_P1_2D, BF_C_T_P2_2D, BF_C_T_P3_2D,
                BF_C_T_P4_2D, BF_C_T_P5_2D, BF_C_T_P6_2D, BF_C_T_P7_2D,
                BF_C_T_P8_2D, BF_C_T_P9_2D };

FEDescriptor_type TTriaIsoparametric::FEDescFromOrder[] = { 
                FE_C_T_P0_2D, FE_C_T_P1_2D, FE_C_T_P2_2D, FE_C_T_P3_2D,
                FE_C_T_P4_2D, FE_C_T_P5_2D, FE_C_T_P6_2D, FE_C_T_P7_2D,
                FE_C_T_P8_2D, FE_C_T_P9_2D };

/** constuctor */
TTriaIsoparametric::TTriaIsoparametric()
{
}

/** transfer from reference element to original element */
void TTriaIsoparametric::GetOrigFromRef(double xi, double eta,
                                        double &X, double &Y) const
{
/*
  int i, j;

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(std::abs(xi-XI[i])<1e-8 && std::abs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TTriaIsoparametric::GetOrigFromRef(1)" << endl;
    return;
  }

  X = xc0 + xc1*xi + xc2*eta;
  Y = yc0 + yc1*xi + yc2*eta;

  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
  }

  // cout << "(" << xi << ", " << eta << ") = ";
  // cout << "(" << X << "," << Y << ")" << endl;
*/
  int N_Points;
  double Xi[1], Eta[1], Xarray[1], Yarray[1];

  N_Points = 1;
  Xi[0] = xi; Eta[0] = eta;
  Xarray[0] = X; Yarray[0] = Y;
  GetOrigFromRef(N_Points, Xi, Eta, Xarray, Yarray);
  X = Xarray[0]; Y = Yarray[0];
}

/** transfer a set of points from reference to original element */
void TTriaIsoparametric::GetOrigFromRef(int N_Points, const double *xi,
                                        const double *eta, double *X, double *Y)
  const
{
/*
  int i, j, k;
  double Xi, Eta, a11, a12, a21, a22;

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];

    a11 = xc1;
    a21 = xc2;
    a12 = yc1;
    a22 = yc2;

    j = -1;
    for(k=0;k<N_QuadPoints;k++)
    {
      if(std::abs(Xi-XI[k])<1e-8 && std::abs(Eta-ETA[k])<1e-8)
      {
        j = k;
        break;
      }
    } // endfor

    if(j==-1)
    {
      cout << "error in TTriaIsoparametric::GetOrigFromRef(2)" << endl;
      return;
    }

    // cout << "j= " << j << endl;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;

    // cout << "(" << Xi << ", " << Eta << ") = ";
    // cout << "(" << X[i] << "," << Y[i] << ") -> ";
    for(k=0;k<N_AuxPoints;k++)
    {
      // cout << "fct: " << FctValues[j][k] << endl;
      X[i] += XDistance[k] * FctValues[j][k];
      Y[i] += YDistance[k] * FctValues[j][k];

      a11 += XDistance[k] * XiDerValues[j][k];
      a21 += XDistance[k] * EtaDerValues[j][k];

      a12 += YDistance[k] * XiDerValues[j][k];
      a22 += YDistance[k] * EtaDerValues[j][k];
    }

    absdetjk[i] = std::abs(a11*a22 - a12*a21);

    // cout << "(" << X[i] << "," << Y[i] << ")" << endl;
  } // endfor i
*/
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

    dx1 = xc1;
    dx2 = xc2;
    
    dy1 = yc1;
    dy2 = yc2;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;

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

    // cout << "2-x: " << X[i] << " y: " << Y[i] << " z: " << Z[i] << endl;
    // cout << "3-x: " << Xi << " " << Eta << " " << Zeta;
    // cout << " " << dx1 << " y: " << dx2 << " z: " << dx3 << endl;

    // cout << endl;
    // cout << Xi << " " << Eta << " " << Zeta << endl;
    // cout << X[i] << " " << Y[i] << " " << Z[i] << endl;
    // cout << dx1 << " " << dx2 << " " << dx3 << endl;
    // cout << dy1 << " " << dy2 << " " << dy3 << endl;
    // cout << dz1 << " " << dz2 << " " << dz3 << endl;
  }
}

void TTriaIsoparametric::GetOrigFromRef(const TQuadFormula& qf_ref,
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
    
    double dx1 = xc1;
    double dx2 = xc2;
    
    double dy1 = yc1;
    double dy2 = yc2;
    
    double x = xc0 + xc1*Xi + xc2*Eta;
    double y = yc0 + yc1*Xi + yc2*Eta;
    
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
    double absdet = std::abs(dx1*dy2- dx2*dy1);
    double weight = qf_ref.get_weight(i) * absdet;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y}});
  }
}

/** transfer from original element to reference element */
void TTriaIsoparametric::GetRefFromOrig(double, double, double &, double &)
  const
{
  cout << "not implemented yet 1 " << endl;
}

/** calculate functions and derivatives from reference element
    to original element */
void TTriaIsoparametric::GetOrigValues(double xi, double eta, int N_BaseFunct,
                                       const double *uref,
                                       const double *uxiref,
                                       const double *uetaref,
                                       double *uorig, double *uxorig,
                                       double *uyorig, int) const
{
/*
  int i, j, k;
  double a11, a12, a21, a22, rec_detjk;

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(std::abs(xi-XI[i])<1e-8 && std::abs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TTriaIsoparametric::GetOrigFromRef(1)" << endl;
    return;
  }

  // D00
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  a11 = xc1;
  a21 = xc2;
  a12 = yc1;
  a22 = yc2;

  for(k=0;k<N_AuxPoints;k++)
  {
    a11 += XDistance[k] * XiDerValues[j][k];
    a21 += XDistance[k] * EtaDerValues[j][k];

    a12 += YDistance[k] * XiDerValues[j][k];
    a22 += YDistance[k] * EtaDerValues[j][k];
  } // endfor k

  rec_detjk = 1 / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
  } // endfor i
*/
  int i, j, k;
  double dx1, dx2;
  double dy1, dy2;
  double rec_detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];


  dx1 = xc1;
  dx2 = xc2;
    
  dy1 = yc1;
  dy2 = yc2;
    
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

void TTriaIsoparametric::GetOrigAllDerivatives(
  double xi, double eta, int N_BaseFunct,
  const double* r00, const double* r10, const double* r01,
  const double* r20, const double* r11, const double* r02,
  double* o00, double* o10, double* o01, double* o20, double* o11, double* o02,
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
    Output::warn("TriaIsoparametric", "transformation of second order "
                 "derivatives is not yet implementd.");
  // avoid compiler warnings
  (void)o20;
  (void)o11;
  (void)o02;
}


/** calculate functions and derivatives from reference element
    to original element */
void TTriaIsoparametric::GetOrigValues(int joint, double zeta, int N_BaseFunct,
                                       const double *uref,const double *uxiref,
                                       const double *uetaref,
                                       double *uorig, double *uxorig,
                                       double *uyorig, int) const
{
  int i, k;
  double a11, a12, a21, a22, rec_detjk;
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  double valxi[MaxN_BaseFunctions2D];
  double valeta[MaxN_BaseFunctions2D];
  double *values[1];

//  switch(joint)
//  {
//    case 0:
//      xi = 0.5*(1+zeta); eta = 0;
//    break;
//
//    case 1:
//      xi = 0.5*(1-zeta); eta = 0.5*(1+zeta);
//    break;
//
//    case 2:
//      xi = 0; eta = 0.5*(1-zeta);
//    break;
//  }

  a11 = xc1;
  a21 = xc2;
  a12 = yc1;
  a22 = yc2;

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

  rec_detjk = 1. / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
    uorig[i] = uref[i];
  } // endfor i
}

void TTriaIsoparametric::SetCell(const TBaseCell *cell)
{
  int i, j;
  const TBoundComp2D *comp;
  JointType type;
  double t0, t1, t, dt;
  double xa, ya, xe, ye, xm, ym, xp, yp, dx, dy;
//  int compid;
  int *JointDOF;
  int N_Vertices;
  TVertex **Vertices;
#ifdef __3D__
  double z[3], zp;
#endif

  N_AuxPoints = 0;

  Cell = cell;

  FEDescriptor fedesc(FEDescFromOrder[ApproximationOrder]);
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);

  // cout << endl;
  for(i=0;i<3;i++)
  {
#ifdef __3D__
    Cell->GetVertex(i)->GetCoords(x[i], y[i], z[i]);
#else
    Cell->GetVertex(i)->GetCoords(x[i], y[i]);
#endif
    // cout << setw(20) << x[i] << setw(20) << y[i] << endl;
  }

  xc0=x[0];
  xc1=x[1]-x[0];
  xc2=x[2]-x[0];

  yc0=y[0];
  yc1=y[1]-y[0];
  yc2=y[2]-y[0];

  for(i=0;i<3;i++)
  {
    // check whether the joint i is curved
    auto joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryEdge || type == InterfaceJoint)
    {
//       cout << "joint: " << i << endl;
      xa = x[i];
      ya = y[i];
      xe = x[(i+1) % 3];
      ye = y[(i+1) % 3];
      if(type == BoundaryEdge)
      {
        auto boundedge = (const TBoundEdge *)(joint);
        comp = boundedge->GetBoundComp();
//        compid = comp->GetID();
        boundedge->GetParameters(t0, t1);
      }
      else
      {
        auto interface = (const TInterfaceJoint *)(joint);
        comp = interface->GetBoundComp();
//        compid = comp->GetID();

        if(Cell == interface->GetNeighbour(0))
          interface->GetParameters(t0, t1); // correct order
        else
          interface->GetParameters(t1, t0); // backward order
      }

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

        // cout << "m: (" << xm << ", " << ym << ")" << endl;
        // cout << "p: (" << xp << ", " << yp << ")" << endl;
        // cout << "d: (" << xp-xm << ", " << yp-ym << ")" << endl;

        if(std::abs(xp-xm) > 1e-8 || std::abs(yp-ym) > 1e-8)
        {
          // curved boundary
          XDistance[N_AuxPoints] = xp - xm;
          YDistance[N_AuxPoints] = yp - ym;
          IntAux[N_AuxPoints] = JointDOF[j];
          // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
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
            N_Vertices = ((const TIsoInterfaceJoint *)joint)->GetN_Vertices();
            Vertices = ((const TIsoInterfaceJoint *)joint)->GetVertices();
          break;

          case IsoBoundEdge:
            N_Vertices = ((const TIsoBoundEdge *)joint)->GetN_Vertices();
            Vertices = ((const TIsoBoundEdge *)joint)->GetVertices();
          break;

          case IsoJointEqN:
            N_Vertices = ((const TIsoJointEqN *)joint)->GetN_Vertices();
            Vertices = ((const TIsoJointEqN *)joint)->GetVertices();
          break;
	  default:
	    break;
        } // endswitch

        // cout << "N_AuxVertices: " << N_Vertices << endl;
        ApproximationOrder = N_Vertices+1;

        JointDOF = fedesc.GetJointDOF(i);

        xa = x[i];
        ya = y[i];
        xe = x[(i+1) % 3];
        ye = y[(i+1) % 3];

//         dt = (t1-t0)/ApproximationOrder;
        dx = (xe-xa)/ApproximationOrder;
        dy = (ye-ya)/ApproximationOrder;
        for(j=1;j<ApproximationOrder;j++)
        {
//           t = t0 + j * dt;
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
	    
	    //if(XDistance[N_AuxPoints]<=0)
	     // cout <<  " xa "  << xa <<    " xe " << xe <<  " xp " << xp<< " XDist " << XDistance[N_AuxPoints] << endl;
	      
            // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
            N_AuxPoints++;
          } // endif 
        } // endfor
      } // endif
    } // endelse
  } // endfor

  if(N_AuxPoints && ApproximationOrder==3)
  {
    xm = 0;
    ym = 0;
    for(i=0;i<N_AuxPoints;i++)
    {
      xm += XDistance[i];
      ym += YDistance[i];
    }
    XDistance[N_AuxPoints] = xm/4;
    YDistance[N_AuxPoints] = ym/4;
    IntAux[N_AuxPoints] = 5;
    N_AuxPoints++;
  }

  if(N_AuxPoints)
  {
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
} // TTriaIsoparametric::SetCell

void TTriaIsoparametric::SetQuadFormula(QuadratureFormula_type formula)
{
  quadrature_formula = std::unique_ptr<TQuadFormula>(new TQuadFormula(formula));
}

void TTriaIsoparametric::PiolaMapOrigFromRef(double, double, int,
                                             const double *, double *) const
{
  ErrThrow("TTriaIsoparametric::PiolaMapOrigFromRef not implemented");
}
    
void TTriaIsoparametric::PiolaMapOrigFromRef(double, double, int,
                                             const double *, const double *,
                                             const double *, double *, double *)
  const
{
  ErrThrow("TTriaIsoparametric::PiolaMapOrigFromRef not implemented");
}

