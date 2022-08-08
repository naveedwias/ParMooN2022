#include "HexaIsoparametric.h"
#include "BoundFace.h"
#include "FEDescriptor.h"
#include "GridCell.h"
#include "IsoBoundFace.h"
#include "IsoJointEqN.h"
#include "MooNMD_Io.h"
#include "QuadFormula.h"

#include <cmath>

constexpr BaseFunction_type BaseFunctFromOrder[] = { 
                BF_C_H_Q0_3D, BF_C_H_Q1_3D, BF_C_H_Q2_3D, BF_C_H_Q3_3D,
                BF_C_H_Q4_3D, BF_C_H_Q00_3D};

constexpr FEDescriptor_type FEDescFromOrder[] = { 
                FE_C_H_Q0_3D, FE_C_H_Q1_3D, FE_C_H_Q2_3D, FE_C_H_Q3_3D,
                FE_C_H_Q4_3D, FE_C_H_Q00_3D};

/** constuctor */
THexaIsoparametric::THexaIsoparametric()
{
}

/** transfer from reference element to original element */
void THexaIsoparametric::GetOrigFromRef(double xi, double eta, double zeta,
                                        double &X, double &Y, double &Z) const
{
  int j;

  j = -1;
  unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
  for(unsigned int i=0;i<N_QuadPoints;i++)
  {
    auto p = quadrature_formula->get_point(i);
    if(std::abs(xi-p.x)<1e-8 && std::abs(eta-p.y)<1e-8
       && std::abs(zeta-p.z)<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in THexaIsoparametric::GetOrigFromRef(1)" << endl;
    exit(-1);
    return;
  }

  X = xc0 + xc1*xi + xc2*eta + xc3*zeta + xc4*xi*eta
          + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta + yc4*xi*eta
          + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta + zc4*xi*eta
          + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;

  for(int i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
    Z += ZDistance[i] * FctValues[j][i];
  }
}

/** transfer a set of point from reference to original element */
void THexaIsoparametric::GetOrigFromRef(int N_Points,
                const double *xi, const double *eta, const double *zeta,
                double *X, double *Y, double *Z) const
{
  int i, j, k;
  double Xi, Eta, Zeta;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double AuxVector[4*MaxN_BaseFunctions3D];

  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];

    dx1=xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
    dx2=xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
    dx3=xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
    
    dy1=yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
    dy2=yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
    dy3=yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
    
    dz1=zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    dz2=zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    dz3=zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;

    bf.GetDerivatives(MultiIndex3D::D000, Xi, Eta, Zeta, AuxVector);
    bf.GetDerivatives(MultiIndex3D::D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
    bf.GetDerivatives(MultiIndex3D::D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
    bf.GetDerivatives(MultiIndex3D::D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);

    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];
      Z[i] += ZDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
    }
  }
}

void THexaIsoparametric::GetOrigFromRef(const TQuadFormula& qf_ref,
                                        TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  double AuxVector[4*MaxN_BaseFunctions3D];
  BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
  for(unsigned int i=0;i<n_points;i++)
  {
    auto p = qf_ref.get_point(i);
    double Xi = p.x;
    double Eta = p.y;
    double Zeta = p.z;

    double dx1 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
    double dx2 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
    double dx3 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
    
    double dy1 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
    double dy2 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
    double dy3 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
    
    double dz1 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    double dz2 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    double dz3 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    
    double x = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    double y = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    double z = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;

    bf.GetDerivatives(MultiIndex3D::D000, Xi, Eta, Zeta, AuxVector);
    bf.GetDerivatives(MultiIndex3D::D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
    bf.GetDerivatives(MultiIndex3D::D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
    bf.GetDerivatives(MultiIndex3D::D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);

    for(int k=0;k<N_AuxPoints;k++)
    {
      int j = IntAux[k];
      x += XDistance[k] * AuxVector[j];
      y += YDistance[k] * AuxVector[j];
      z += ZDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
    }

    double absdet = std::abs(dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                           - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2);
    double weight = qf_ref.get_weight(i) * absdet;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y, z}});
  }
}

/** transfer from original element to reference element */
void THexaIsoparametric::GetRefFromOrig(double X, double Y, double Z,
                                        double &xi, double &eta, double &zeta)
  const
{
  double xt = X - xc0;
  double yt = Y - yc0;
  double zt = Z - zc0;
  double t0, t1, t2;
  double xi0, eta0, zeta0;
  double recdetaffine;
  double eps=1e-14;

  recdetaffine = 1/( xc1*yc2*zc3-xc1*yc3*zc2-yc1*xc2*zc3
                    +yc1*xc3*zc2+zc1*xc2*yc3-zc1*xc3*yc2 );

  t0 = (-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
  t1 = ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
  t2 = (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;

  xi0 = t0;
  eta0 = t1;
  zeta0 = t2;

  xi = xi0+10;
  eta = eta0+10;
  zeta = zeta0+10;
  while( std::abs(xi-xi0)+std::abs(eta-eta0)+std::abs(zeta-zeta0) > eps)
  {
    xi = xi0;
    eta = eta0;
    zeta = zeta0;

    xt = xc4*xi*eta + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
    yt = yc4*xi*eta + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
    zt = zc4*xi*eta + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;

    xi0   = -t0+(-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
    eta0  = -t1+ ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
    zeta0  = -t2+ (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;
  }
  xi = xi0;
  eta = eta0;
  zeta = zeta0;
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaIsoparametric::GetOrigValues(double xi, double eta, double zeta,
                                       int N_BaseFunct, const double *uref,
                                       const double *uxiref,
                                       const double *uetaref,
                                       const double *uzetaref,
                                       double *uorig, double *uxorig,
                                       double *uyorig, double *uzorig,
                                       int BaseVectDim) const
{
  if(BaseVectDim != 1)
  {
    ErrThrow("HexaIsoparametric for vector valued functions not "
             "implemented. (Piola transform)");
  }
  int j,k;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  // D000
  for(int i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  j = -1;
  unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
  for(unsigned int i=0;i<N_QuadPoints;i++)
  {
    auto p = quadrature_formula->get_point(i);
    if(std::abs(xi-p.x)<1e-8 && std::abs(eta-p.y)<1e-8
       && std::abs(zeta-p.z)<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    ErrThrow("THexaIsoparametric::GetOrigFromRef(3)");
  }

  dx1=xc1 + xc4*eta + xc5*zeta + xc7*eta*zeta;
  dx2=xc2 + xc4*xi + xc6*zeta + xc7*xi*zeta;
  dx3=xc3 + xc5*xi + xc6*eta + xc7*xi*eta;
              
  dy1=yc1 + yc4*eta + yc5*zeta + yc7*eta*zeta;
  dy2=yc2 + yc4*xi + yc6*zeta + yc7*xi*zeta;
  dy3=yc3 + yc5*xi + yc6*eta + yc7*xi*eta;
  
  dz1=zc1 + zc4*eta + zc5*zeta + zc7*eta*zeta;
  dz2=zc2 + zc4*xi + zc6*zeta + zc7*xi*zeta;
  dz3=zc3 + zc5*xi + zc6*eta + zc7*xi*eta;
  
  for(k=0;k<N_AuxPoints;k++)
  {
    dx1 += XDistance[k] * XiDerValues[j][k];
    dx2 += XDistance[k] * EtaDerValues[j][k];
    dx3 += XDistance[k] * ZetaDerValues[j][k];

    dy1 += YDistance[k] * XiDerValues[j][k];
    dy2 += YDistance[k] * EtaDerValues[j][k];
    dy3 += YDistance[k] * ZetaDerValues[j][k];

    dz1 += ZDistance[k] * XiDerValues[j][k];
    dz2 += ZDistance[k] * EtaDerValues[j][k];
    dz3 += ZDistance[k] * ZetaDerValues[j][k];
  }
  double detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  double rec_detjk = 1/detjk;

  // D100, D010 and D001
  for(int i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((dy2*dz3 - dy3*dz2)*uxiref[i] +
                 (dy3*dz1 - dy1*dz3)*uetaref[i] +
                 (dy1*dz2 - dy2*dz1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((dx3*dz2 - dx2*dz3)*uxiref[i] +
                 (dx1*dz3 - dx3*dz1)*uetaref[i] +
                 (dx2*dz1 - dx1*dz2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((dx2*dy3 - dx3*dy2)*uxiref[i] +
                 (dx3*dy1 - dx1*dy3)*uetaref[i] + 
                 (dx1*dy2 - dx2*dy1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void THexaIsoparametric::GetOrigAllDerivatives(
    double xi, double eta, double zeta, int N_BaseFunct,
    const double* refD000, const double* refD100,
    const double* refD010, const double* refD001,
    const double* refD200, const double* refD110,
    const double* refD101, const double* refD020,
    const double* refD011, const double* refD002,
    double* origD000, double* origD100, double* origD010, double* origD001,
    double* origD200, double* origD110, double* origD101, double* origD020,
    double* origD011, double* origD002, int BaseVectDim) const
{
  GetOrigValues(xi, eta, zeta, N_BaseFunct, refD000, refD100, refD010, refD001,
                origD000, origD100, origD010, origD001, BaseVectDim);
  bool all_zero_2nd_derivatives = true;
  for(int i=0;i<N_BaseFunct;i++)
  {
    if(refD200[i] != 0. || refD110[i] != 0. || refD101[i] != 0. 
       || refD020[i] != 0. || refD011[i] != 0. || refD002[i] != 0.)
    {
      all_zero_2nd_derivatives = false;
      break;
    }
  } // endfor i
  if(!all_zero_2nd_derivatives)
    Output::warn("THexaIsoparametric", "transformation of second order "
                 "derivatives is not yet implemented.");
  // avoid compiler warnings
  (void)origD200;
  (void)origD110;
  (void)origD101;
  (void)origD020;
  (void)origD011;
  (void)origD002;
}

void THexaIsoparametric::SetCell(const TBaseCell *cell)
{
  int i, j, k, l;
  const TJoint *joint;
  JointType type;

  const TShapeDesc *ShapeDesc;
  const int *TmpFaceVertex, *TmpLen;
  int MaxLen;
  int *JointDOF;
  double Param1[4], Param2[4];
  TBoundFace *bdface;
  double dt, ds;
  double T, S, xm, ym, zm, xp, yp, zp;
  const TVertex * const * Vertices;
  const TVertex * const * AuxVertices=nullptr;
  double x[4], y[4], z[4], factor;
  double LinComb[4];
  int CurvedJoint = 0;

  N_AuxPoints = 0;
  Cell = cell;
  Vertices = ((TGridCell*)Cell)->GetVertices();
  ShapeDesc = Cell->GetShapeDesc();
  ShapeDesc->GetFaceVertex(TmpFaceVertex, TmpLen, MaxLen);

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
  Cell->GetVertex(4)->GetCoords(x4, y4, z4);
  Cell->GetVertex(5)->GetCoords(x5, y5, z5);
  Cell->GetVertex(6)->GetCoords(x6, y6, z6);
  Cell->GetVertex(7)->GetCoords(x7, y7, z7);

  xc0 = ( x0 + x1 + x3 + x2 + x4 + x5 + x7 + x6) * 0.125;
  xc1 = (-x0 + x1 - x3 + x2 - x4 + x5 - x7 + x6) * 0.125;
  xc2 = (-x0 - x1 + x3 + x2 - x4 - x5 + x7 + x6) * 0.125;
  xc3 = (-x0 - x1 - x3 - x2 + x4 + x5 + x7 + x6) * 0.125;
  xc4 = ( x0 - x1 - x3 + x2 + x4 - x5 - x7 + x6) * 0.125;
  xc5 = ( x0 - x1 + x3 - x2 - x4 + x5 - x7 + x6) * 0.125;
  xc6 = ( x0 + x1 - x3 - x2 - x4 - x5 + x7 + x6) * 0.125;
  xc7 = (-x0 + x1 + x3 - x2 + x4 - x5 - x7 + x6) * 0.125;

  yc0 = ( y0 + y1 + y3 + y2 + y4 + y5 + y7 + y6) * 0.125;
  yc1 = (-y0 + y1 - y3 + y2 - y4 + y5 - y7 + y6) * 0.125;
  yc2 = (-y0 - y1 + y3 + y2 - y4 - y5 + y7 + y6) * 0.125;
  yc3 = (-y0 - y1 - y3 - y2 + y4 + y5 + y7 + y6) * 0.125;
  yc4 = ( y0 - y1 - y3 + y2 + y4 - y5 - y7 + y6) * 0.125;
  yc5 = ( y0 - y1 + y3 - y2 - y4 + y5 - y7 + y6) * 0.125;
  yc6 = ( y0 + y1 - y3 - y2 - y4 - y5 + y7 + y6) * 0.125;
  yc7 = (-y0 + y1 + y3 - y2 + y4 - y5 - y7 + y6) * 0.125;

  zc0 = ( z0 + z1 + z3 + z2 + z4 + z5 + z7 + z6) * 0.125;
  zc1 = (-z0 + z1 - z3 + z2 - z4 + z5 - z7 + z6) * 0.125;
  zc2 = (-z0 - z1 + z3 + z2 - z4 - z5 + z7 + z6) * 0.125;
  zc3 = (-z0 - z1 - z3 - z2 + z4 + z5 + z7 + z6) * 0.125;
  zc4 = ( z0 - z1 - z3 + z2 + z4 - z5 - z7 + z6) * 0.125;
  zc5 = ( z0 - z1 + z3 - z2 - z4 + z5 - z7 + z6) * 0.125;
  zc6 = ( z0 + z1 - z3 - z2 - z4 - z5 + z7 + z6) * 0.125;
  zc7 = (-z0 + z1 + z3 - z2 + z4 - z5 - z7 + z6) * 0.125;

  // Output::print(Cell->GetClipBoard());
  for(i=0;i<6;i++)
  {
    joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryFace || type == IsoBoundFace || type == IsoJointEqN)
    {
      FEDescriptor fedesc(FEDescFromOrder[ApproximationOrder]);
      JointDOF = fedesc.GetJointDOF(i);

      if(type == BoundaryFace || type == IsoBoundFace)
      {
        bdface = (TBoundFace*)joint;
        bdface->GetParameters(Param1, Param2);
      }
      
      dt = ds = 1.0/ApproximationOrder;
      Vertices[TmpFaceVertex[i*MaxLen+0]]->GetCoords(x[0],y[0],z[0]);
      Vertices[TmpFaceVertex[i*MaxLen+1]]->GetCoords(x[1],y[1],z[1]);
      Vertices[TmpFaceVertex[i*MaxLen+2]]->GetCoords(x[2],y[2],z[2]);
      Vertices[TmpFaceVertex[i*MaxLen+3]]->GetCoords(x[3],y[3],z[3]);

      if(type == IsoBoundFace)
        AuxVertices = ((TIsoBoundFace *)joint)->GetVertices();

      if(type == IsoJointEqN)
        AuxVertices = ((TIsoJointEqN *)joint)->GetVertices();

      for(j=0;j<=ApproximationOrder;j++)
      {
        for(k=0;k<=ApproximationOrder;k++)
        {
          LinComb[0] = (1-j*dt)*(1-k*ds);
          LinComb[1] = (j*dt)*(1-k*ds);
          LinComb[2] = (j*dt)*(k*ds);
          LinComb[3] = (1-j*dt)*(k*ds);

          xm = 0; ym = 0; zm = 0;
          for(l=0;l<4;l++)
          {
            factor = LinComb[l];
            // cout << "LinComb: " << LinComb[l] << endl;
            xm += factor*x[l]; 
            ym += factor*y[l]; 
            zm += factor*z[l]; 
            // cout << x[l] << " " << y[l] << " " << z[l] << endl;
          }

          if(type == BoundaryFace || ApproximationOrder == 1)
          {
            if(type == BoundaryFace || type == IsoBoundFace)
              bdface->GetBoundComp()->GetXYZandTS(4, LinComb, x, y, z,
                                                  Param1, Param2,
                                                  xp, yp, zp, T, S);
            else
            {
              xp = xm; yp = ym; zp = zm;
            }
          }
          else
          {
            if(type == IsoJointEqN)
            {
              if(joint->GetNeighbour(0) == cell)
              {
                AuxVertices[k*(ApproximationOrder+1)+j]->GetCoords(xp, yp, zp);
              }
              else
              {
                switch(joint->GetMapType())
                {
                  case 0:
                    AuxVertices[j*(ApproximationOrder+1)
                                +k]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 1:
                    AuxVertices[k*(ApproximationOrder+1)
                                +(ApproximationOrder-j)]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 2:
                    AuxVertices[(ApproximationOrder-j)*(ApproximationOrder+1)
                                +(ApproximationOrder-k)]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 3:
                    AuxVertices[(ApproximationOrder-k)*(ApproximationOrder+1)
                                +j]
                                ->GetCoords(xp, yp, zp);
                  break;
                } // endswitch
              } // endelse
            } // endelse
            else
            {
              AuxVertices[k*(ApproximationOrder+1)+j]->GetCoords(xp, yp, zp);
            }
          } // endelse

          if(std::abs(xp-xm) > 1e-8 || std::abs(yp-ym) > 1e-8 || std::abs(zp-zm) > 1e-8)
          {
            // cout << "j: " << j << " k: " << k << endl;
            // cout << "I: " << N_AuxPoints << " T:" << T << " S:" << S << endl;
            // cout << "Straight: " << xm << " " << ym << " " << zm << endl;
            // cout << "Boundary: " << xp << " " << yp << " " << zp << endl;
            // cout << "Diff:" << xp-xm << " " << yp-ym << " " << zp-zm << endl;
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            ZDistance[N_AuxPoints] = zp - zm;
            IntAux[N_AuxPoints] = JointDOF[k*(ApproximationOrder+1)+j];
            // cout << "jd: " << JointDOF[k*(ApproximationOrder+1)+j] << endl;
            N_AuxPoints++;

            CurvedJoint = i;
          } // endif
        } // endfor k
      } // endfor j
    } // endif
  } // endfor i

  // cout << "cell nr: " << cell->GetClipBoard() << " " << N_AuxPoints << endl;

  if(N_AuxPoints)
  {
    if(ApproximationOrder == 2)
    {
      XDistance[N_AuxPoints] = XDistance[0]*0.5;
      YDistance[N_AuxPoints] = YDistance[0]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[0]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*0.5;
      YDistance[N_AuxPoints] = YDistance[1]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[1]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*0.5;
      YDistance[N_AuxPoints] = YDistance[2]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[2]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*0.5;
      YDistance[N_AuxPoints] = YDistance[3]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[3]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*0.5;
      YDistance[N_AuxPoints] = YDistance[4]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[4]*0.5;
      N_AuxPoints++;

      N_AuxPoints -= 5;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
        break;

        case 5:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
        break;

        /*
        default:
          Output::print("This case has not been implemented yet!");
        break;
        */
      } 
    } // ApproximationOrder == 2

    if(ApproximationOrder == 3)
    {
      XDistance[N_AuxPoints] = XDistance[0]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[0]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[1]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[2]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[3]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[4]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[4]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[5]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[5]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[6]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[6]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[7]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[7]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[8]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[8]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[9]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[9]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[10]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[10]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[11]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[11]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[0]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[1]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[2]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[3]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[4]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[4]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[5]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[5]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[6]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[6]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[7]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[7]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[8]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[8]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[9]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[9]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[10]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[10]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[11]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[11]*1.0/3;
      N_AuxPoints++;

      N_AuxPoints -= 24;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
        break;
        
        default:
          Output::print("This case has not been implemented yet!");
        break;
      }
    }  // ApproximationOrder == 3

    if(ApproximationOrder == 4)
    {
      XDistance[N_AuxPoints] = XDistance[0]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*1.0/4;
      N_AuxPoints++;

      N_AuxPoints -= 63;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

        break;
        
        default:
          Output::print("This case has not been implemented yet!");
        break;
      }
    }  // ApproximationOrder == 4
  } // N_AuxPoints > 0

  if(N_AuxPoints)
  {
    BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
    unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
    for(unsigned int i=0;i<N_QuadPoints;i++)
    {
      auto p = quadrature_formula->get_point(i);
      bf.GetDerivatives(MultiIndex3D::D000, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]]; 

      bf.GetDerivatives(MultiIndex3D::D100, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf.GetDerivatives(MultiIndex3D::D010, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf.GetDerivatives(MultiIndex3D::D001, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        ZetaDerValues[i][j] = DoubleAux[IntAux[j]]; 
    } // endfor i
  } // endif
}

void THexaIsoparametric::SetQuadFormula(QuadratureFormula_type formula)
{
  QuadFormula_type = formula;
  quadrature_formula = std::unique_ptr<TQuadFormula>(new TQuadFormula(formula));
}

void THexaIsoparametric::PiolaMapOrigFromRef(double, double, double, int,
                                             const double *, double *) const
{
  ErrThrow("THexaIsoparametric::PiolaMapOrigFromRef not implemented");
}
