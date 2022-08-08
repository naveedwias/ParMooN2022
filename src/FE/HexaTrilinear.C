#include "HexaTrilinear.h"
#include "BaseCell.h"
#include "MooNMD_Io.h"
#include "QuadFormula.h"

#include <cmath>

/** constuctor */
THexaTrilinear::THexaTrilinear()
{
}

/** transfer from reference element to original element */
void THexaTrilinear::GetOrigFromRef(double xi, double eta, double zeta,
                                    double &X, double &Y, double &Z) const
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta + xc4*xi*eta
          + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta + yc4*xi*eta
          + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta + zc4*xi*eta
          + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;
}

/** transfer a set of point from reference to original element */
void THexaTrilinear::GetOrigFromRef(int N_Points, const double *xi,
                                    const double *eta, const double *zeta,
                                    double *X, double *Y, double *Z) const
{
  double Xi, Eta, Zeta;
  for(int i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];

    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;
  }
}

void THexaTrilinear::GetOrigFromRef(const TQuadFormula& qf_ref,
                                    TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  for(unsigned int i = 0; i < n_points; ++i)
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
    
    double absdet = std::abs(dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                           - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2);
    
    double x = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    double y = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    double z = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;
    double weight = qf_ref.get_weight(i) * absdet;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y, z}});
  }
}

/** transfer from original element to reference element */
void THexaTrilinear::GetRefFromOrig(double X, double Y, double Z,
                                    double &xi, double &eta, double &zeta) const
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
void THexaTrilinear::GetOrigValues(double xi, double eta, double zeta,
                                   int N_BaseFunct,
                                   const double *uref, const double *uxiref,
                                   const double *uetaref,
                                   const double *uzetaref,
                                   double *uorig, double *uxorig,
                                   double *uyorig, double *uzorig,
                                   int _BaseVectDim) const
{
  int i;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  if(_BaseVectDim != 1)
  {
    ErrThrow("mixed finite elements with trilinear reference transformation "
             "not available");
  }
  
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  dx1=xc1 + xc4*eta + xc5*zeta + xc7*eta*zeta;
  dx2=xc2 + xc4*xi + xc6*zeta + xc7*xi*zeta;
  dx3=xc3 + xc5*xi + xc6*eta + xc7*xi*eta;
              
  dy1=yc1 + yc4*eta + yc5*zeta + yc7*eta*zeta;
  dy2=yc2 + yc4*xi + yc6*zeta + yc7*xi*zeta;
  dy3=yc3 + yc5*xi + yc6*eta + yc7*xi*eta;
  
  dz1=zc1 + zc4*eta + zc5*zeta + zc7*eta*zeta;
  dz2=zc2 + zc4*xi + zc6*zeta + zc7*xi*zeta;
  dz3=zc3 + zc5*xi + zc6*eta + zc7*xi*eta;

  double detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2 
               - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  double rec_detjk = 1/detjk;

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((dy2*dz3 - dy3*dz2)*uxiref[i] + (dy3*dz1 - dy1*dz3)*uetaref[i] + (dy1*dz2 - dy2*dz1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((dx3*dz2 - dx2*dz3)*uxiref[i] + (dx1*dz3 - dx3*dz1)*uetaref[i] + (dx2*dz1 - dx1*dz2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((dx2*dy3 - dx3*dy2)*uxiref[i] + (dx3*dy1 - dx1*dy3)*uetaref[i] + (dx1*dy2 - dx2*dy1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void THexaTrilinear::GetOrigAllDerivatives(double xi, double eta, double zeta,
                                           int N_BaseFunct,
                                           const double* refD000,
                                           const double* refD100,
                                           const double* refD010,
                                           const double* refD001,
                                           const double* refD200,
                                           const double* refD110,
                                           const double* refD101,
                                           const double* refD020,
                                           const double* refD011,
                                           const double* refD002,
                                           double* origD000, double* origD100,
                                           double* origD010, double* origD001,
                                           double* origD200, double* origD110,
                                           double* origD101, double* origD020,
                                           double* origD011, double* origD002,
                                           int BaseVectDim) const
{
  GetOrigValues(xi, eta, zeta, N_BaseFunct, refD000, refD100, refD010, refD001,
                origD000, origD100, origD010, origD001, BaseVectDim);
  // avoid compiler warnings
  (void)refD200;
  (void)refD110;
  (void)refD101;
  (void)refD020;
  (void)refD011;
  (void)refD002;
  (void)origD200;
  (void)origD110;
  (void)origD101;
  (void)origD020;
  (void)origD011;
  (void)origD002;
  ErrThrow("TTetraIsoparametric: transformation of second order derivatives is "
           "not yet implemented.");
}

void THexaTrilinear::SetCell(const TBaseCell *cell)
{

  Cell = cell;

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
  Cell->GetVertex(4)->GetCoords(x4, y4, z4);
  Cell->GetVertex(5)->GetCoords(x5, y5, z5);
  Cell->GetVertex(6)->GetCoords(x6, y6, z6);
  Cell->GetVertex(7)->GetCoords(x7, y7, z7);

/*
  Output::print("0: ", Cell->GetVertex(0));
  Output::print("1: ", Cell->GetVertex(1));
  Output::print("2: ", Cell->GetVertex(2));
  Output::print("3: ", Cell->GetVertex(3));
  Output::print("4: ", Cell->GetVertex(4));
  Output::print("5: ", Cell->GetVertex(5));
  Output::print("6: ", Cell->GetVertex(6));
  Output::print("7: ", Cell->GetVertex(7));
*/

  xc0 = (x0 + x1 + x3 + x2 + x4 + x5 + x7 + x6) * 0.125;
  xc1 = (-x0 + x1 - x3 + x2 - x4 + x5 - x7 + x6) * 0.125;
  xc2 = (-x0 - x1 + x3 + x2 - x4 - x5 + x7 + x6) * 0.125;
  xc3 = (-x0 - x1 - x3 - x2 + x4 + x5 + x7 + x6) * 0.125;
  xc4 = (x0 - x1 - x3 + x2 + x4 - x5 - x7 + x6) * 0.125;
  xc5 = (x0 - x1 + x3 - x2 - x4 + x5 - x7 + x6) * 0.125;
  xc6 = (x0 + x1 - x3 - x2 - x4 - x5 + x7 + x6) * 0.125;
  xc7 = (-x0 + x1 + x3 - x2 + x4 - x5 - x7 + x6) * 0.125;

  // cout << "x:" << xc4 << xc5 << xc6 << xc7 << endl;

  yc0 = (y0 + y1 + y3 + y2 + y4 + y5 + y7 + y6) * 0.125;
  yc1 = (-y0 + y1 - y3 + y2 - y4 + y5 - y7 + y6) * 0.125;
  yc2 = (-y0 - y1 + y3 + y2 - y4 - y5 + y7 + y6) * 0.125;
  yc3 = (-y0 - y1 - y3 - y2 + y4 + y5 + y7 + y6) * 0.125;
  yc4 = (y0 - y1 - y3 + y2 + y4 - y5 - y7 + y6) * 0.125;
  yc5 = (y0 - y1 + y3 - y2 - y4 + y5 - y7 + y6) * 0.125;
  yc6 = (y0 + y1 - y3 - y2 - y4 - y5 + y7 + y6) * 0.125;
  yc7 = (-y0 + y1 + y3 - y2 + y4 - y5 - y7 + y6) * 0.125;

  // cout << "y:" << yc4 << yc5 << yc6 << yc7 << endl;

  zc0 = (z0 + z1 + z3 + z2 + z4 + z5 + z7 + z6) * 0.125;
  zc1 = (-z0 + z1 - z3 + z2 - z4 + z5 - z7 + z6) * 0.125;
  zc2 = (-z0 - z1 + z3 + z2 - z4 - z5 + z7 + z6) * 0.125;
  zc3 = (-z0 - z1 - z3 - z2 + z4 + z5 + z7 + z6) * 0.125;
  zc4 = (z0 - z1 - z3 + z2 + z4 - z5 - z7 + z6) * 0.125;
  zc5 = (z0 - z1 + z3 - z2 - z4 + z5 - z7 + z6) * 0.125;
  zc6 = (z0 + z1 - z3 - z2 - z4 - z5 + z7 + z6) * 0.125;
  zc7 = (-z0 + z1 + z3 - z2 + z4 - z5 - z7 + z6) * 0.125;

  // cout << "z:" << zc4 << zc5 << zc6 << zc7 << endl;
}

void THexaTrilinear::PiolaMapOrigFromRef(double, double, double, int,
                                         const double *, double *) const
{
  ErrThrow("Piola Map for HexaTrilinear reference map not yet implemented");
}
