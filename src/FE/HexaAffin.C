#include <cmath>

#include "HexaAffin.h"
#include "BaseCell.h"
#include "LinAlg.h"
#include "MooNMD_Io.h"
#include "QuadFormula.h"

/** constuctor */
THexaAffin::THexaAffin()
{
}

/** transfer from reference element to original element */
void THexaAffin::GetOrigFromRef(double xi, double eta, double zeta, double &X,
                                double &Y, double &Z) const
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
}

/** transfer a set of point from reference to original element */
void THexaAffin::GetOrigFromRef(int N_Points, const double *xi,
                                const double *eta, const double *zeta,
                                double *X, double *Y, double *Z) const
{
  double Xi, Eta, Zeta;
  for(int i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;
  }
}

void THexaAffin::GetOrigFromRef(const TQuadFormula& qf_ref,
                                TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  double absdet = std::abs(detjk);
  for(unsigned int i = 0; i < n_points; ++i)
  {
    auto p = qf_ref.get_point(i);
    double weight = qf_ref.get_weight(i) * absdet;
    double x = xc0 + xc1*p.x + xc2*p.y + xc3*p.z;
    double y = yc0 + yc1*p.x + yc2*p.y + yc3*p.z;
    double z = zc0 + zc1*p.x + zc2*p.y + zc3*p.z;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y, z}});
  }
}

/** transfer from original element to reference element */
void THexaAffin::GetRefFromOrig(double X, double Y, double Z, double &xi,
                                double &eta, double &zeta) const
{
  double xt=(X - xc0)/detjk;
  double yt=(Y - yc0)/detjk;
  double zt=(Z - zc0)/detjk;

  xi  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
  eta = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
  zeta = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaAffin::GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref,const  double *uzetaref,
                               double *uorig, double *uxorig, double *uyorig, double *uzorig,
                               int _BaseVectDim) const
{
  if(_BaseVectDim == 1)
  {
    // D000
    for(int i=0;i<N_BaseFunct;i++)
      uorig[i] = uref[i];
    
    // D100, D010 and D001
    for(int i=0;i<N_BaseFunct;i++)
    {
      uxorig[i] = ( (yc2*zc3 - yc3*zc2)*uxiref[i] 
                  + (yc3*zc1 - yc1*zc3)*uetaref[i] 
                  + (yc1*zc2 - yc2*zc1)*uzetaref[i] ) / detjk;
      uyorig[i] = ( (xc3*zc2 - xc2*zc3)*uxiref[i] 
                  + (xc1*zc3 - xc3*zc1)*uetaref[i] 
                  + (xc2*zc1 - xc1*zc2)*uzetaref[i] ) / detjk;
      uzorig[i] = ( (xc2*yc3 - xc3*yc2)*uxiref[i] 
                  + (xc3*yc1 - xc1*yc3)*uetaref[i] 
                  + (xc1*yc2 - xc2*yc1)*uzetaref[i] ) / detjk;
    } // endfor i
  }
  else if(_BaseVectDim == 3)
  {
    // D000
    PiolaMapOrigFromRef(xi, eta, zeta, N_BaseFunct, uref, uxiref, uetaref,
                        uzetaref, uorig, uxorig, uyorig, uzorig);
  }
  else
  {
    ErrThrow("unknown basis function dimension");
  }
}

void THexaAffin::GetOrigAllDerivatives(double xi, double eta, double zeta,
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
  if(!ready_for_transforming_2nd_derivatives)
  {
    double GeoData[6][6];
    // reset matrices
    std::fill((double*)&GeoData[0], (double*)&GeoData[0]+36, 0.);
    std::fill((double*)&Eye[0], (double*)&Eye[0]+36, 0.);
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    Eye[3][3] = 1;
    Eye[4][4] = 1;
    Eye[5][5] = 1;
    
    GeoData[0][0] = xc1*xc1;
    GeoData[0][1] = 2*xc1*yc1;
    GeoData[0][2] = 2*xc1*zc1;
    GeoData[0][3] = yc1*yc1;
    GeoData[0][4] = 2*yc1*zc1;
    GeoData[0][5] = zc1*zc1;
    
    GeoData[1][0] = xc1*xc2;
    GeoData[1][1] = yc1*xc2 + xc1*yc2; 
    GeoData[1][2] = xc1*zc2 + xc2*zc1;
    GeoData[1][3] = yc1*yc2;
    GeoData[1][4] = yc1*zc2 + yc2*zc1;
    GeoData[1][5] = zc1*zc2;
    
    GeoData[2][0] = xc1*xc3;
    GeoData[2][1] = xc1*yc3 + xc3*yc1;
    GeoData[2][2] = xc1*zc3 + xc3*zc1;
    GeoData[2][3] = yc1*yc3;
    GeoData[2][4] = yc1*zc3 + yc3*zc1;
    GeoData[2][5] = zc1*zc3;
    
    GeoData[3][0] = xc2*xc2;
    GeoData[3][1] = 2*xc2*yc2;
    GeoData[3][2] = 2*xc2*zc2;
    GeoData[3][3] = yc2*yc2;
    GeoData[3][4] = 2*yc2*zc2;
    GeoData[3][5] = zc2*zc2;
    
    GeoData[4][0] = xc2*xc3;
    GeoData[4][1] = xc2*yc3 + xc3*yc2;
    GeoData[4][2] = xc2*zc3 + xc3*zc2;
    GeoData[4][3] = yc2*yc3;
    GeoData[4][4] = yc2*zc3 + yc3*zc2;
    GeoData[4][5] = zc2*zc3;
    
    GeoData[5][0] = xc3*xc3;
    GeoData[5][1] = 2*xc3*yc3;
    GeoData[5][2] = 2*xc3*zc3;
    GeoData[5][3] = yc3*yc3;
    GeoData[5][4] = 2*yc3*zc3;
    GeoData[5][5] = zc3*zc3;
    
    // invert GeoData and store result in Eye
    SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 6, 6, 6, 6);
    ready_for_transforming_2nd_derivatives = true;
  }
  if(BaseVectDim != 1)
  {
    ErrThrow("transformation of second derivatives of vector valued basis "
             "functions not supported yet");
  }
  for(int k = 0; k < N_BaseFunct; k++)
  {
    double r200 = refD200[k];
    double r110 = refD110[k];
    double r101 = refD101[k];
    double r020 = refD020[k];
    double r011 = refD011[k];
    double r002 = refD002[k];
    
    origD200[k] =  Eye[0][0]*r200 + Eye[0][1]*r110 + Eye[0][2]*r101
                 + Eye[0][3]*r020 + Eye[0][4]*r011 + Eye[0][5]*r002;
    origD110[k] =  Eye[1][0]*r200 + Eye[1][1]*r110 + Eye[1][2]*r101
                 + Eye[1][3]*r020 + Eye[1][4]*r011 + Eye[1][5]*r002;
    origD101[k] =  Eye[2][0]*r200 + Eye[2][1]*r110 + Eye[2][2]*r101
                 + Eye[2][3]*r020 + Eye[2][4]*r011 + Eye[2][5]*r002;
    origD020[k] =  Eye[3][0]*r200 + Eye[3][1]*r110 + Eye[3][2]*r101
                 + Eye[3][3]*r020 + Eye[3][4]*r011 + Eye[3][5]*r002;
    origD011[k] =  Eye[4][0]*r200 + Eye[4][1]*r110 + Eye[4][2]*r101
                 + Eye[4][3]*r020 + Eye[4][4]*r011 + Eye[4][5]*r002;
    origD002[k] =  Eye[5][0]*r200 + Eye[5][1]*r110 + Eye[5][2]*r101
                 + Eye[5][3]*r020 + Eye[5][4]*r011 + Eye[5][5]*r002;
  } // endfor k
}

void THexaAffin::SetCell(const TBaseCell *cell)
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

  xc0 = (x1 + x3 + x4 - x0) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;
  xc3 = (x4 - x0) * 0.5;

  yc0 = (y1 + y3 + y4 - y0) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;
  yc3 = (y4 - y0) * 0.5;

  zc0 = (z1 + z3 + z4 - z0) * 0.5;
  zc1 = (z1 - z0) * 0.5;
  zc2 = (z3 - z0) * 0.5;
  zc3 = (z4 - z0) * 0.5;

  detjk=xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2 - xc3*yc2*zc1 - xc2*yc1*zc3 - xc1*yc3*zc2;
  ready_for_transforming_2nd_derivatives = false;
}

/** Piola transformation for vectorial basis functions */
void THexaAffin::PiolaMapOrigFromRef(double , double , double , int N_Functs,
                                     const double *refD000,
                                     const double *refD100,
                                     const double *refD010,
                                     const double *refD001,
                                     double *origD000, double *origD100,
                                     double *origD010, double *origD001) const
{
  double a11 = xc1 / detjk;
  double a12 = xc2 / detjk;
  double a13 = xc3 / detjk;
  double a21 = yc1 / detjk;
  double a22 = yc2 / detjk;
  double a23 = yc3 / detjk;
  double a31 = zc1 / detjk;
  double a32 = zc2 / detjk;
  double a33 = zc3 / detjk;
  double i11 = (yc2*zc3 - yc3*zc2) / detjk;
  double i12 = (yc3*zc1 - yc1*zc3) / detjk;
  double i13 = (yc1*zc2 - yc2*zc1) / detjk;
  double i21 = (xc3*zc2 - xc2*zc3) / detjk;
  double i22 = (xc1*zc3 - xc3*zc1) / detjk;
  double i23 = (xc2*zc1 - xc1*zc2) / detjk;
  double i31 = (xc2*yc3 - xc3*yc2) / detjk;
  double i32 = (xc3*yc1 - xc1*yc3) / detjk;
  double i33 = (xc1*yc2 - xc2*yc1) / detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // three components:
    origD000[k             ] = a11 * refD000[k             ] 
                             + a12 * refD000[k +   N_Functs] 
                             + a13 * refD000[k + 2*N_Functs]; 
    origD000[k + N_Functs  ] = a21 * refD000[k             ] 
                             + a22 * refD000[k +   N_Functs] 
                             + a23 * refD000[k + 2*N_Functs]; 
    origD000[k + 2*N_Functs] = a31 * refD000[k             ] 
                             + a32 * refD000[k +   N_Functs] 
                             + a33 * refD000[k + 2*N_Functs];

    double p11 = a11 * refD100[k             ] 
               + a12 * refD100[k + N_Functs  ] 
               + a13 * refD100[k + 2*N_Functs];
    double p21 = a11 * refD010[k             ] 
               + a12 * refD010[k + N_Functs  ] 
               + a13 * refD010[k + 2*N_Functs];
    double p31 = a11 * refD001[k             ] 
               + a12 * refD001[k + N_Functs  ] 
               + a13 * refD001[k + 2*N_Functs];
    origD100[k] = i11*p11 + i12*p21 + i13*p31;
    origD010[k] = i21*p11 + i22*p21 + i23*p31;
    origD001[k] = i31*p11 + i32*p21 + i33*p31;
    
    double p12 = a21 * refD100[k             ] 
               + a22 * refD100[k + N_Functs  ] 
               + a23 * refD100[k + 2*N_Functs];
    double p22 = a21 * refD010[k             ] 
               + a22 * refD010[k + N_Functs  ] 
               + a23 * refD010[k + 2*N_Functs];
    double p32 = a21 * refD001[k             ] 
               + a22 * refD001[k + N_Functs  ] 
               + a23 * refD001[k + 2*N_Functs];
    origD100[k + N_Functs] = i11*p12 + i12*p22 + i13*p32;
    origD010[k + N_Functs] = i21*p12 + i22*p22 + i23*p32;
    origD001[k + N_Functs] = i31*p12 + i32*p22 + i33*p32;
    
    double p13 = a31 * refD100[k             ] 
               + a32 * refD100[k + N_Functs  ] 
               + a33 * refD100[k + 2*N_Functs];
    double p23 = a31 * refD010[k             ] 
               + a32 * refD010[k + N_Functs  ] 
               + a33 * refD010[k + 2*N_Functs];
    double p33 = a31 * refD001[k             ] 
               + a32 * refD001[k + N_Functs  ] 
               + a33 * refD001[k + 2*N_Functs];
    origD100[k + 2*N_Functs] = i11*p13 + i12*p23 + i13*p33;
    origD010[k + 2*N_Functs] = i21*p13 + i22*p23 + i23*p33;
    origD001[k + 2*N_Functs] = i31*p13 + i32*p23 + i33*p33;
  }
}


