#include <cmath>

#include "QuadAffin.h"
#include "BaseCell.h"
#include "LinAlg.h"
#include "QuadFormula.h"

/** constuctor */
TQuadAffin::TQuadAffin()
{
}

/** transfer from reference element to original element */
void TQuadAffin::GetOrigFromRef(double xi, double eta, double &X, double &Y)
  const
{
  X = xc0 + xc1*xi + xc2*eta;
  Y = yc0 + yc1*xi + yc2*eta;
}

/** transfer a set of point from reference to original element */
void TQuadAffin::GetOrigFromRef(int N_Points, const double *xi,
                                const double *eta, double *X, double *Y) const
{
  double Xi, Eta;
  for(int i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;
  }
}

void TQuadAffin::GetOrigFromRef(const TQuadFormula& qf_ref,
                                TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  double absdet = std::abs(detjk);
  for(unsigned int i = 0; i < n_points; ++i)
  {
    auto p_ref = qf_ref.get_point(i);
    double weight = qf_ref.get_weight(i) * absdet;
    double x = xc0 + xc1*p_ref.x + xc2*p_ref.y;
    double y = yc0 + yc1*p_ref.x + yc2*p_ref.y;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y}});
  }
}

/** transfer from original element to reference element */
void TQuadAffin::GetRefFromOrig(double X, double Y, double &xi, double &eta)
  const
{
  double xt=(X - xc0)/detjk;
  double yt=(Y - yc0)/detjk;

  xi  =  yc2*xt - xc2*yt;
  eta = -yc1*xt + xc1*yt;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadAffin::GetOrigValues(double xi, double eta, int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref, 
                               double *uorig, double *uxorig, double *uyorig,
                               int _BaseVectDim) const
{
  if(_BaseVectDim == 1)
  {
    for(int i = 0; i < N_BaseFunct; i++) 
    {
      // D00
      uorig[i] = uref[i];
      // D10 and D01
      if (uxiref!=nullptr && uetaref!=nullptr && uxorig!=nullptr && uyorig!=nullptr)
      {
        uxorig[i]=(yc2*uxiref[i]-yc1*uetaref[i]) / detjk;
        uyorig[i]=(-xc2*uxiref[i]+xc1*uetaref[i]) / detjk;
      }
    }
  }
  else if(_BaseVectDim == 2)
  {
    // D00
    this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, uref, uorig);
    // D10 and D01
    if (uxiref!=nullptr && uetaref!=nullptr && uxorig!=nullptr && uyorig!=nullptr)
    {
      this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, uref, uxiref, uetaref,
                                uxorig, uyorig);
    }
  }
}

void TQuadAffin::GetOrigAllDerivatives(double xi, double eta, int N_BaseFunct,
                                       const double *r00, const double *r10,
                                       const double *r01, const double *r20,
                                       const double *r11, const double *r02,
                                       double *o00, double *o10, double *o01,
                                       double *o20, double *o11, double *o02,
                                       int BaseVectDim) const
{
  GetOrigValues(xi, eta, N_BaseFunct, r00, r10, r01, o00, o10, o01);
  if(!ready_for_transforming_2nd_derivatives)
  {
    double GeoData[3][3];
    // reset matrices
    std::fill((double*)&GeoData[0], (double*)&GeoData[0]+9, 0.);
    std::fill((double*)&Eye[0], (double*)&Eye[0]+9, 0.);
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    
    GeoData[0][0] = xc1*xc1;
    GeoData[0][1] = 2*xc1*yc1;
    GeoData[0][2] = yc1*yc1;
      
    GeoData[1][0] = xc1*xc2;
    GeoData[1][1] = yc1*xc2+xc1*yc2; 
    GeoData[1][2] = yc1*yc2;
          
    GeoData[2][0] = xc2*xc2;
    GeoData[2][1] = 2*xc2*yc2;
    GeoData[2][2] = yc2*yc2;
    
    // subroutine for solving a multiple systems of linear equations
    SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 3, 3, 3, 3);
    ready_for_transforming_2nd_derivatives = true;
  }
  
  for(int k=0;k<N_BaseFunct*BaseVectDim;k++)
  {
    o20[k] = Eye[0][0]*r20[k] + Eye[0][1]*r11[k] + Eye[0][2]*r02[k];
    o11[k] = Eye[1][0]*r20[k] + Eye[1][1]*r11[k] + Eye[1][2]*r02[k];
    o02[k] = Eye[2][0]*r20[k] + Eye[2][1]*r11[k] + Eye[2][2]*r02[k];
  } // endfor k
}


void TQuadAffin::GetOrigValues(int joint, double zeta, int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref,
                               double *uorig, double *uxorig, double *uyorig,
                               int _BaseVectDim) const
{
  double xi, eta;
  switch(joint)
  {
    case 0:
      xi = zeta; eta = -1;
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
    default:
      ErrThrow("wrong joint number ", joint, " in TQuadAffin::GetOrigValues.");
      break;
   }
   this->GetOrigValues(xi, eta, N_BaseFunct, uref, uxiref, uetaref, 
                       uorig, uxorig, uyorig, _BaseVectDim);
}

void TQuadAffin::SetCell(const TBaseCell *cell)
{

#ifdef __3D__
  double z0, z1, z2, z3;
#endif

  Cell = cell;

#ifdef __3D__
  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
#else
  Cell->GetVertex(0)->GetCoords(x0, y0);
  Cell->GetVertex(1)->GetCoords(x1, y1);
  Cell->GetVertex(2)->GetCoords(x2, y2);
  Cell->GetVertex(3)->GetCoords(x3, y3);
#endif

  xc0 = (x1 + x3) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;

  yc0 = (y1 + y3) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;

  detjk = xc1*yc2-xc2*yc1;
  
  ready_for_transforming_2nd_derivatives = false;
}

/** Piola transformation for vectorial basis functions */
void TQuadAffin::PiolaMapOrigFromRef(double, double, int N_Functs,
                                     const double *refD00, double *origD00)
  const
{
  double a11 = xc1 / detjk;
  double a12 = xc2 / detjk;
  double a21 = yc1 / detjk;
  double a22 = yc2 / detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    origD00[k] = a11*refD00[k] + a12*refD00[N_Functs+k];
    origD00[N_Functs+k] = a21*refD00[k] + a22*refD00[N_Functs+k];

  }
}
   
void TQuadAffin::PiolaMapOrigFromRef(double, double, int N_Functs,
                                     const double *, const double *refD10,
                                     const double *refD01, double *origD10,
                                     double *origD01) const
{
  double a11 = xc1 / detjk;
  double a12 = xc2 / detjk;
  double a21 = yc1 / detjk;
  double a22 = yc2 / detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    
    // x-component (k=0,N_Functs-1)
    double refPiolaD10_k = a11*refD10[k] + a12*refD10[N_Functs+k];
    double refPiolaD01_k = a11*refD01[k] + a12*refD01[N_Functs+k];
    origD10[k] = ( yc2*refPiolaD10_k - yc1*refPiolaD01_k) / detjk;
    origD01[k] = (-xc2*refPiolaD10_k + xc1*refPiolaD01_k) / detjk;
    
    // y-component (k=N_Functs,BaseVectDim*N_Functs-1)
    refPiolaD10_k = a21*refD10[k] + a22*refD10[N_Functs+k];
    refPiolaD01_k = a21*refD01[k] + a22*refD01[N_Functs+k];
    origD10[k+N_Functs] = ( yc2*refPiolaD10_k - yc1*refPiolaD01_k) / detjk;
    origD01[k+N_Functs] = (-xc2*refPiolaD10_k + xc1*refPiolaD01_k) / detjk;
  } // endfor k
}


