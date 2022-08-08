#include <cmath>

#include "QuadBilinear.h"
#include "BaseCell.h"
#include "LinAlg.h"
#include "QuadFormula.h"

/** constuctor */
TQuadBilinear::TQuadBilinear()
{
}

/** transfer from reference element to original element */
void TQuadBilinear::GetOrigFromRef(double xi, double eta, double &X, double &Y)
  const
{
  X = xc0 + xc1*xi + xc2*eta + xc3*xi*eta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*xi*eta;
}

/** transfer a set of point from reference to original element */
void TQuadBilinear::GetOrigFromRef(int N_Points, const double *xi,
                                   const double *eta, double *X, double *Y)
  const
{
  double Xi, Eta;
  for(int i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta +xc3*Eta*Xi;
    Y[i] = yc0 + yc1*Xi + yc2*Eta +yc3*Eta*Xi;
  }
}

void TQuadBilinear::GetOrigFromRef(const TQuadFormula& qf_ref,
                                   TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  for(unsigned int i = 0; i < n_points; ++i)
  {
    auto p_ref = qf_ref.get_point(i);
    double absdetjk = std::abs( (xc1+xc3*p_ref.y)*(yc2+yc3*p_ref.x)
                               -(xc2+xc3*p_ref.x)*(yc1+yc3*p_ref.y));
    double weight = qf_ref.get_weight(i) * absdetjk;
    double x = xc0 + xc1*p_ref.x + xc2*p_ref.y +xc3*p_ref.y*p_ref.x;
    double y = yc0 + yc1*p_ref.x + yc2*p_ref.y +yc3*p_ref.y*p_ref.x;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y}});
  }
}

/** transfer from original element to reference element */
void TQuadBilinear::GetRefFromOrig(double x, double y, double &xi, double &eta)
  const
{
  double z0, z1, xi0, eta0;
  double recdetaffine;
  double eps=1e-14;

  recdetaffine = 1/(xc1*yc2-xc2*yc1);

  z0 = ( yc2 * (x-xc0) - xc2 * (y-yc0)) * recdetaffine;
  z1 = (-yc1 * (x-xc0) + xc1 * (y-yc0)) * recdetaffine;
  
  xi0 = z0;
  eta0 = z1;

  xi = xi0+10;
  eta = eta0+10;
  while( std::abs(xi-xi0)+std::abs(eta-eta0) > eps )
  {
    xi  = xi0;
    eta = eta0;
    xi0  = z0 - ( yc2 * xc3*xi*eta - xc2 * yc3*xi*eta) * recdetaffine;
    eta0 = z1 - (-yc1 * xc3*xi*eta + xc1 * yc3*xi*eta) * recdetaffine;
  }
  xi = xi0;
  eta = eta0;

  return;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadBilinear::GetOrigValues(double xi, double eta, int N_BaseFunct,
                                  const double *uref, const double *uxiref,
                                  const double *uetaref,
                                  double *uorig, double *uxorig,
                                  double *uyorig, int _BaseVectDim) const
{
  int i;
  double rec_detjk;

  if(_BaseVectDim==1) // standard case
  {
    // D00
    for(i=0;i<N_BaseFunct;i++)
    {
      uorig[i] = uref[i];
    } // endfor i
    
    // D10 and D01
    if (uxiref!=nullptr && uetaref!=nullptr && uxorig!=nullptr && uyorig!=nullptr)
    {
      rec_detjk = 1/( (xc1+xc3*eta)*(yc2+yc3*xi) - (xc2+xc3*xi)*(yc1+yc3*eta) );
      for(i=0;i<N_BaseFunct;i++)
      {
        uxorig[i]= ( (yc2+yc3*xi) * uxiref[i] - (yc1+yc3*eta) * uetaref[i] )
                  *rec_detjk;
        uyorig[i]= (-(xc2+xc3*xi) * uxiref[i] + (xc1+xc3*eta) * uetaref[i] )
                  *rec_detjk;
      } // endfor i
    }
  }
  else
  {
    // Piola transformation
    // D00
    this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, uref, uorig);
    // D10, D01
    if (uxiref!=nullptr && uetaref!=nullptr && uxorig!=nullptr && uyorig!=nullptr)
    {
      this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, uref, uxiref, uetaref,
                                uxorig, uyorig);
    }
  }
}

void TQuadBilinear::GetOrigAllDerivatives(double xi, double eta,
                                          int N_BaseFunct,
                                          const double *r00, const double *r10,
                                          const double *r01, const double *r20,
                                          const double *r11, const double *r02,
                                          double *o00, double *o10, double *o01,
                                          double *o20, double *o11, double *o02,
                                          int BaseVectDim) const
{
  if(!ready_for_transforming_2nd_derivatives)
  {
    double GeoData[5][5]; // double GeoData1[5][5];
    // reset matrices
    std::fill((double*)&GeoData[0], (double*)&GeoData[0]+25, 0.);
    std::fill((double*)&Eye[0], (double*)&Eye[0]+25, 0.);
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    Eye[3][3] = 1;
    Eye[4][4] = 1;
    
    GeoData[0][0] = xc1+xc3*eta;
    GeoData[0][1] = yc1+yc3*eta;
    
    GeoData[1][0] = xc2+xc3*xi;
    GeoData[1][1] = yc2+yc3*xi;
    
    GeoData[2][0] = 0;
    GeoData[2][1] = 0;
    GeoData[2][2] = GeoData[0][0]*GeoData[0][0];
    GeoData[2][3] = 2*GeoData[0][0]*GeoData[0][1];
    GeoData[2][4] = GeoData[0][1]*GeoData[0][1];
    
    GeoData[3][0] = xc3;
    GeoData[3][1] = yc3;
    GeoData[3][2] = GeoData[0][0]*GeoData[1][0];
    GeoData[3][3] = GeoData[0][1]*GeoData[1][0]+GeoData[0][0]*GeoData[1][1];
    GeoData[3][4] = GeoData[0][1]*GeoData[1][1];
    
    GeoData[4][0] = 0;
    GeoData[4][1] = 0;
    GeoData[4][2] = GeoData[1][0]*GeoData[1][0];
    GeoData[4][3] = 2*GeoData[1][0]*GeoData[1][1];
    GeoData[4][4] = GeoData[1][1]*GeoData[1][1];
    
    // subroutine for solving a multiple systems of linear equations
    SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 5, 5, 5, 5);
    ready_for_transforming_2nd_derivatives = true;
  }
  
  if(BaseVectDim==1) // standard case
  {
    for(int k = 0; k < N_BaseFunct; k++)
    {
      o00[k] = r00[k];
      o10[k] = Eye[0][0]*r10[k] + Eye[0][1]*r01[k] + Eye[0][2]*r20[k]
              +Eye[0][3]*r11[k] + Eye[0][4]*r02[k];
      o01[k] = Eye[1][0]*r10[k] + Eye[1][1]*r01[k] + Eye[1][2]*r20[k]
              +Eye[1][3]*r11[k] + Eye[1][4]*r02[k];
      o20[k] = Eye[2][0]*r10[k] + Eye[2][1]*r01[k] + Eye[2][2]*r20[k]
              +Eye[2][3]*r11[k] + Eye[2][4]*r02[k];
      o11[k] = Eye[3][0]*r10[k] + Eye[3][1]*r01[k] + Eye[3][2]*r20[k]
              +Eye[3][3]*r11[k] + Eye[3][4]*r02[k];
      o02[k] = Eye[4][0]*r10[k] + Eye[4][1]*r01[k] + Eye[4][2]*r20[k]
              +Eye[4][3]*r11[k] + Eye[4][4]*r02[k];
    } // endfor k
  }
  else
  {
    this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, r00, o00);
    // D10, D01
    if(r10!=nullptr && r01!=nullptr && o10!=nullptr && o01!=nullptr)
    {
      this->PiolaMapOrigFromRef(xi, eta, N_BaseFunct, r00, r10, r01, o10, o01);
    }
    ErrThrow("Piola transformation is not implemented for second order "
             "derivatives");
  }
}


void TQuadBilinear::GetOrigValues(int joint, double zeta, int N_BaseFunct,
                                  const double *uref, const double *uxiref,
                                  const double *uetaref,
                                  double *uorig, double *uxorig,
                                  double *uyorig,
                                  int _BaseVectDim) const
{
  double xi=0, eta=0;
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
  this->GetOrigValues(xi, eta, N_BaseFunct, uref, uxiref, uetaref, 
                      uorig, uxorig, uyorig, _BaseVectDim);
}


void TQuadBilinear::SetCell(const TBaseCell *cell)
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

  xc0=( x0 + x1 + x2 + x3) * 0.25;
  xc1=(-x0 + x1 + x2 - x3) * 0.25;
  xc2=(-x0 - x1 + x2 + x3) * 0.25;
  xc3=( x0 - x1 + x2 - x3) * 0.25;

  yc0=( y0 + y1 + y2 + y3) * 0.25;
  yc1=(-y0 + y1 + y2 - y3) * 0.25;
  yc2=(-y0 - y1 + y2 + y3) * 0.25;
  yc3=( y0 - y1 + y2 - y3) * 0.25;
  ready_for_transforming_2nd_derivatives = false;
}

void TQuadBilinear::PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                                        const double *refD00, double *origD00)
  const
{
  double rec_detjk = 1./((xc1+xc3*eta)*(yc2+yc3*xi)-(xc2+xc3*xi)*(yc1+yc3*eta));
  double a11 = (xc1+xc3*eta)*rec_detjk;
  double a12 = (xc2+xc3*xi )*rec_detjk;
  double a21 = (yc1+yc3*eta)*rec_detjk;
  double a22 = (yc2+yc3*xi )*rec_detjk;
  
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    origD00[k] = a11*refD00[k]+a12*refD00[N_Functs+k];
    origD00[N_Functs+k] = a21*refD00[k]+a22*refD00[N_Functs+k];
  }
}

void TQuadBilinear::PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                                        const double *refD00,
                                        const double *refD10,
                                        const double *refD01,
                                        double *origD10, double *origD01) const
{
  // Piola transformation
  // phi = 1/|J| DF phi_hat
  // 
  // this is \hat{D}F
  double a11 = (xc1+xc3*eta);
  double a12 = (xc2+xc3* xi);
  double a21 = (yc1+yc3*eta);
  double a22 = (yc2+yc3* xi);
  double rec_detjk = 1./(a11*a22 - a12*a21);
  // this is DF^{-1}
  double z11 =  a22*rec_detjk;
  double z12 = -a12*rec_detjk;
  double z21 = -a21*rec_detjk;
  double z22 =  a11*rec_detjk;
  // D(det(DF)) (row vector)
  double b1 = xc1*yc3 - yc1*xc3;
  double b2 = yc2*xc3 - xc2*yc3;
  b1 *= -sgn(rec_detjk)*rec_detjk*rec_detjk;
  b2 *= -sgn(rec_detjk)*rec_detjk*rec_detjk;
  
  // the derivative has three parts (since the Piola transform is a 
  // product of three factors, namely 1/|J|, DF, and phi_hat) 
  // according to the product rule
  
  // refPiola_(xy)_D.. is the derivative with respect to the 
  // reference variables xi and eta. In the end we apply the chain 
  // rule to get the derivative with respect to x and y.
  double refPiola_x_D10,refPiola_x_D01,refPiola_y_D10,refPiola_y_D01;
  
  for(int k = 0; k < N_Functs; k++)
  {
    // first part (differentiate 1/|J|)
    // \hat{D}F \hat{v}
    double DFv1 = a11*refD00[k] + a12*refD00[k+N_Functs];
    double DFv2 = a21*refD00[k] + a22*refD00[k+N_Functs];
    
    // Dfv (-sign(J)/(J^2))DJ
    refPiola_x_D10 = DFv1*b1;
    refPiola_x_D01 = DFv1*b2;
    refPiola_y_D10 = DFv2*b1;
    refPiola_y_D01 = DFv2*b2;
    
    // second part (differentiate DF)
    // D^2F v /det(DF)
    refPiola_x_D10 += refD00[k+N_Functs]*xc3*rec_detjk;
    refPiola_x_D01 += refD00[k]         *xc3*rec_detjk;
    refPiola_y_D10 += refD00[k+N_Functs]*yc3*rec_detjk;
    refPiola_y_D01 += refD00[k]         *yc3*rec_detjk;
    
    // third part (differentiate phi_hat), similar to the affine linear 
    // case
    // \hat{D}F \hat{D}\hat{v} /|J|
    refPiola_x_D10 += (a11*refD10[k]+a12*refD10[N_Functs+k])*rec_detjk;
    refPiola_x_D01 += (a11*refD01[k]+a12*refD01[N_Functs+k])*rec_detjk;
    refPiola_y_D10 += (a21*refD10[k]+a22*refD10[N_Functs+k])*rec_detjk;
    refPiola_y_D01 += (a21*refD01[k]+a22*refD01[N_Functs+k])*rec_detjk;
    
    // apply DF^{-1} from right (chain rule)
    origD10[k] = z11*refPiola_x_D10 + z21*refPiola_x_D01;
    origD01[k] = z12*refPiola_x_D10 + z22*refPiola_x_D01;
    origD10[k+N_Functs] = z11*refPiola_y_D10 + z21*refPiola_y_D01;
    origD01[k+N_Functs] = z12*refPiola_y_D10 + z22*refPiola_y_D01;
    
  } // endfor k
}

