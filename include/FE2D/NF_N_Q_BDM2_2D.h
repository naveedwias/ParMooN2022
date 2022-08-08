// Second order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// points for 1D Gauss quadrature with four points (two are symmetric)
// { -0.86114, -0.33998, 0.33998, 0.86114 }
static double NF_N_Q_BDM2_2D_q[4] = 
 {-std::sqrt(3./7. + (2./7.)*std::sqrt(6./5.)), -std::sqrt(3./7. - (2./7.)*std::sqrt(6./5.)),
   std::sqrt(3./7. - (2./7.)*std::sqrt(6./5.)),  std::sqrt(3./7. + (2./7.)*std::sqrt(6./5.)) };
// weights for the 1D Gauss quadrature with four points
// { 0.34785, 0.65215, 0.65215, 0.34785 }
static double NF_N_Q_BDM2_2D_w[4] =
{ (18. - std::sqrt(30.)) / 36., (18. + std::sqrt(30.)) / 36.,
  (18. + std::sqrt(30.)) / 36., (18. - std::sqrt(30.)) / 36. }; 

// P_2(x) with P_2 being the second Legendre polynomial and x the quad points 
// from above
static double NF_N_Q_BDM2_2D_p2[4] = 
  { 0.5*(3*NF_N_Q_BDM2_2D_q[0]*NF_N_Q_BDM2_2D_q[0] - 1.),
    0.5*(3*NF_N_Q_BDM2_2D_q[1]*NF_N_Q_BDM2_2D_q[1] - 1.),
    0.5*(3*NF_N_Q_BDM2_2D_q[2]*NF_N_Q_BDM2_2D_q[2] - 1.),
    0.5*(3*NF_N_Q_BDM2_2D_q[3]*NF_N_Q_BDM2_2D_q[3] - 1.) };

static double NF_N_Q_BDM2_2D_Xi[] = 
{ NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
   1, 1, 1, 1,
  NF_N_Q_BDM2_2D_q[3], NF_N_Q_BDM2_2D_q[2] , NF_N_Q_BDM2_2D_q[1], 
  NF_N_Q_BDM2_2D_q[0],
  -1, -1, -1, -1,
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
};
static double NF_N_Q_BDM2_2D_Eta[] = 
{-1, -1, -1, -1,
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[3],
  1, 1, 1, 1,
  NF_N_Q_BDM2_2D_q[3], NF_N_Q_BDM2_2D_q[2] , NF_N_Q_BDM2_2D_q[1], 
  NF_N_Q_BDM2_2D_q[0],
  NF_N_Q_BDM2_2D_q[0], NF_N_Q_BDM2_2D_q[0] , NF_N_Q_BDM2_2D_q[0], 
  NF_N_Q_BDM2_2D_q[0],
  NF_N_Q_BDM2_2D_q[1], NF_N_Q_BDM2_2D_q[1] , NF_N_Q_BDM2_2D_q[1], 
  NF_N_Q_BDM2_2D_q[1],
  NF_N_Q_BDM2_2D_q[2], NF_N_Q_BDM2_2D_q[2] , NF_N_Q_BDM2_2D_q[2], 
  NF_N_Q_BDM2_2D_q[2],
  NF_N_Q_BDM2_2D_q[3], NF_N_Q_BDM2_2D_q[3] , NF_N_Q_BDM2_2D_q[3], 
  NF_N_Q_BDM2_2D_q[3],
};



void NF_N_Q_BDM2_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                            const double *PointValues, double *Functionals)
{
  // short names
  const double * q = NF_N_Q_BDM2_2D_q;
  const double * w = NF_N_Q_BDM2_2D_w;
  const double * p2 = NF_N_Q_BDM2_2D_p2;
  
  // set all functionals to zero first
  for(unsigned int i = 0; i < 14; ++i)
    Functionals[i] = 0;
  
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[0] -= PointValues[32+i]*w[i];
      Functionals[1] -= PointValues[32+i]*w[i] * q[i];
      Functionals[2] -= PointValues[32+i]*w[i] * p2[i];
      
      Functionals[3] += PointValues[4+i]*w[i];
      Functionals[4] += PointValues[4+i]*w[i] * q[i];
      Functionals[5] += PointValues[4+i]*w[i] * p2[i];
      
      Functionals[6] += PointValues[40+i]*w[i];
      Functionals[7] += PointValues[40+i]*w[i] * q[i];
      Functionals[8] += PointValues[40+i]*w[i] * p2[i];
      
      Functionals[9] -= PointValues[12+i]*w[i];
      Functionals[10]-= PointValues[12+i]*w[i] * q[i];
      Functionals[11]-= PointValues[12+i]*w[i] * p2[i];
      
      for(unsigned int j =0; j < 4; ++j)
      {
        Functionals[12] += PointValues[16+4*i+j] * w[i]*w[j];
        Functionals[13] += PointValues[48+4*i+j] * w[i]*w[j];
      }
    }
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_BDM2_2D_EvalAll not tested for non affine ",
               "reference transformations");
    }
    double x0, x1, x2, x3, y0, y1, y2, y3, z; // z remains zero in 2D
    Cell->GetVertex(0)->GetCoords(x0, y0, z);
    Cell->GetVertex(1)->GetCoords(x1, y1, z);
    Cell->GetVertex(2)->GetCoords(x2, y2, z);
    Cell->GetVertex(3)->GetCoords(x3, y3, z);
    
    // outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[32+i]*ny)*w[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[32+i]*ny)*w[i] * q[i];
      Functionals[2] += (PointValues[i]*nx + PointValues[32+i]*ny)*w[i] * p2[i];
    }
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[1] *= 0.5; // Cell->GetNormalOrientation(0);
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 4; ++i)
    {
      double pvx = PointValues[4+i], pvy = PointValues[36+i];
      Functionals[3] += (pvx*nx + pvy*ny)*w[i];
      Functionals[4] += (pvx*nx + pvy*ny)*w[i] * q[i];
      Functionals[5] += (pvx*nx + pvy*ny)*w[i] * p2[i];
    }
    Functionals[3] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[4] *= 0.5; // Cell->GetNormalOrientation(1);
    Functionals[5] *= 0.5*Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    for(unsigned int i = 0; i < 4; ++i)
    {
      double pvx = PointValues[8+i], pvy = PointValues[40+i];
      Functionals[6] += (pvx*nx + pvy*ny)*w[i];
      Functionals[7] += (pvx*nx + pvy*ny)*w[i] * q[i];
      Functionals[8] += (pvx*nx + pvy*ny)*w[i] * p2[i];
    }
    Functionals[6] *= 0.5*Cell->GetNormalOrientation(2);
    Functionals[7] *= 0.5; // Cell->GetNormalOrientation(2);
    Functionals[8] *= 0.5*Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    for(unsigned int i = 0; i < 4; ++i)
    {
      double pvx = PointValues[12+i], pvy = PointValues[44+i];
      Functionals[9] += (pvx*nx + pvy*ny)*w[i];
      Functionals[10]+= (pvx*nx + pvy*ny)*w[i] * q[i];
      Functionals[11]+= (pvx*nx + pvy*ny)*w[i] * p2[i];
    }
    Functionals[9] *= 0.5*Cell->GetNormalOrientation(3);
    Functionals[10]*= 0.5; // Cell->GetNormalOrientation(3);
    Functionals[11]*= 0.5*Cell->GetNormalOrientation(3);
    
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 0.25*Cell->GetMeasure();
    
    // we could use TQuadAffin if the Cell is a parallelogram, but this also
    // works for general quads
    TQuadBilinear referenceTransform;
    referenceTransform.SetCell(Cell);
    
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    // more short names 
    double * xi = NF_N_Q_BDM2_2D_Xi;
    double *eta = NF_N_Q_BDM2_2D_Eta;
    
    // int_cell v . (1 0)^T
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i],  eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[12] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i%4] * w[i/4];
    }
    Functionals[12] *= measure;
    
    // int_cell v . (0 1)^T
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = 1., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[13] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i%4] * w[i/4];
    }
    Functionals[13] *= measure;
  }
}

void NF_N_Q_BDM2_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int,
                             const double *PointValues, double *Functionals)
{
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1, z); // 4=number of edges
  // length of joint, 0.5 due to 1D-reference cell having measure 2
  double l = 0.5*std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); 
  Functionals[0] = ( NF_N_Q_BDM2_2D_w[0]*PointValues[0]
                    +NF_N_Q_BDM2_2D_w[1]*PointValues[1]
                    +NF_N_Q_BDM2_2D_w[2]*PointValues[2]
                    +NF_N_Q_BDM2_2D_w[3]*PointValues[3] )*l;
  Functionals[1] = ( NF_N_Q_BDM2_2D_w[0]*NF_N_Q_BDM2_2D_q[0]*PointValues[0]
                    +NF_N_Q_BDM2_2D_w[1]*NF_N_Q_BDM2_2D_q[1]*PointValues[1]
                    +NF_N_Q_BDM2_2D_w[2]*NF_N_Q_BDM2_2D_q[2]*PointValues[2]
                    +NF_N_Q_BDM2_2D_w[3]*NF_N_Q_BDM2_2D_q[3]*PointValues[3] )*l;
  Functionals[2] = ( NF_N_Q_BDM2_2D_w[0]*NF_N_Q_BDM2_2D_p2[0]*PointValues[0]
                    +NF_N_Q_BDM2_2D_w[1]*NF_N_Q_BDM2_2D_p2[1]*PointValues[1]
                    +NF_N_Q_BDM2_2D_w[2]*NF_N_Q_BDM2_2D_p2[2]*PointValues[2]
                    +NF_N_Q_BDM2_2D_w[3]*NF_N_Q_BDM2_2D_p2[3]*PointValues[3])*l;
}
