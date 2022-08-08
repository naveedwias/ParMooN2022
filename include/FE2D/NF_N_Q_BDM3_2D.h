// Third order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// points for 1D Gauss quadrature with five points (symmetric)
static double NF_N_Q_BDM3_2D_q[5] =
{ -std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3., -std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3., 0.,
   std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3.,  std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3. };
// weights for the 1D Gauss quadrature with five points
static double NF_N_Q_BDM3_2D_w[5] = 
{ (322. - 13.*std::sqrt(70.))/900., (322. + 13.*std::sqrt(70.))/900., 128./225.,
  (322. + 13.*std::sqrt(70.))/900., (322. - 13.*std::sqrt(70.))/900. };
// P_2(x) with P_2 being the second Legendre polynomial and x the Gauss 
// quadrature // points from above
static double NF_N_Q_BDM3_2D_p2[5] = 
 { 0.5*(3*NF_N_Q_BDM3_2D_q[0]*NF_N_Q_BDM3_2D_q[0] - 1.), 
   0.5*(3*NF_N_Q_BDM3_2D_q[1]*NF_N_Q_BDM3_2D_q[1] - 1.),
   0.5*(3*NF_N_Q_BDM3_2D_q[2]*NF_N_Q_BDM3_2D_q[2] - 1.),
   0.5*(3*NF_N_Q_BDM3_2D_q[3]*NF_N_Q_BDM3_2D_q[3] - 1.),
   0.5*(3*NF_N_Q_BDM3_2D_q[4]*NF_N_Q_BDM3_2D_q[4] - 1.) };
// P_3(x) with P_3 being the third Legendre polynomial and x the Gauss 
// quadrature // points from above
static double NF_N_Q_BDM3_2D_p3[5] = 
{ 0.5*( 5.*NF_N_Q_BDM3_2D_q[0]*NF_N_Q_BDM3_2D_q[0]*NF_N_Q_BDM3_2D_q[0] 
       -3.*NF_N_Q_BDM3_2D_q[0]),
  0.5*( 5.*NF_N_Q_BDM3_2D_q[1]*NF_N_Q_BDM3_2D_q[1]*NF_N_Q_BDM3_2D_q[1] 
       -3.*NF_N_Q_BDM3_2D_q[1]),
  0.5*( 5.*NF_N_Q_BDM3_2D_q[2]*NF_N_Q_BDM3_2D_q[2]*NF_N_Q_BDM3_2D_q[2] 
       -3.*NF_N_Q_BDM3_2D_q[2]),
  0.5*( 5.*NF_N_Q_BDM3_2D_q[3]*NF_N_Q_BDM3_2D_q[3]*NF_N_Q_BDM3_2D_q[3] 
       -3.*NF_N_Q_BDM3_2D_q[3]),
  0.5*( 5.*NF_N_Q_BDM3_2D_q[4]*NF_N_Q_BDM3_2D_q[4]*NF_N_Q_BDM3_2D_q[4] 
       -3.*NF_N_Q_BDM3_2D_q[4]),
};

static double NF_N_Q_BDM3_2D_Xi[] = 
{ 
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  1, 1, 1, 1, 1,
  
NF_N_Q_BDM3_2D_q[4],NF_N_Q_BDM3_2D_q[3],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[1],
  NF_N_Q_BDM3_2D_q[0],
  -1, -1, -1, -1, -1,
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
};
static double NF_N_Q_BDM3_2D_Eta[] = 
{ -1, -1, -1, -1, -1,
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[4],
  1, 1, 1, 1, 1,
  
NF_N_Q_BDM3_2D_q[4],NF_N_Q_BDM3_2D_q[3],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[1],
  NF_N_Q_BDM3_2D_q[0],
  
NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[0],NF_N_Q_BDM3_2D_q[0],
  NF_N_Q_BDM3_2D_q[0],
  
NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[1],NF_N_Q_BDM3_2D_q[1],
  NF_N_Q_BDM3_2D_q[1],
  
NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[2],NF_N_Q_BDM3_2D_q[2],
  NF_N_Q_BDM3_2D_q[2],
  
NF_N_Q_BDM3_2D_q[3],NF_N_Q_BDM3_2D_q[3],NF_N_Q_BDM3_2D_q[3],NF_N_Q_BDM3_2D_q[3],
  NF_N_Q_BDM3_2D_q[3],
  
NF_N_Q_BDM3_2D_q[4],NF_N_Q_BDM3_2D_q[4],NF_N_Q_BDM3_2D_q[4],NF_N_Q_BDM3_2D_q[4],
  NF_N_Q_BDM3_2D_q[4],
};


void NF_N_Q_BDM3_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                            const double *PointValues, double *Functionals)
{
  // short names
  const double * q = NF_N_Q_BDM3_2D_q;
  const double * w = NF_N_Q_BDM3_2D_w;
  const double * p2 = NF_N_Q_BDM3_2D_p2;
  const double * p3 = NF_N_Q_BDM3_2D_p3;
  
  // set all Functionals to zero intitially
  for(unsigned int i = 0; i < 22; ++i)
    Functionals[i] = 0.;
  
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    for(unsigned int i = 0; i < 5; ++i)
    {
      // first edge
      Functionals[0] -= PointValues[45+i] * w[i];
      Functionals[1] -= PointValues[45+i] * w[i] * q[i];
      Functionals[2] -= PointValues[45+i] * w[i] * p2[i];
      Functionals[3] -= PointValues[45+i] * w[i] * p3[i];
      // second edge
      Functionals[4] += PointValues[5+i] * w[i];
      Functionals[5] += PointValues[5+i] * w[i] * q[i];
      Functionals[6] += PointValues[5+i] * w[i] * p2[i];
      Functionals[7] += PointValues[5+i] * w[i] * p3[i];
      // third edge
      Functionals[8]  += PointValues[55+i] * w[i];
      Functionals[9]  += PointValues[55+i] * w[i] * q[i];
      Functionals[10] += PointValues[55+i] * w[i] * p2[i];
      Functionals[11] += PointValues[55+i] * w[i] * p3[i];
      // fourth edge
      Functionals[12] -= PointValues[15+i] * w[i];
      Functionals[13] -= PointValues[15+i] * w[i] * q[i];
      Functionals[14] -= PointValues[15+i] * w[i] * p2[i];
      Functionals[15] -= PointValues[15+i] * w[i] * p3[i];
      
      // inner dofs
      double pvx, pvy; // PointValues, x and y component
      for(unsigned int j = 0; j < 5; ++j)
      {
        pvx = PointValues[20+i*5+j];
        pvy = PointValues[65+i*5+j];
        // test with (1 0)^T and (0 1)^T
        Functionals[16] += pvx * w[i]*w[j];
        Functionals[17] += pvy * w[i]*w[j];
        // test with (x 0)^T and (0 x)^T
        Functionals[18] += pvx * w[i]*w[j] * q[j];
        Functionals[19] += pvy * w[i]*w[j] * q[j];
        // test with (y 0)^T and (0 y)^T
        Functionals[20] += pvx * w[i]*w[j] * q[i];
        Functionals[21] += pvy * w[i]*w[j] * q[i];
      }
    }
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_BDM3_2D_EvalAll not tested for non affine ",
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
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[45+i]*ny)*w[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[45+i]*ny)*w[i] * q[i];
      Functionals[2] += (PointValues[i]*nx + PointValues[45+i]*ny)*w[i] * p2[i];
      Functionals[3] += (PointValues[i]*nx + PointValues[45+i]*ny)*w[i] * p3[i];
    }
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[1] *= 0.5; // Cell->GetNormalOrientation(0);
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[3] *= 0.5; // Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[4] += (PointValues[5+i]*nx + PointValues[50+i]*ny)*w[i];
      Functionals[5] += (PointValues[5+i]*nx + PointValues[50+i]*ny)*w[i]*q[i];
      Functionals[6] += (PointValues[5+i]*nx + PointValues[50+i]*ny)*w[i]*p2[i];
      Functionals[7] += (PointValues[5+i]*nx + PointValues[50+i]*ny)*w[i]*p3[i];
    }
    Functionals[4] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[5] *= 0.5; // Cell->GetNormalOrientation(1);
    Functionals[6] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[7] *= 0.5; // Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[8] += (PointValues[10+i]*nx+ PointValues[55+i]*ny)*w[i];
      Functionals[9] += (PointValues[10+i]*nx+ PointValues[55+i]*ny)*w[i]*q[i];
      Functionals[10]+= (PointValues[10+i]*nx+ PointValues[55+i]*ny)*w[i]*p2[i];
      Functionals[11]+= (PointValues[10+i]*nx+ PointValues[55+i]*ny)*w[i]*p3[i];
    }
    Functionals[8] *= 0.5*Cell->GetNormalOrientation(2);
    Functionals[9] *= 0.5; // Cell->GetNormalOrientation(2);
    Functionals[10]*= 0.5*Cell->GetNormalOrientation(2);
    Functionals[11]*= 0.5; // Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[12]+= (PointValues[15+i]*nx+ PointValues[60+i]*ny)*w[i];
      Functionals[13]+= (PointValues[15+i]*nx+ PointValues[60+i]*ny)*w[i]*q[i];
      Functionals[14]+= (PointValues[15+i]*nx+ PointValues[60+i]*ny)*w[i]*p2[i];
      Functionals[15]+= (PointValues[15+i]*nx+ PointValues[60+i]*ny)*w[i]*p3[i];
    }
    Functionals[12]*= 0.5*Cell->GetNormalOrientation(3);
    Functionals[13]*= 0.5; // Cell->GetNormalOrientation(3);
    Functionals[14]*= 0.5*Cell->GetNormalOrientation(3);
    Functionals[15]*= 0.5; // Cell->GetNormalOrientation(3);
    
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
    
    // the dofs:
    // int_cell v . (1 0)^T 
    // int_cell v . (0 1)^T 
    // int_cell v . (x 0)^T
    // int_cell v . (0 x)^T
    // int_cell v . (y 0)^T
    // int_cell v . (0 y)^T
    for(unsigned int i = 0; i < 25; ++i) // loop over inner quadrature points
    {
      double x = NF_N_Q_BDM3_2D_Xi[20+i], y = NF_N_Q_BDM3_2D_Eta[20+i];
      double uxi[6] = {  1,0,  x,0,  y,0 };
      double ueta[6] = { 0,1,  0,x,  0,y };
      for(unsigned int d = 0; d < 6; ++d) // inner dofs
      {
        double uref = 0., uxiref = uxi[d], uetaref = ueta[d];
        double uorig, uxorig, uyorig;
        referenceTransform.GetOrigValues(x, y, 1, &uref, &uxiref, &uetaref, 
                                         &uorig, &uxorig, &uyorig);
        Functionals[16+d] += ( PointValues[20+i]*uxorig
                              +PointValues[65+i]*uyorig ) * w[i%5] * w[i/5];
      }
    }
    for(unsigned int d = 0; d < 6; ++d) 
    {
      Functionals[16+d] *= measure;
    }
  }
}

void NF_N_Q_BDM3_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int,
                             const double *PointValues, double *Functionals)
{
  Functionals[0] = 0.;
  Functionals[1] = 0.;
  Functionals[2] = 0.;
  Functionals[3] = 0.;
  for(unsigned int i = 0; i < 5; ++i)
  {
    double pv = PointValues[i];
    Functionals[0] += pv * NF_N_Q_BDM3_2D_w[i];
    Functionals[1] += pv * NF_N_Q_BDM3_2D_w[i] * NF_N_Q_BDM3_2D_q[i];
    Functionals[2] += pv * NF_N_Q_BDM3_2D_w[i] * NF_N_Q_BDM3_2D_p2[i];
    Functionals[3] += pv * NF_N_Q_BDM3_2D_w[i] * NF_N_Q_BDM3_2D_p3[i];
  }
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1, z);// 4=number of edges
  // length of joint, 0.5 due to 1D-reference cell having measure 2
  double l = 0.5*std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] *= l;
  Functionals[1] *= l;
  Functionals[2] *= l;
  Functionals[3] *= l;
}
