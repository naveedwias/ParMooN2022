// First order Raviart-Thomas vector element on quads, nonconforming, 2D

static double NF_N_Q_RT1_2D_a = std::sqrt(3./5.);
static double NF_N_Q_RT1_2D_q[5] =
{ -std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3., -std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3., 0.,
   std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3.,  std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3. };

static double NF_N_Q_RT1_2D_Xi[37] =
{-NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
 1, 1, 1,
 NF_N_Q_RT1_2D_a, 0, -NF_N_Q_RT1_2D_a,
 -1, -1, -1,
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4],
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4],
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4],
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4],
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4]
};
static double NF_N_Q_RT1_2D_Eta[37]  =
{-1, -1, -1,
 -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
  1, 1, 1,
  NF_N_Q_RT1_2D_a, 0, -NF_N_Q_RT1_2D_a,
  NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[0],NF_N_Q_RT1_2D_q[0],
  NF_N_Q_RT1_2D_q[0],
  NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[1],NF_N_Q_RT1_2D_q[1],
  NF_N_Q_RT1_2D_q[1],
  NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[2],NF_N_Q_RT1_2D_q[2],
  NF_N_Q_RT1_2D_q[2],
  NF_N_Q_RT1_2D_q[3],NF_N_Q_RT1_2D_q[3],NF_N_Q_RT1_2D_q[3],NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[3],
  NF_N_Q_RT1_2D_q[4],NF_N_Q_RT1_2D_q[4],NF_N_Q_RT1_2D_q[4],NF_N_Q_RT1_2D_q[4],
  NF_N_Q_RT1_2D_q[4]
};

static double NF_N_Q_RT1_2D_T[] = { -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a };

void NF_N_Q_RT1_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  for(unsigned int i = 0; i < 12; ++i)
    Functionals[i] = 0;
  
  double w[3] = { 5./9., 8./9., 5./9. }; // weights
  double iw[5] = 
   { (322. - 13.*std::sqrt(70.))/900., (322. + 13.*std::sqrt(70.))/900., 128./225.,
     (322. + 13.*std::sqrt(70.))/900., (322. - 13.*std::sqrt(70.))/900. };
  double * q = NF_N_Q_RT1_2D_T;
  double * iq = NF_N_Q_RT1_2D_q;
  
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[0] -= PointValues[37+i] * w[i];
      Functionals[1] -= PointValues[37+i] * w[i] * q[i];
      
      Functionals[2] += PointValues[3+i] * w[i];
      Functionals[3] += PointValues[3+i] * w[i] * q[i];
      
      Functionals[4] += PointValues[43+i] * w[i];
      Functionals[5] += PointValues[43+i] * w[i] * q[i];
      
      Functionals[6] -= PointValues[9+i] * w[i];
      Functionals[7] -= PointValues[9+i] * w[i] * q[i];
    }
    
    for(unsigned int i = 0; i < 25; ++i)
    {
      Functionals[8] += PointValues[12+i] * iw[i%5] * iw[i/5];
      Functionals[9] += PointValues[49+i] * iw[i%5] * iw[i/5];
      
      Functionals[10] += PointValues[12+i] * iw[i%5] * iw[i/5] * iq[i/5];
      Functionals[11] += PointValues[49+i] * iw[i%5] * iw[i/5] * iq[i%5];
    }
  }
  else // on a real cell
  {
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
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[37+i]*ny)*w[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[37+i]*ny)*w[i] * q[i];
    }
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[1] *= 0.5; // Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[2] += (PointValues[3+i]*nx + PointValues[40+i]*ny)*w[i];
      Functionals[3] += (PointValues[3+i]*nx + PointValues[40+i]*ny)*w[i]*q[i];
    }
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[3] *= 0.5; // Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[4] += (PointValues[6+i]*nx + PointValues[43+i]*ny)*w[i];
      Functionals[5] += (PointValues[6+i]*nx + PointValues[43+i]*ny)*w[i]*q[i];
    }
    Functionals[4] *= 0.5*Cell->GetNormalOrientation(2);
    Functionals[5] *= 0.5; // Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[6] += (PointValues[9+i]*nx + PointValues[46+i]*ny)*w[i];
      Functionals[7] += (PointValues[9+i]*nx + PointValues[46+i]*ny)*w[i]*q[i];
    }
    Functionals[6] *= 0.5*Cell->GetNormalOrientation(3);
    Functionals[7] *= 0.5; // Cell->GetNormalOrientation(3);
    
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
    
    for(unsigned int i = 0; i < 25; ++i) // loop over inner quadrature points
    {
      double x = NF_N_Q_RT1_2D_Xi[12+i], y = NF_N_Q_RT1_2D_Eta[12+i];
      double uxi[4] =  { 1,0,  y,0 };
      double ueta[4] = { 0,1,  0,x };
      double u[4] = {0, 0, 0, 0}; // not used
      double uorig[4], ux[4], uy[4];
      referenceTransform.GetOrigValues(x, y, 4, u, uxi, ueta, uorig, ux, uy);
      for(unsigned int d = 0; d < 4; ++d) // inner dofs
      {
        Functionals[8+d] += (  PointValues[12+i]*ux[d]
                              +PointValues[49+i]*uy[d] ) * iw[i%5] * iw[i/5];
      }
    }
    for(unsigned int d = 0; d < 4; ++d) 
    {
      Functionals[8+d] *= measure;
    }
  }
}

void NF_N_Q_RT1_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int,
                            const double *PointValues, double *Functionals)
{
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1, z); // 4=number of edges
  double l = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); // length of joint
  Functionals[0] = (5*PointValues[0]+8*PointValues[1]+5*PointValues[2])*l/18.;
  Functionals[1] = (-PointValues[0] + PointValues[2])*NF_N_Q_RT1_2D_a*l*5/18.;
}
