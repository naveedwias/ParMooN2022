// First order Brezzi-Douglas-Marini vector element, nonconforming, 2D

static double NF_N_Q_BDM1_2D_a = std::sqrt(3./5.);
static double NF_N_Q_BDM1_2D_Xi[] =
{-NF_N_Q_BDM1_2D_a, 0, NF_N_Q_BDM1_2D_a,
 1, 1, 1,
 NF_N_Q_BDM1_2D_a, 0, -NF_N_Q_BDM1_2D_a,
 -1, -1, -1
};
static double NF_N_Q_BDM1_2D_Eta[]  =
{-1, -1, -1,
 -NF_N_Q_BDM1_2D_a, 0, NF_N_Q_BDM1_2D_a,
  1, 1, 1,
  NF_N_Q_BDM1_2D_a, 0, -NF_N_Q_BDM1_2D_a
};

static double NF_N_Q_BDM1_2D_T[] = { -NF_N_Q_BDM1_2D_a, 0, NF_N_Q_BDM1_2D_a };


void NF_N_Q_BDM1_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                            const double *PointValues, double *Functionals)
{
  // short names
  const double w[3] = { 5./9., 8./9., 5./9. };
  const double * q = NF_N_Q_BDM1_2D_T;
  
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -( PointValues[12]*w[0] + PointValues[13]*w[1]
                       +PointValues[14]*w[2]);
    Functionals[1] = -( PointValues[12]*w[0]*q[0] + PointValues[13]*w[1]*q[1]
                       +PointValues[14]*w[2]*q[2]);
    
    Functionals[2] = PointValues[3]*w[0] + PointValues[4]*w[1]
                    +PointValues[5]*w[2];
    Functionals[3] = PointValues[3]*w[0]*q[0] + PointValues[4]*w[1]*q[1]
                    +PointValues[5]*w[2]*q[2];
    
    Functionals[4] = PointValues[18]*w[0] + PointValues[19]*w[1]
                    +PointValues[20]*w[2];
    Functionals[5] = PointValues[18]*w[0]*q[0] + PointValues[19]*w[1]*q[1]
                    +PointValues[20]*w[2]*q[2];
    
    Functionals[6] = -( PointValues[9]*w[0] + PointValues[10]*w[1]
                       +PointValues[11]*w[2]);
    Functionals[7] = -( PointValues[9]*w[0]*q[0] + PointValues[10]*w[1]*q[1]
                       +PointValues[11]*w[2]*q[2]);
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_BDM1_2D_EvalAll not tested for non affine ",
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
    Functionals[0] = ( PointValues[0]*w[0] + PointValues[1]*w[1]
                      +PointValues[2]*w[2] )*nx
                    +( PointValues[12]*w[0] + PointValues[13]*w[1]
                      +PointValues[14]*w[2] )*ny;
    Functionals[1] = ( PointValues[0]*w[0]*q[0] + PointValues[1]*w[1]*q[1]
                      +PointValues[2]*w[2]*q[2] )*nx
                    +( PointValues[12]*w[0]*q[0] + PointValues[13]*w[1]*q[1]
                      +PointValues[14]*w[2]*q[2] )*ny;
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[1] *= 0.5; // Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[2] = ( PointValues[3]*w[0] + PointValues[4]*w[1]
                      +PointValues[5]*w[2] )*nx
                    +( PointValues[15]*w[0] + PointValues[16]*w[1]
                      +PointValues[17]*w[2] )*ny;
    Functionals[3] = ( PointValues[3]*w[0]*q[0] + PointValues[4]*w[1]*q[1]
                      +PointValues[5]*w[2]*q[2] )*nx
                    +( PointValues[15]*w[0]*q[0] + PointValues[16]*w[1]*q[1]
                      +PointValues[17]*w[2]*q[2] )*ny;
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[3] *= 0.5; // Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[4] = ( PointValues[6]*w[0] + PointValues[7]*w[1]
                      +PointValues[8]*w[2] )*nx
                    +( PointValues[18]*w[0] + PointValues[19]*w[1]
                      +PointValues[20]*w[2] )*ny;
    Functionals[5] = ( PointValues[6]*w[0]*q[0] + PointValues[7]*w[1]*q[1]
                      +PointValues[8]*w[2]*q[2] )*nx
                    +( PointValues[18]*w[0]*q[0] + PointValues[19]*w[1]*q[1]
                      +PointValues[20]*w[2]*q[2] )*ny;
    Functionals[4] *= 0.5*Cell->GetNormalOrientation(2);
    Functionals[5] *= 0.5; // Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[6]= ( PointValues[9]*w[0] + PointValues[10]*w[1]
                      +PointValues[11]*w[2] )*nx
                    +( PointValues[21]*w[0] + PointValues[22]*w[1]
                      +PointValues[23]*w[2] )*ny;
    Functionals[7]= ( PointValues[9]*w[0]*q[0] + PointValues[10]*w[1]*q[1]
                      +PointValues[11]*w[2]*q[2] )*nx
                    +( PointValues[21]*w[0]*q[0] + PointValues[22]*w[1]*q[1]
                      +PointValues[23]*w[2]*q[2] )*ny;
    Functionals[6] *= 0.5*Cell->GetNormalOrientation(3);
    Functionals[7] *= 0.5; // Cell->GetNormalOrientation(3);
  }
}

void NF_N_Q_BDM1_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int,
                             const double *PointValues,double *Functionals)
{
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1, z); // 4=number of edges
  // length of joint, 0.5 due to 1D-reference cell having measure 2
  double l = 0.5*std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] = (5*PointValues[0]+8*PointValues[1]+5*PointValues[2])*l/9.;
  Functionals[1] = (-PointValues[0] + PointValues[2])*NF_N_Q_BDM1_2D_a*l*5/9.;
}
