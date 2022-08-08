// ***********************************************************************
// Raviart-Thomas element of zero-th order on hexahedra, 3D
// ***********************************************************************

// we use a tensor product of a 2-point Gauss quadrature on each face (which is 
// exact for polynomials up to order 3)
const double s13 = std::sqrt(1./3.);
/* for all functionals, these should use (higher order) quadrature formulas as 
well. We need to find a much easier way to implement this. */
static double NF_N_H_RT0_3D_Xi[] =
 { -s13, s13, -s13, s13, -s13, s13, -s13, s13, 1, 1, 1, 1, -s13, s13, -s13, s13,
   -1, -1, -1, -1, -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_Eta[] = 
  { -s13, -s13, s13, s13, -1, -1, -1, -1, -s13, s13, -s13, s13, 1, 1, 1, 1,
    -s13, s13, -s13, s13, -s13, -s13, s13, s13};
static double NF_N_H_RT0_3D_Zeta[] = 
  { -1, -1, -1, -1, -s13, -s13, s13, s13, -s13, -s13, s13, s13,
    -s13, -s13, s13, s13, -s13, -s13, s13, s13, 1, 1, 1, 1 };

/* face 0                               0 */
static double NF_N_H_RT0_3D_F0_Xi[]   = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F0_Eta[]  = { -s13, -s13, s13, s13 };
static double NF_N_H_RT0_3D_F0_Zeta[] = { -1, -1, -1, -1 };

/* face 1                               1 */
static double NF_N_H_RT0_3D_F1_Xi[]   = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F1_Eta[]  = { -1, -1, -1, -1 };
static double NF_N_H_RT0_3D_F1_Zeta[] = { -s13, -s13, s13, s13 };

/* face 2                               2 */
static double NF_N_H_RT0_3D_F2_Xi[]   = { 1, 1, 1, 1 };
static double NF_N_H_RT0_3D_F2_Eta[]  = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F2_Zeta[] = { -s13, -s13, s13, s13 };

/* face 3                               3 */
static double NF_N_H_RT0_3D_F3_Xi[]   = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F3_Eta[]  = { 1, 1, 1, 1 };
static double NF_N_H_RT0_3D_F3_Zeta[] = { -s13, -s13, s13, s13 };

/* face 4                               4 */
static double NF_N_H_RT0_3D_F4_Xi[]   = { -1, -1, -1, -1 };
static double NF_N_H_RT0_3D_F4_Eta[]  = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F4_Zeta[] = { -s13, -s13, s13, s13 };

/* face 5                               5 */
static double NF_N_H_RT0_3D_F5_Xi[]   = { -s13, s13, -s13, s13 };
static double NF_N_H_RT0_3D_F5_Eta[]  = { -s13, -s13, s13, s13 };
static double NF_N_H_RT0_3D_F5_Zeta[] = { 1, 1, 1, 1 };

static double *NF_N_H_RT0_3D_XiArray[6] = { 
                        NF_N_H_RT0_3D_F0_Xi,
                        NF_N_H_RT0_3D_F1_Xi,
                        NF_N_H_RT0_3D_F2_Xi,
                        NF_N_H_RT0_3D_F3_Xi,
                        NF_N_H_RT0_3D_F4_Xi,
                        NF_N_H_RT0_3D_F5_Xi };

static double *NF_N_H_RT0_3D_EtaArray[6] = { 
                        NF_N_H_RT0_3D_F0_Eta,
                        NF_N_H_RT0_3D_F1_Eta,
                        NF_N_H_RT0_3D_F2_Eta,
                        NF_N_H_RT0_3D_F3_Eta,
                        NF_N_H_RT0_3D_F4_Eta,
                        NF_N_H_RT0_3D_F5_Eta };

static double *NF_N_H_RT0_3D_ZetaArray[6] = { 
                        NF_N_H_RT0_3D_F0_Zeta,
                        NF_N_H_RT0_3D_F1_Zeta,
                        NF_N_H_RT0_3D_F2_Zeta,
                        NF_N_H_RT0_3D_F3_Zeta,
                        NF_N_H_RT0_3D_F4_Zeta,
                        NF_N_H_RT0_3D_F5_Zeta };

static double NF_N_H_RT0_3D_T[1] = {0};// ???
static double NF_N_H_RT0_3D_S[1] = {0};// ???

// forward declaration
void NF_N_H_RT0_3D_EvalFace(const TCollection *Coll, const TBaseCell *Cell, int face,
                            const double *PointValues, double *Functionals);

void NF_N_H_RT0_3D_EvalAll(const TCollection *Coll, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  if(Coll != nullptr && Cell != nullptr)
  {
    for(unsigned int i = 0; i < 6; ++i)
      Functionals[i] = 0;
    
//     double x0, x1, x2, x3, x4, x5;
//     double y0, y1, y2, y3, y4, y5;
//     double z0, z1, z2, z3, z4, z5;
//     Cell->GetVertex(0)->GetCoords(x0, y0, z0);
//     Cell->GetVertex(1)->GetCoords(x1, y1, z1);
//     Cell->GetVertex(2)->GetCoords(x2, y2, z2);
//     Cell->GetVertex(3)->GetCoords(x3, y3, z3);
//     Cell->GetVertex(4)->GetCoords(x4, y4, z4);
//     Cell->GetVertex(5)->GetCoords(x5, y5, z5);
    ErrThrow("NF_N_H_RT0_3D_EvalAll not implemented");
  }
  else
  {
    Functionals[0] = -(PointValues[48] + PointValues[49]
                     + PointValues[50] + PointValues[51]);
    Functionals[1] = -(PointValues[28] + PointValues[29]
                     + PointValues[30] + PointValues[31]);
    Functionals[2] = PointValues[8] + PointValues[9]
                     + PointValues[10] + PointValues[11];
    Functionals[3] = PointValues[36] + PointValues[37]
                     + PointValues[38] + PointValues[39];
    Functionals[4] = -(PointValues[16] + PointValues[17]
                     + PointValues[18] + PointValues[19]);
    Functionals[3] = PointValues[68] + PointValues[69]
                     + PointValues[70] + PointValues[71];
  }
}

void NF_N_H_RT0_3D_EvalFace(const TCollection *, const TBaseCell *Cell, int face,
                            const double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 4, length == {4,4,4,4}
  Cell->GetVertex(faceVertex[4*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[4*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[4*face + 2])->GetCoords(x2,y2,z2);
  // compute measure of this face
  s = std::sqrt( std::pow((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + std::pow((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + std::pow((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = (PointValues[0] + PointValues[1] + PointValues[2] 
                    + PointValues[3])*s/4.;
}

static int NF_N_H_RT0_3D_N_AllFunctionals = 6;
static int NF_N_H_RT0_3D_N_PointsAll = 24;
static int NF_N_H_RT0_3D_N_FaceFunctionals[] = { 1, 1, 1, 1, 1, 1 };
static int NF_N_H_RT0_3D_N_PointsFace[] = { 4, 4, 4, 4, 4, 4 };
