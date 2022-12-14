// ***********************************************************************
// Raviart-Thomas element of zero-th order on tetrahedra, 3D
// ***********************************************************************

/* for all functionals */
static double NF_N_T_RT0_3D_Xi[]  = { 0.33333333333333333333,
                                      0.33333333333333333333,
                                      0.33333333333333333333,
                                      0 };
static double NF_N_T_RT0_3D_Eta[] = { 0.33333333333333333333,
                                      0,
                                      0.33333333333333333333,
                                      0.33333333333333333333 };
static double NF_N_T_RT0_3D_Zeta[]= { 0,
                                      0.33333333333333333333,
                                      0.33333333333333333333,
                                      0.33333333333333333333 };

/* face 0                               0 */
static double NF_N_T_RT0_3D_F0_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F0_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F0_Zeta[] = { 0 };

/* face 1                               1 */
static double NF_N_T_RT0_3D_F1_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F1_Eta[]  = { 0 };
static double NF_N_T_RT0_3D_F1_Zeta[] = { 0.33333333333333333333 };

/* face 2                               2 */
static double NF_N_T_RT0_3D_F2_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F2_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F2_Zeta[] = { 0.33333333333333333333 };

/* face 3                               3 */
static double NF_N_T_RT0_3D_F3_Xi[]   = { 0 };
static double NF_N_T_RT0_3D_F3_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_RT0_3D_F3_Zeta[] = { 0.33333333333333333333 };

static double *NF_N_T_RT0_3D_XiArray[4] = { 
                        NF_N_T_RT0_3D_F0_Xi,
                        NF_N_T_RT0_3D_F1_Xi,
                        NF_N_T_RT0_3D_F2_Xi,
                        NF_N_T_RT0_3D_F3_Xi };

static double *NF_N_T_RT0_3D_EtaArray[4] = { 
                        NF_N_T_RT0_3D_F0_Eta,
                        NF_N_T_RT0_3D_F1_Eta,
                        NF_N_T_RT0_3D_F2_Eta,
                        NF_N_T_RT0_3D_F3_Eta };

static double *NF_N_T_RT0_3D_ZetaArray[4] = { 
                        NF_N_T_RT0_3D_F0_Zeta,
                        NF_N_T_RT0_3D_F1_Zeta,
                        NF_N_T_RT0_3D_F2_Zeta,
                        NF_N_T_RT0_3D_F3_Zeta };

static double NF_N_T_RT0_3D_T[1] = {0.33333333333333333333};// ???
static double NF_N_T_RT0_3D_S[1] = {0.33333333333333333333};// ???

void NF_N_T_RT0_3D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  // PointValues[4*i + j] means i-th component (i=0 for x, i=1 for y, i=2 for z)
  // at j-th evaluation point (see NF_N_T_RT0_3D_Xi, ...Eta, ...Zeta)
  Functionals[0] = -PointValues[8];
  Functionals[1] = -PointValues[5];
  Functionals[2] =  (PointValues[2] + PointValues[6] + PointValues[10]);
  Functionals[3] = -PointValues[3];
}

void NF_N_T_RT0_3D_EvalFace(const TCollection *, const TBaseCell *Cell, int face,
                            const double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 3, length == {3,3,3,3}
  Cell->GetVertex(faceVertex[3*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[3*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[3*face + 2])->GetCoords(x2,y2,z2);
  // compute measure of this face times 2 because functional is 2 times Int_F q.n
  s = std::sqrt( std::pow((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + std::pow((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + std::pow((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  //Output::print("face: ", face, " s: ", s);
  Functionals[0] = PointValues[0]*s;
}


static int NF_N_T_RT0_3D_N_AllFunctionals = 4;
static int NF_N_T_RT0_3D_N_PointsAll = 4;
static int NF_N_T_RT0_3D_N_FaceFunctionals[] = { 1, 1, 1, 1 }; 
static int NF_N_T_RT0_3D_N_PointsFace[] = { 1, 1, 1, 1 };
