// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order on hexahedra, 3D
// ***********************************************************************

//tschebyscheff point
static double BDDF1_tp = 0.707106781186547;

static double NF_N_H_BDDF1_3D_Xi[] = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp, -BDDF1_tp, -BDDF1_tp, BDDF1_tp, 1,1,1,
    -BDDF1_tp, -BDDF1_tp, BDDF1_tp, -1,-1,-1, -BDDF1_tp, -BDDF1_tp, BDDF1_tp};
static double NF_N_H_BDDF1_3D_Eta[]  = { -BDDF1_tp, -BDDF1_tp, BDDF1_tp, -1,-1,-1, -BDDF1_tp, -BDDF1_tp, BDDF1_tp,
    1,1,1, -BDDF1_tp, BDDF1_tp, -BDDF1_tp, -BDDF1_tp, BDDF1_tp, -BDDF1_tp};
static double NF_N_H_BDDF1_3D_Zeta[] = {-1,-1,-1, -BDDF1_tp, BDDF1_tp, -BDDF1_tp, -BDDF1_tp, BDDF1_tp, -BDDF1_tp,
    -BDDF1_tp, BDDF1_tp, -BDDF1_tp, -BDDF1_tp, -BDDF1_tp, BDDF1_tp, 1,1,1};

// face 0
static double NF_N_H_BDDF1_3D_F0_Xi[]   = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};
static double NF_N_H_BDDF1_3D_F0_Eta[]  = {-BDDF1_tp, -BDDF1_tp, BDDF1_tp};
static double NF_N_H_BDDF1_3D_F0_Zeta[] = { -1, -1, -1};

// face 1                               1
static double NF_N_H_BDDF1_3D_F1_Xi[]   = {-BDDF1_tp, -BDDF1_tp, BDDF1_tp};
static double NF_N_H_BDDF1_3D_F1_Eta[]  = { -1, -1, -1};
static double NF_N_H_BDDF1_3D_F1_Zeta[] = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};

// face 2                               2
static double NF_N_H_BDDF1_3D_F2_Xi[]   = { 1,1,1 };
static double NF_N_H_BDDF1_3D_F2_Eta[]  = {-BDDF1_tp, -BDDF1_tp, BDDF1_tp };
static double NF_N_H_BDDF1_3D_F2_Zeta[] = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};

 //face 3                               3
static double NF_N_H_BDDF1_3D_F3_Xi[]   = { -BDDF1_tp, -BDDF1_tp, BDDF1_tp };
static double NF_N_H_BDDF1_3D_F3_Eta[]  = { 1,1,1 };
static double NF_N_H_BDDF1_3D_F3_Zeta[] = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};

// face 4                               4
static double NF_N_H_BDDF1_3D_F4_Xi[]   = { -1, -1, -1};
static double NF_N_H_BDDF1_3D_F4_Eta[]  = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};
static double NF_N_H_BDDF1_3D_F4_Zeta[] = {-BDDF1_tp, -BDDF1_tp, BDDF1_tp };

// face 5                               5
static double NF_N_H_BDDF1_3D_F5_Xi[]   = {-BDDF1_tp, -BDDF1_tp, BDDF1_tp};
static double NF_N_H_BDDF1_3D_F5_Eta[]  = {-BDDF1_tp, BDDF1_tp, -BDDF1_tp};
static double NF_N_H_BDDF1_3D_F5_Zeta[] = { 1,1,1 };

static double *NF_N_H_BDDF1_3D_XiArray[6] = {
                        NF_N_H_BDDF1_3D_F0_Xi,
                        NF_N_H_BDDF1_3D_F1_Xi,
                        NF_N_H_BDDF1_3D_F2_Xi,
                        NF_N_H_BDDF1_3D_F3_Xi,
                        NF_N_H_BDDF1_3D_F4_Xi,
                        NF_N_H_BDDF1_3D_F5_Xi };

static double *NF_N_H_BDDF1_3D_EtaArray[6] = {
                        NF_N_H_BDDF1_3D_F0_Eta,
                        NF_N_H_BDDF1_3D_F1_Eta,
                        NF_N_H_BDDF1_3D_F2_Eta,
                        NF_N_H_BDDF1_3D_F3_Eta,
                        NF_N_H_BDDF1_3D_F4_Eta,
                        NF_N_H_BDDF1_3D_F5_Eta };

static double *NF_N_H_BDDF1_3D_ZetaArray[6] = {
                        NF_N_H_BDDF1_3D_F0_Zeta,
                        NF_N_H_BDDF1_3D_F1_Zeta,
                        NF_N_H_BDDF1_3D_F2_Zeta,
                        NF_N_H_BDDF1_3D_F3_Zeta,
                        NF_N_H_BDDF1_3D_F4_Zeta,
                        NF_N_H_BDDF1_3D_F5_Zeta };

static double NF_N_H_BDDF1_3D_T[] = {-100};//???
static double NF_N_H_BDDF1_3D_S[] = {-100};//???



void NF_N_H_BDDF1_3D_EvalAll(const TCollection *, const TBaseCell *,
                             const double *PointValues, double *Functionals)
{
  //face 0
  Functionals[0]  = -PointValues[36] * 4.0;
  Functionals[1]  = -PointValues[37] * 4.0;
  Functionals[2]  = -PointValues[38] * 4.0;
  //face 1
  Functionals[3]  = -PointValues[21] * 4.0;
  Functionals[4]  = -PointValues[22] * 4.0;
  Functionals[5]  = -PointValues[23] * 4.0;
  //face 2
  Functionals[6]  =  PointValues[6]  * 4.0;
  Functionals[7]  =  PointValues[7]  * 4.0;
  Functionals[8]  =  PointValues[8]  * 4.0;
  //face 3
  Functionals[9]  =  PointValues[27] * 4.0;
  Functionals[10] =  PointValues[28] * 4.0;
  Functionals[11] =  PointValues[29] * 4.0;
  //face 4
  Functionals[12] = -PointValues[12] * 4.0;
  Functionals[13] = -PointValues[13] * 4.0;
  Functionals[14] = -PointValues[14] * 4.0;
  //face 5
  Functionals[15] =  PointValues[51] * 4.0;
  Functionals[16] =  PointValues[52] * 4.0;
  Functionals[17] =  PointValues[53] * 4.0;
}

void NF_N_H_BDDF1_3D_EvalFace(const TCollection *, const TBaseCell *Cell, int face,
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
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;
}

static int NF_N_H_BDDF1_3D_N_AllFunctionals = 18;
static int NF_N_H_BDDF1_3D_N_PointsAll = 18;
static int NF_N_H_BDDF1_3D_N_FaceFunctionals[] = { 3, 3, 3, 3, 3, 3 };
static int NF_N_H_BDDF1_3D_N_PointsFace[] = { 3, 3, 3, 3, 3, 3 };
