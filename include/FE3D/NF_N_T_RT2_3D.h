// ***********************************************************************
// Raviart-Thomas element of second order on hexahedra, 3D
// ***********************************************************************

/* for all functionals */
static double NF_N_T_RT2_3D_Xi[]  = {
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0,0,0,0,0,0,
    //inner dof
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0,
    0.333333333333333333333333333333333,
    0.25,
    0.909090909090909090909090909090909e-1,
    0.909090909090909090909090909090909e-1,
    0.727272727272727272727272727272727,
    0.909090909090909090909090909090909e-1,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697
};
static double NF_N_T_RT2_3D_Eta[] = {
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0,0,0,0,0,0,
    0.6, 0.4, 0.2, 0.4, 0.2, 0.2,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    //inner dof
    0.333333333333333333333333333333333,
    0,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.25,
    0.909090909090909090909090909090909e-1,
    0.727272727272727272727272727272727,
    0.909090909090909090909090909090909e-1,
    0.909090909090909090909090909090909e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1
};
static double NF_N_T_RT2_3D_Zeta[]= {
    0,0,0,0,0,0,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    //inner dof
    0,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.25,
    0.727272727272727272727272727272727,
    0.90909090909090909090909090909091e-1,
    0.909090909090909090909090909090911e-1,
    0.909090909090909090909090909090911e-1,
    0.433449846426335701760119535773699,
    0.433449846426335701760119535773699,
    0.66550153573664298239880464226304e-1,
    0.433449846426335701760119535773699,
    0.66550153573664298239880464226304e-1,
    0.665501535736642982398804642263035e-1
};

static double NF_N_T_RT2_3D_Weights[]= {
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.302836780970891758063769725577305e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1};

/* face 0                               0 */
static double NF_N_T_RT2_3D_F0_Xi[]   = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_RT2_3D_F0_Eta[]  = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};
static double NF_N_T_RT2_3D_F0_Zeta[] = { 0,0,0,0,0,0 };

/* face 1                               1 */
static double NF_N_T_RT2_3D_F1_Xi[]   = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};
static double NF_N_T_RT2_3D_F1_Eta[]  = { 0,0,0,0,0,0 };
static double NF_N_T_RT2_3D_F1_Zeta[] = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};

/* face 2                               2 */
static double NF_N_T_RT2_3D_F2_Xi[]   = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_RT2_3D_F2_Eta[]  = {0.6, 0.4, 0.2, 0.4, 0.2, 0.2};
static double NF_N_T_RT2_3D_F2_Zeta[] = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};

/* face 3                               3 */
static double NF_N_T_RT2_3D_F3_Xi[]   = { 0,0,0,0,0,0 };
static double NF_N_T_RT2_3D_F3_Eta[]  = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_RT2_3D_F3_Zeta[] = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};

static double *NF_N_T_RT2_3D_XiArray[4] = {
                        NF_N_T_RT2_3D_F0_Xi,
                        NF_N_T_RT2_3D_F1_Xi,
                        NF_N_T_RT2_3D_F2_Xi,
                        NF_N_T_RT2_3D_F3_Xi };

static double *NF_N_T_RT2_3D_EtaArray[4] = {
                        NF_N_T_RT2_3D_F0_Eta,
                        NF_N_T_RT2_3D_F1_Eta,
                        NF_N_T_RT2_3D_F2_Eta,
                        NF_N_T_RT2_3D_F3_Eta };

static double *NF_N_T_RT2_3D_ZetaArray[4] = {
                        NF_N_T_RT2_3D_F0_Zeta,
                        NF_N_T_RT2_3D_F1_Zeta,
                        NF_N_T_RT2_3D_F2_Zeta,
                        NF_N_T_RT2_3D_F3_Zeta };

static double NF_N_T_RT2_3D_T[1] = {};// ???
static double NF_N_T_RT2_3D_S[1] = {};// ???

void NF_N_T_RT2_3D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  // PointValues[4*i + j] means i-th component (i=0 for x, i=1 for y, i=2 for z)
  // at j-th evaluation point (see NF_N_T_RT2_3D_Xi, ...Eta, ...Zeta)
  //face 0
  Functionals[0]  = -PointValues[78];
  Functionals[1]  = -PointValues[79];
  Functionals[2]  = -PointValues[80];
  Functionals[3]  = -PointValues[81];
  Functionals[4]  = -PointValues[82];
  Functionals[5]  = -PointValues[83];
  //face 1
  Functionals[6]  = -PointValues[45];
  Functionals[7]  = -PointValues[46];
  Functionals[8]  = -PointValues[47];
  Functionals[9]  = -PointValues[48];
  Functionals[10] = -PointValues[49];
  Functionals[11] = -PointValues[50];
  //face 2
  Functionals[12] =  PointValues[12]+PointValues[51]+PointValues[90];
  Functionals[13] =  PointValues[13]+PointValues[52]+PointValues[91];
  Functionals[14] =  PointValues[14]+PointValues[53]+PointValues[92];
  Functionals[15] =  PointValues[15]+PointValues[54]+PointValues[93];
  Functionals[16] =  PointValues[16]+PointValues[55]+PointValues[94];
  Functionals[17] =  PointValues[17]+PointValues[56]+PointValues[95];
  //face 3
  Functionals[18] = -PointValues[18];
  Functionals[19] = -PointValues[19];
  Functionals[20] = -PointValues[20];
  Functionals[21] = -PointValues[21];
  Functionals[22] = -PointValues[22];
  Functionals[23] = -PointValues[23];

  //inner dofs
  int i;
  double s;

  //x-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[24] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+24] * NF_N_T_RT2_3D_Xi[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[25] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+24] * NF_N_T_RT2_3D_Eta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[26] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+24] * NF_N_T_RT2_3D_Zeta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[27] = s;

  //y-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+63] * NF_N_T_RT2_3D_Weights[i];
  Functionals[28] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+63] * NF_N_T_RT2_3D_Xi[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[29] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+63] * NF_N_T_RT2_3D_Eta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[30] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+63] * NF_N_T_RT2_3D_Zeta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[31] = s;

  //z-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+102] * NF_N_T_RT2_3D_Weights[i];
  Functionals[32] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+102] * NF_N_T_RT2_3D_Xi[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[33] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+102] * NF_N_T_RT2_3D_Eta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[34] = s;
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+102] * NF_N_T_RT2_3D_Zeta[i+24] * NF_N_T_RT2_3D_Weights[i];
  Functionals[35] = s;
}

void NF_N_T_RT2_3D_EvalFace(const TCollection *, const TBaseCell *Cell, int face,
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
  // compute measure of this face
  s = std::sqrt( std::pow((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + std::pow((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + std::pow((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  for(int i=0; i<6; i++)
    Functionals[i] = PointValues[i]*s;
}

static int NF_N_T_RT2_3D_N_AllFunctionals = 36;
static int NF_N_T_RT2_3D_N_PointsAll = 39;
static int NF_N_T_RT2_3D_N_FaceFunctionals[] = { 6, 6, 6, 6 };
static int NF_N_T_RT2_3D_N_PointsFace[] = { 6, 6, 6, 6 };
