// ***********************************************************************
// P2 element, discontinuous, 3D
// ***********************************************************************
/* for all functionals */
static double NF_D_T_P2_3D_Xi[]   = {
  0.25,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462};
static double NF_D_T_P2_3D_Eta[]  = {
  0.25,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.0714285714285714285714285714286,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538};
static double NF_D_T_P2_3D_Zeta[] = {
  0.25,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538 };

static double NF_D_T_P2_3D_Weights[]= {
   -0.0131555555555555555555555555556,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889
};

/* face 0                               0 */
static double *NF_D_T_P2_3D_F0_Xi = nullptr;
static double *NF_D_T_P2_3D_F0_Eta = nullptr;
static double *NF_D_T_P2_3D_F0_Zeta = nullptr;

/* face 1                               1 */
static double *NF_D_T_P2_3D_F1_Xi = nullptr;
static double *NF_D_T_P2_3D_F1_Eta = nullptr;
static double *NF_D_T_P2_3D_F1_Zeta = nullptr;

/* face 2                               2 */
static double *NF_D_T_P2_3D_F2_Xi = nullptr;
static double *NF_D_T_P2_3D_F2_Eta = nullptr;
static double *NF_D_T_P2_3D_F2_Zeta = nullptr;

/* face 3                               3 */
static double *NF_D_T_P2_3D_F3_Xi = nullptr;
static double *NF_D_T_P2_3D_F3_Eta = nullptr;
static double *NF_D_T_P2_3D_F3_Zeta = nullptr;

static double *NF_D_T_P2_3D_XiArray[4] = {
                        NF_D_T_P2_3D_F0_Xi,
                        NF_D_T_P2_3D_F1_Xi,
                        NF_D_T_P2_3D_F2_Xi,
                        NF_D_T_P2_3D_F3_Xi };

static double *NF_D_T_P2_3D_EtaArray[4] = {
                        NF_D_T_P2_3D_F0_Eta,
                        NF_D_T_P2_3D_F1_Eta,
                        NF_D_T_P2_3D_F2_Eta,
                        NF_D_T_P2_3D_F3_Eta };

static double *NF_D_T_P2_3D_ZetaArray[4] = {
                        NF_D_T_P2_3D_F0_Zeta,
                        NF_D_T_P2_3D_F1_Zeta,
                        NF_D_T_P2_3D_F2_Zeta,
                        NF_D_T_P2_3D_F3_Zeta };

static double *NF_D_T_P2_3D_T = nullptr;
static double *NF_D_T_P2_3D_S = nullptr;

void NF_D_T_P2_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  int i;
  double s;

  s = 0;
  //constant
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[0] = s;

  s = 0;
  //x
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[1] = s;

  s = 0;
  //y
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[2] = s;

  s = 0;
  //z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[3] = s;

  s = 0;
  //x*x
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[4] = s;

  s = 0;
  //x*y
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[5] = s;

  s = 0;
  //x*z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[6] = s;

  s = 0;
  //yy
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[7] = s;

  s = 0;
  //y*z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[8] = s;

  s = 0;
  //zz
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[9] = s;
}

void NF_D_T_P2_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}

static int NF_D_T_P2_3D_N_AllFunctionals = 10;
static int NF_D_T_P2_3D_N_PointsAll = 11;
static int NF_D_T_P2_3D_N_FaceFunctionals[] = { 0, 0, 0, 0 };
static int NF_D_T_P2_3D_N_PointsFace[] = { 0, 0, 0, 0 };
