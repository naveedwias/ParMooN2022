// ***********************************************************************
// P1 element, discontinuous, 3D
// ***********************************************************************
/* for all functionals */
static double NF_D_T_P1_3D_Xi[]   = {0, 1, 0, 0, 0.25 };
static double NF_D_T_P1_3D_Eta[]  = {0, 0, 1, 0, 0.25 };
static double NF_D_T_P1_3D_Zeta[] = {0, 0, 0, 1, 0.25 };

/* face 0                               0 */
static double *NF_D_T_P1_3D_F0_Xi = nullptr;
static double *NF_D_T_P1_3D_F0_Eta = nullptr;
static double *NF_D_T_P1_3D_F0_Zeta = nullptr;

/* face 1                               1 */
static double *NF_D_T_P1_3D_F1_Xi = nullptr;
static double *NF_D_T_P1_3D_F1_Eta = nullptr;
static double *NF_D_T_P1_3D_F1_Zeta = nullptr;

/* face 2                               2 */
static double *NF_D_T_P1_3D_F2_Xi = nullptr;
static double *NF_D_T_P1_3D_F2_Eta = nullptr;
static double *NF_D_T_P1_3D_F2_Zeta = nullptr;

/* face 3                               3 */
static double *NF_D_T_P1_3D_F3_Xi = nullptr;
static double *NF_D_T_P1_3D_F3_Eta = nullptr;
static double *NF_D_T_P1_3D_F3_Zeta = nullptr;

static double *NF_D_T_P1_3D_XiArray[4] = {
                        NF_D_T_P1_3D_F0_Xi,
                        NF_D_T_P1_3D_F1_Xi,
                        NF_D_T_P1_3D_F2_Xi,
                        NF_D_T_P1_3D_F3_Xi };

static double *NF_D_T_P1_3D_EtaArray[4] = {
                        NF_D_T_P1_3D_F0_Eta,
                        NF_D_T_P1_3D_F1_Eta,
                        NF_D_T_P1_3D_F2_Eta,
                        NF_D_T_P1_3D_F3_Eta };

static double *NF_D_T_P1_3D_ZetaArray[4] = {
                        NF_D_T_P1_3D_F0_Zeta,
                        NF_D_T_P1_3D_F1_Zeta,
                        NF_D_T_P1_3D_F2_Zeta,
                        NF_D_T_P1_3D_F3_Zeta };

static double *NF_D_T_P1_3D_T = nullptr;
static double *NF_D_T_P1_3D_S = nullptr;

void NF_D_T_P1_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = ( 16*PointValues[4] + PointValues[0] + PointValues[1]
                     + PointValues[2] + PointValues[3] )/120;
  Functionals[1] =( -30*PointValues[0] +90*PointValues[1]
                    -30*PointValues[2] -30*PointValues[3] )/120;
  Functionals[2] =( -30*PointValues[0] -30*PointValues[1]
                    +90*PointValues[2] -30*PointValues[3] )/120;
  Functionals[3] =( -30*PointValues[0] -30*PointValues[1]
                    -30*PointValues[2] +90*PointValues[3] )/120;
}

void NF_D_T_P1_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}

static int NF_D_T_P1_3D_N_AllFunctionals = 4;
static int NF_D_T_P1_3D_N_PointsAll = 5;
static int NF_D_T_P1_3D_N_FaceFunctionals[] = { 0, 0, 0, 0 };
static int NF_D_T_P1_3D_N_PointsFace[] = { 0, 0, 0, 0 };
