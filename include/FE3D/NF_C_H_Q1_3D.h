/* for all functionals */
static double NF_C_H_Q1_3D_Xi[]   = { -1,  1, -1,  1, -1,  1, -1,  1 };
static double NF_C_H_Q1_3D_Eta[]  = { -1, -1,  1,  1, -1, -1,  1,  1 };
static double NF_C_H_Q1_3D_Zeta[] = { -1, -1, -1, -1,  1,  1,  1,  1 };

/* face 0                               0   1   2   3 */
static double NF_C_H_Q1_3D_F0_Xi[]   = { -1,  1, -1,  1 };
static double NF_C_H_Q1_3D_F0_Eta[]  = { -1, -1,  1,  1 };
static double NF_C_H_Q1_3D_F0_Zeta[] = { -1, -1, -1, -1 };

/* face 1                               0   4   1   5 */
static double NF_C_H_Q1_3D_F1_Xi[]   = { -1, -1,  1,  1 };
static double NF_C_H_Q1_3D_F1_Eta[]  = { -1, -1, -1, -1 };
static double NF_C_H_Q1_3D_F1_Zeta[] = { -1,  1, -1,  1 };

/* face 2                               1   5   3   7 */
static double NF_C_H_Q1_3D_F2_Xi[]   = {  1,  1,  1,  1 };
static double NF_C_H_Q1_3D_F2_Eta[]  = { -1, -1,  1,  1 };
static double NF_C_H_Q1_3D_F2_Zeta[] = { -1,  1, -1,  1 };

/* face 3                               3   7   2   6 */
static double NF_C_H_Q1_3D_F3_Xi[]   = {  1,  1, -1, -1 };
static double NF_C_H_Q1_3D_F3_Eta[]  = {  1,  1,  1,  1 };
static double NF_C_H_Q1_3D_F3_Zeta[] = { -1,  1, -1,  1 };

/* face 4                               0   2   4   6 */
static double NF_C_H_Q1_3D_F4_Xi[]   = { -1, -1, -1, -1 };
static double NF_C_H_Q1_3D_F4_Eta[]  = { -1,  1, -1,  1 };
static double NF_C_H_Q1_3D_F4_Zeta[] = { -1, -1,  1,  1 };

/* face 5                               4   6   5   7 */
static double NF_C_H_Q1_3D_F5_Xi[]   = { -1, -1,  1,  1 };
static double NF_C_H_Q1_3D_F5_Eta[]  = { -1,  1, -1,  1 };
static double NF_C_H_Q1_3D_F5_Zeta[] = {  1,  1,  1,  1 };

static double *NF_C_H_Q1_3D_XiArray[6] = { 
                        NF_C_H_Q1_3D_F0_Xi,
                        NF_C_H_Q1_3D_F1_Xi,
                        NF_C_H_Q1_3D_F2_Xi,
                        NF_C_H_Q1_3D_F3_Xi,
                        NF_C_H_Q1_3D_F4_Xi,
                        NF_C_H_Q1_3D_F5_Xi };

static double *NF_C_H_Q1_3D_EtaArray[6] = { 
                        NF_C_H_Q1_3D_F0_Eta,
                        NF_C_H_Q1_3D_F1_Eta,
                        NF_C_H_Q1_3D_F2_Eta,
                        NF_C_H_Q1_3D_F3_Eta,
                        NF_C_H_Q1_3D_F4_Eta,
                        NF_C_H_Q1_3D_F5_Eta };

static double *NF_C_H_Q1_3D_ZetaArray[6] = { 
                        NF_C_H_Q1_3D_F0_Zeta,
                        NF_C_H_Q1_3D_F1_Zeta,
                        NF_C_H_Q1_3D_F2_Zeta,
                        NF_C_H_Q1_3D_F3_Zeta,
                        NF_C_H_Q1_3D_F4_Zeta,
                        NF_C_H_Q1_3D_F5_Zeta };

static double NF_C_H_Q1_3D_T[4] = {  0, 1, 0, 1 };
static double NF_C_H_Q1_3D_S[4] = {  0, 0, 1, 1 };

void NF_C_H_Q1_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];
}

void NF_C_H_Q1_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
}

static int NF_C_H_Q1_3D_N_AllFunctionals = 8;
static int NF_C_H_Q1_3D_N_PointsAll = 8;
static int NF_C_H_Q1_3D_N_FaceFunctionals[] = { 4, 4, 4, 4, 4, 4 };
static int NF_C_H_Q1_3D_N_PointsFace[] = { 4, 4, 4, 4, 4, 4 };
