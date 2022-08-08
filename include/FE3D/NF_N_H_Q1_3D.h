/* for all functionals */
static double NF_N_H_Q1_3D_Xi[]   = {  0,  0,  1,  0, -1,  0 };
static double NF_N_H_Q1_3D_Eta[]  = {  0, -1,  0,  1,  0,  0 };
static double NF_N_H_Q1_3D_Zeta[] = { -1,  0,  0,  0,  0,  1 };

/* face 0                               0 */
static double NF_N_H_Q1_3D_F0_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F0_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F0_Zeta[] = { -1 };

/* face 1                               1 */
static double NF_N_H_Q1_3D_F1_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F1_Eta[]  = { -1 };
static double NF_N_H_Q1_3D_F1_Zeta[] = {  0 };

/* face 2                               2 */
static double NF_N_H_Q1_3D_F2_Xi[]   = {  1 };
static double NF_N_H_Q1_3D_F2_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F2_Zeta[] = {  0 };

/* face 3                               3 */
static double NF_N_H_Q1_3D_F3_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F3_Eta[]  = {  1 };
static double NF_N_H_Q1_3D_F3_Zeta[] = {  0 };

/* face 4                               4 */
static double NF_N_H_Q1_3D_F4_Xi[]   = { -1 };
static double NF_N_H_Q1_3D_F4_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F4_Zeta[] = {  0 };

/* face 5                               5 */
static double NF_N_H_Q1_3D_F5_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F5_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F5_Zeta[] = {  1 };

static double *NF_N_H_Q1_3D_XiArray[6] = { 
                        NF_N_H_Q1_3D_F0_Xi,
                        NF_N_H_Q1_3D_F1_Xi,
                        NF_N_H_Q1_3D_F2_Xi,
                        NF_N_H_Q1_3D_F3_Xi,
                        NF_N_H_Q1_3D_F4_Xi,
                        NF_N_H_Q1_3D_F5_Xi };

static double *NF_N_H_Q1_3D_EtaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Eta,
                        NF_N_H_Q1_3D_F1_Eta,
                        NF_N_H_Q1_3D_F2_Eta,
                        NF_N_H_Q1_3D_F3_Eta,
                        NF_N_H_Q1_3D_F4_Eta,
                        NF_N_H_Q1_3D_F5_Eta };

static double *NF_N_H_Q1_3D_ZetaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Zeta,
                        NF_N_H_Q1_3D_F1_Zeta,
                        NF_N_H_Q1_3D_F2_Zeta,
                        NF_N_H_Q1_3D_F3_Zeta,
                        NF_N_H_Q1_3D_F4_Zeta,
                        NF_N_H_Q1_3D_F5_Zeta };

static double NF_N_H_Q1_3D_T[1] = { 0.5 };
static double NF_N_H_Q1_3D_S[1] = { 0.5 };

void NF_N_H_Q1_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
}

void NF_N_H_Q1_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

static int NF_N_H_Q1_3D_N_AllFunctionals = 6;
static int NF_N_H_Q1_3D_N_PointsAll = 6;
static int NF_N_H_Q1_3D_N_FaceFunctionals[] = { 1, 1, 1, 1, 1, 1 };
static int NF_N_H_Q1_3D_N_PointsFace[] = { 1, 1, 1, 1, 1, 1 };
