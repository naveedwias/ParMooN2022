/* for all functionals */
static double NF_C_T_B2_3D_Xi[]   = { 0, 0.5, 1, 0, 0.5, 0,
                                      0, 0.5, 0,
                                      0,
       0.33333333333333333333, 0.33333333333333333333, 0.33333333333333333333, 0,
                                      0.25  };
static double NF_C_T_B2_3D_Eta[]  = { 0, 0, 0, 0.5, 0.5, 1,
                                      0, 0, 0.5,
                                      0, 
    0.33333333333333333333, 0, 0.33333333333333333333, 0.33333333333333333333,   
                                      0.25 };
static double NF_C_T_B2_3D_Zeta[] = { 0, 0, 0, 0, 0, 0,
                                      0.5, 0.5, 0.5,
                                      1, 
    0, 0.33333333333333333333, 0.33333333333333333333, 0.33333333333333333333,    
                                      0.25 };

/* face 0                               0,   1,  2,  3,  4,  5, 10*/
static double NF_C_T_B2_3D_F0_Xi[]   = {  0, 0.5, 1, 0, 0.5, 0, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F0_Eta[]  = {  0, 0, 0, 0.5, 0.5, 1, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F0_Zeta[] = {  0, 0, 0, 0, 0, 0, 0 };

/* face 1                               0,  6,  9,  1,  7,  2,  11 */
static double NF_C_T_B2_3D_F1_Xi[]   = {  0, 0, 0, 0.5, 0.5, 1, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F1_Eta[]  = {  0, 0, 0, 0, 0, 0, 0 };
static double NF_C_T_B2_3D_F1_Zeta[] = {  0, 0.5, 1, 0, 0.5, 0, 0.33333333333333333333 };

/* face 2                               5,  4,  2,  8,  7,  9, 12 */
static double NF_C_T_B2_3D_F2_Xi[]   = {  0, 0.5, 1, 0, 0.5, 0, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F2_Eta[]  = {  1, 0.5, 0, 0.5, 0, 0, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F2_Zeta[] = {  0, 0, 0, 0.5, 0.5, 1, 0.33333333333333333333 };

/* face 3                               0,  3,  5,  6,  8,  9, 13 */
static double NF_C_T_B2_3D_F3_Xi[]   = {  0, 0, 0, 0, 0, 0, 0 };
static double NF_C_T_B2_3D_F3_Eta[]  = {  0, 0.5, 1, 0, 0.5, 0, 0.33333333333333333333 };
static double NF_C_T_B2_3D_F3_Zeta[] = {  0, 0, 0, 0.5, 0.5, 1, 0.33333333333333333333 };

static double *NF_C_T_B2_3D_XiArray[4] = { 
                        NF_C_T_B2_3D_F0_Xi,
                        NF_C_T_B2_3D_F1_Xi,
                        NF_C_T_B2_3D_F2_Xi,
                        NF_C_T_B2_3D_F3_Xi };

static double *NF_C_T_B2_3D_EtaArray[4] = { 
                        NF_C_T_B2_3D_F0_Eta,
                        NF_C_T_B2_3D_F1_Eta,
                        NF_C_T_B2_3D_F2_Eta,
                        NF_C_T_B2_3D_F3_Eta };

static double *NF_C_T_B2_3D_ZetaArray[4] = { 
                        NF_C_T_B2_3D_F0_Zeta,
                        NF_C_T_B2_3D_F1_Zeta,
                        NF_C_T_B2_3D_F2_Zeta,
                        NF_C_T_B2_3D_F3_Zeta };

static double NF_C_T_B2_3D_T[7] = { 0, 0.5, 1,   0, 0.5, 0, 0.33333333333333333333 };
static double NF_C_T_B2_3D_S[7] = { 0,   0, 0, 0.5, 0.5, 1, 0.33333333333333333333 };

void NF_C_T_B2_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 15*sizeof(double));
}

void NF_C_T_B2_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 7*sizeof(double));
}

static int NF_C_T_B2_3D_N_AllFunctionals = 15;
static int NF_C_T_B2_3D_N_PointsAll = 15;
static int NF_C_T_B2_3D_N_FaceFunctionals[] = { 7, 7, 7, 7, 7, 7, 7};
static int NF_C_T_B2_3D_N_PointsFace[] = {  7, 7, 7, 7, 7, 7, 7};
