static double NF_C_T_B3_2D_Xi[] = {
                0.0000000000000000, 0.3333333333333333,
                0.6666666666666667, 1.0000000000000000,
                0.0000000000000000, 0.6666666666666667,
                0.0000000000000000, 0.3333333333333333,
                0.0000000000000000,
                0.3333333333333333, 
                0.7974269853530873, 0.1012865073234563, 
                0.1012865073234563, 
                0.05971587178976982, 0.4701420641051151, 
                0.4701420641051151 };
static double NF_C_T_B3_2D_Eta[] = {
                0.0000000000000000, 0.0000000000000000,
                0.0000000000000000, 0.0000000000000000,
                0.3333333333333333, 0.3333333333333333,
                0.6666666666666667, 0.6666666666666667,
                1.0000000000000000,
                0.3333333333333333, 
                0.1012865073234563, 0.7974269853530873, 
                0.1012865073234563, 
                0.4701420641051151, 0.4701420641051151, 
                0.05971587178976982 };
static double NF_C_T_B3_2D_T[] = {
               -1.0000000000000000, -0.3333333333333333,
                0.3333333333333333,  1.0000000000000000 };

static double NF_C_T_B3_2D_W[] = {
                0.1125, 
                0.0629695902724136, 0.0629695902724136,
                0.0629695902724136, 
                0.0661970763942531, 0.0661970763942531, 
                0.0661970763942531 };

void NF_C_T_B3_2D_EvalAll(const TCollection *, const TBaseCell *,
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
  Functionals[8] = PointValues[8];
  Functionals[9] =(NF_C_T_B3_2D_W[0]*PointValues[9]
                  +NF_C_T_B3_2D_W[1]*PointValues[10]
                  +NF_C_T_B3_2D_W[2]*PointValues[11]
                  +NF_C_T_B3_2D_W[3]*PointValues[12]
                  +NF_C_T_B3_2D_W[4]*PointValues[13]
                  +NF_C_T_B3_2D_W[5]*PointValues[14]
                  +NF_C_T_B3_2D_W[6]*PointValues[15])*360;
  Functionals[10] =(NF_C_T_B3_2D_W[0]*PointValues[9]*NF_C_T_B3_2D_Xi[9]
                   +NF_C_T_B3_2D_W[1]*PointValues[10]*NF_C_T_B3_2D_Xi[10]
                   +NF_C_T_B3_2D_W[2]*PointValues[11]*NF_C_T_B3_2D_Xi[11]
                   +NF_C_T_B3_2D_W[3]*PointValues[12]*NF_C_T_B3_2D_Xi[12]
                   +NF_C_T_B3_2D_W[4]*PointValues[13]*NF_C_T_B3_2D_Xi[13]
                   +NF_C_T_B3_2D_W[5]*PointValues[14]*NF_C_T_B3_2D_Xi[14]
                   +NF_C_T_B3_2D_W[6]*PointValues[15]*NF_C_T_B3_2D_Xi[15])*360;
  Functionals[11] =(NF_C_T_B3_2D_W[0]*PointValues[9]*NF_C_T_B3_2D_Eta[9]
                   +NF_C_T_B3_2D_W[1]*PointValues[10]*NF_C_T_B3_2D_Eta[10]
                   +NF_C_T_B3_2D_W[2]*PointValues[11]*NF_C_T_B3_2D_Eta[11]
                   +NF_C_T_B3_2D_W[3]*PointValues[12]*NF_C_T_B3_2D_Eta[12]
                   +NF_C_T_B3_2D_W[4]*PointValues[13]*NF_C_T_B3_2D_Eta[13]
                   +NF_C_T_B3_2D_W[5]*PointValues[14]*NF_C_T_B3_2D_Eta[14]
                   +NF_C_T_B3_2D_W[6]*PointValues[15]*NF_C_T_B3_2D_Eta[15])*360;
}

void NF_C_T_B3_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
}
