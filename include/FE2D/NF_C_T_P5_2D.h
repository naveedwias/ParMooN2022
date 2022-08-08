static double NF_C_T_P5_2D_Xi[] = { 0, 0.2, 0.4, 0.6, 0.8, 1,
                                    0, 0.2, 0.4, 0.6, 0.8,
                                    0, 0.2, 0.4, 0.6,
                                    0, 0.2, 0.4,
                                    0, 0.2,
                                    0 };
static double NF_C_T_P5_2D_Eta[] = { 0,   0,   0,   0,   0,   0,
                                     0.2, 0.2, 0.2, 0.2, 0.2,
                                     0.4, 0.4, 0.4, 0.4,
                                     0.6, 0.6, 0.6,
                                     0.8, 0.8,
                                     1 };
static double NF_C_T_P5_2D_T[] = { -1, -0.6, -0.2, 0.2, 0.6, 1 };

void NF_C_T_P5_2D_EvalAll(const TCollection *, const TBaseCell *,
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
  Functionals[9] = PointValues[9];
  Functionals[10] = PointValues[10];
  Functionals[11] = PointValues[11];
  Functionals[12] = PointValues[12];
  Functionals[13] = PointValues[13];
  Functionals[14] = PointValues[14];
  Functionals[15] = PointValues[15];
  Functionals[16] = PointValues[16];
  Functionals[17] = PointValues[17];
  Functionals[18] = PointValues[18];
  Functionals[19] = PointValues[19];
  Functionals[20] = PointValues[20];
}

void NF_C_T_P5_2D_EvalEdge(const TCollection *, const TBaseCell *, int ,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
}
