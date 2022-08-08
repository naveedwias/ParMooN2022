static double NF_C_T_SV2_2D_Xi[] = { 0, 0.5, 1, 0.5, 0, 0, 
                                     0.16666666666666666667,
                                     0.66666666666666666667,
                                     0.16666666666666666667,
                                     0.33333333333333333333
                                   };
static double NF_C_T_SV2_2D_Eta[] = { 0, 0, 0, 0.5, 1, 0.5,
                                      0.16666666666666666667,
                                      0.16666666666666666667,
                                      0.66666666666666666667,
                                      0.33333333333333333333
                                    };
static double NF_C_T_SV2_2D_T[] = { -1, 0, 1 };

void NF_C_T_SV2_2D_EvalAll(const TCollection *, const TBaseCell *,
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
}

void NF_C_T_SV2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}
