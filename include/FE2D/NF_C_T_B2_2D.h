static double NF_C_T_B2_2D_Xi[] = { 0, 0.5, 1,
                                  0, 0.5,
                                  0, 
                                  0.33333333333333333333 };
static double NF_C_T_B2_2D_Eta[] = { 0, 0, 0,
                                   0.5, 0.5,
                                   1,
                                   0.33333333333333333333 };
static double NF_C_T_B2_2D_T[] = { -1, 0, 1 };

void NF_C_T_B2_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = 1.35*PointValues[6]
                  +0.15*(PointValues[0]+PointValues[2]+PointValues[5])
                  +0.40*(PointValues[1]+PointValues[3]+PointValues[4]);
}

void NF_C_T_B2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}
