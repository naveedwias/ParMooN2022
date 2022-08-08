static double NF_C_T_P1MINI_2D_Xi[7] = { 0, 0.5, 1, 0, 0.5, 0,
                                    0.33333333333333333333333 };
static double NF_C_T_P1MINI_2D_Eta[7] = { 0, 0, 0, 0.5, 0.5, 1,
                                     0.33333333333333333333333 };
static double NF_C_T_P1MINI_2D_T[2] = { -1, 1 };

void NF_C_T_P1MINI_2D_EvalAll(const TCollection *, const TBaseCell *,
                              const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[2];
  Functionals[2] = PointValues[5];

  Functionals[3] = (27*PointValues[6]
                  + 8*(PointValues[1]+PointValues[3]+PointValues[4])
                  + 3*(PointValues[0]+PointValues[2]+PointValues[5]))/27;
}

void NF_C_T_P1MINI_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                               const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}
