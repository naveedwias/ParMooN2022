static double NF_C_T_P2_2D_Xi[] = { 0, 0.5, 1, 0, 0.5, 0 };
static double NF_C_T_P2_2D_Eta[] = { 0, 0, 0, 0.5, 0.5, 1 };
static double NF_C_T_P2_2D_T[] = { -1, 0, 1 };

void NF_C_T_P2_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
}

void NF_C_T_P2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}
