static double NF_C_Q_M3_2D_Xi[] = {
-1.0, -1.0/3.0, 1.0/3.0, 1.0, 1.0, 1.0, 1.0, 1.0/3.0, -1.0/3.0, -1.0, -1.0, -1.0 };

static double NF_C_Q_M3_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -1.0/3.0, 1.0/3.0, 1.0, 1.0, 1.0, 1.0, 1.0/3.0, -1.0/3.0 };

static double NF_C_Q_M3_2D_T[] = { -1.0, -1.0/3.0, 1.0/3.0, 1.0 };

void NF_C_Q_M3_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 12*sizeof(double));
}

void NF_C_Q_M3_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 4*sizeof(double));
}
