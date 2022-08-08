static double NF_C_Q_UL3S_2D_Xi[] = {
-1.0, -0.3333333333333333, 0.3333333333333333, 1.0, 1.0, 1.0, 1.0,
0.3333333333333333, -0.3333333333333333, -1.0, -1.0, -1.0,
-0.3333333333333333, 0.3333333333333333, -0.3333333333333333 };

static double NF_C_Q_UL3S_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -0.3333333333333333, 0.3333333333333333, 1.0,
1.0, 1.0, 1.0, 0.3333333333333333, -0.3333333333333333,
-0.3333333333333333, -0.3333333333333333, 0.3333333333333333 };

static double NF_C_Q_UL3S_2D_T[] = { -1, -0.33333333333333333333,
                                 0.33333333333333333333, 1 };

void NF_C_Q_UL3S_2D_EvalAll(const TCollection *, const TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 15*sizeof(double));
}

void NF_C_Q_UL3S_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 4*sizeof(double));
}
