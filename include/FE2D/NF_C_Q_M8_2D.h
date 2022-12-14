static double NF_C_Q_M8_2D_Xi[] = {
-1.0, -3.0/4.0, -1.0/2.0, -1.0/4.0, 0.0, 1.0/4.0, 1.0/2.0, 3.0/4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0/4.0, 1.0/2.0, 1.0/4.0, 0.0, -1.0/4.0, -1.0/2.0, -3.0/4.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0/3.0, -1.0/3.0, -2.0/3.0, 0.0, -1.0/3.0, -2.0/3.0, 1.0/3.0, 0.0, -1.0/3.0, -2.0/3.0, 2.0/3.0, 1.0/3.0, 0.0, -1.0/3.0, -2.0/3.0
};

static double NF_C_Q_M8_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -3.0/4.0, -1.0/2.0, -1.0/4.0, 0.0, 1.0/4.0, 1.0/2.0, 3.0/4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0/4.0, 1.0/2.0, 1.0/4.0, 0.0, -1.0/4.0, -1.0/2.0, -3.0/4.0, -2.0/3.0, -2.0/3.0, -1.0/3.0, -2.0/3.0, -1.0/3.0, 0.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0, 2.0/3.0
};

static double NF_C_Q_M8_2D_T[] = { -1.0, -3.0/4.0, -1.0/2.0, -1.0/4.0, 0.0, 1.0/4.0, 1.0/2.0, 3.0/4.0, 1.0 };

void NF_C_Q_M8_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 47*sizeof(double));
}

void NF_C_Q_M8_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 9*sizeof(double));
}
