static double NF_C_Q_UL5S_2D_Xi[] = {
-1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.2,
-0.2, -0.6, -1.0, -1.0, -1.0, -1.0, -1.0, -0.6, -0.6, -0.6, -0.2, -0.2,
-0.2, 0.2, 0.2, 0.2, 0.6, -0.6
};

static double NF_C_Q_UL5S_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 0.6, 0.2, -0.2, -0.6, -0.6, -0.2, 0.2, -0.6, -0.2, 0.2,
-0.6, -0.2, 0.2, -0.6, 0.6
};

static double NF_C_Q_UL5S_2D_T[] = {
        -1.000000000000000e+00,
        -6.000000000000000e-01,
        -2.000000000000000e-01,
         2.000000000000000e-01,
         6.000000000000000e-01,
         1.000000000000000e+00
};

void NF_C_Q_UL5S_2D_EvalAll(const TCollection *, const TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 31*sizeof(double));
}

void NF_C_Q_UL5S_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 6*sizeof(double));
}
