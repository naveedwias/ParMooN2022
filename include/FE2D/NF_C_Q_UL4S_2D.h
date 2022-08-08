static double NF_C_Q_UL4S_2D_Xi[] = {
-1.0, -0.5, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.0, -0.5, -1.0,
-1.0, -1.0, -1.0, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5
};

static double NF_C_Q_UL4S_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0,
0.5, 0.0, -0.5, -0.5, 0.0, -0.5, 0.0, -0.5, 0.5
};

static double NF_C_Q_UL4S_2D_T[] = {
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00
};

void NF_C_Q_UL4S_2D_EvalAll(const TCollection *, const TBaseCell *,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 22*sizeof(double));
}

void NF_C_Q_UL4S_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                             const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 5*sizeof(double));
}
