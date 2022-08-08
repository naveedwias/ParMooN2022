static double NF_C_Q_Q6_2D_Xi[] = {
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00
};

static double NF_C_Q_Q6_2D_Eta[] = {
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         3.333333333333333e-01,
         3.333333333333333e-01,
         3.333333333333333e-01,
         3.333333333333333e-01,
         3.333333333333333e-01,
         3.333333333333333e-01,
         3.333333333333333e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         6.666666666666666e-01,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00
};

static double NF_C_Q_Q6_2D_T[] = {
        -1.000000000000000e+00,
        -6.666666666666666e-01,
        -3.333333333333333e-01,
         0.000000000000000e+00,
         3.333333333333333e-01,
         6.666666666666666e-01,
         1.000000000000000e+00
};

void NF_C_Q_Q6_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 49*sizeof(double));
}

void NF_C_Q_Q6_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 7*sizeof(double));
}
