static double NF_C_Q_Q7_2D_Xi[] = {
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00
};

static double NF_C_Q_Q7_2D_Eta[] = {
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         7.142857142857143e-01,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00
};

static double NF_C_Q_Q7_2D_T[] = {
        -1.000000000000000e+00,
        -7.142857142857143e-01,
        -4.285714285714285e-01,
        -1.428571428571428e-01,
         1.428571428571428e-01,
         4.285714285714285e-01,
         7.142857142857143e-01,
         1.000000000000000e+00
};

void NF_C_Q_Q7_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 64*sizeof(double));
}

void NF_C_Q_Q7_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 8*sizeof(double));
}