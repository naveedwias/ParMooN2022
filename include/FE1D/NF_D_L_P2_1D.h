
static double NF_D_L_P2_1D_Xi[] = { -0.77459666924148337703585307995647992,
                                     0,
                                     0.77459666924148337703585307995647992 };
static double NF_D_L_P2_1D_Eta[] = { 0, 0, 0 };
static double *NF_D_L_P2_1D_T = nullptr;

void NF_D_L_P2_1D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{

  Functionals[0] =0.5*(0.55555555555555560*PointValues[0]
                       + 0.88888888888888880*PointValues[1]
                       + 0.55555555555555560*PointValues[2]   );

  Functionals[1] =0.5*(-1.27459666924148337703585307995647992*PointValues[0]
                       + 0.5*PointValues[1]
                       + 0.27459666924148337703585307995647992 *PointValues[2]  );

  Functionals[2] = 0.5*(- 0.27459666924148337703585307995647992*PointValues[0]
                       + 0.5*PointValues[1]
                       + 1.27459666924148337703585307995647992*PointValues[2]  );
}

void NF_D_L_P2_1D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{

}
