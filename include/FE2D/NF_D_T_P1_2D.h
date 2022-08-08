static double NF_D_T_P1_2D_Xi[] = { 0.5, 0.5, 0 };
static double NF_D_T_P1_2D_Eta[] = { 0, 0.5, 0.5 };
static double NF_D_T_P1_2D_T_P[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

/*
   weighting functions: 2, 24*xi-8, 24*eta-8
*/

void NF_D_T_P1_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = 2*(PointValues[0]+PointValues[1]+PointValues[2])/6;
  Functionals[1] = ( 4*PointValues[0]
                    +4*PointValues[1]
                    -8*PointValues[2])/6;
  Functionals[2] = (-8*PointValues[0]
                    +4*PointValues[1]
                    +4*PointValues[2])/6;
}

void NF_D_T_P1_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}
