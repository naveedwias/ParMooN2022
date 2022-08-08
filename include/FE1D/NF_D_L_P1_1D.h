
static double NF_D_L_P1_1D_Xi[] = { -0.577350269189625764520, 0.577350269189625764520};
static double NF_D_L_P1_1D_Eta[]= { 0,  0 };
static double *NF_D_L_P1_1D_T = nullptr;

void NF_D_L_P1_1D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = 0.5*( PointValues[0] + PointValues[1] );

  Functionals[1] = 0.5*( -1.077350269189625764520*PointValues[0]
                         + 0.077350269189625764520*PointValues[1] );
}

void NF_D_L_P1_1D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{

}
