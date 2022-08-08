

static double NF_C_L_P0_1D_Xi[] = { 0 };
static double NF_C_L_P0_1D_Eta[] = { 0 };
static double NF_C_L_P0_1D_T[] = { -1, 1 };

void NF_C_L_P0_1D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_C_L_P0_1D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{

}
