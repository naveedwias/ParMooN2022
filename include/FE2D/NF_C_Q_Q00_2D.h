static double NF_C_Q_Q00_2D_Xi[] = { 0 };
static double NF_C_Q_Q00_2D_Eta[] = { 0 };
static double *NF_C_Q_Q00_2D_T = nullptr;

void NF_C_Q_Q00_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_C_Q_Q00_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *, double *)
{
}
