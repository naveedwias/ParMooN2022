// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************
static double NF_B_Q_IB2_2D_Xi[] = { 0 };
static double NF_B_Q_IB2_2D_Eta[] = { 0 };
static double *NF_B_Q_IB2_2D_T = nullptr;

void NF_B_Q_IB2_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_B_Q_IB2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *, double *)
{
}
