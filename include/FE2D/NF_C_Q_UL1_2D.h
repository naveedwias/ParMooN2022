static double NF_C_Q_UL1_2D_Xi[] = { -1, 1, 1, -1,
 -.5773502691896257645091486, .5773502691896257645091486,
 -.5773502691896257645091486, .5773502691896257645091486  };

static double NF_C_Q_UL1_2D_Eta[] = { -1, -1, 1, 1,
  -.5773502691896257645091486, -.5773502691896257645091486,
   .5773502691896257645091486,  .5773502691896257645091486 };

static double NF_C_Q_UL1_2D_T[] = { -1, 1 };

void NF_C_Q_UL1_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4]+PointValues[5]+PointValues[6]+PointValues[7];
}

void NF_C_Q_UL1_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}
