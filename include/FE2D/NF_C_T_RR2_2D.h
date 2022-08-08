// red refined triangle, second order, used for LPS
static double NF_C_T_RR2_2D_Xi[] = { 0, 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25,
                                     0, 0, 0, 0, 0.25, 0.5, 0.25 };

static double NF_C_T_RR2_2D_Eta[] = { 0, 0, 0, 0, 0, 0.25, 0.5, 0.75,
                                      1, 0.75, 0.5, 0.25, 0.25, 0.25, 0.5 };

static double NF_C_T_RR2_2D_T[] = { -1, -0.5, 0, 0.5, 1 };

void NF_C_T_RR2_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 15*SizeOfDouble);
}

void NF_C_T_RR2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 5*SizeOfDouble);
}

// TNodalFunctional2D *NF_C_T_RR2_2D_Obj = new TNodalFunctional2D
//         (NF_C_T_B2_2D, 15, 5, 15, 5, NF_C_T_RR2_2D_Xi, NF_C_T_RR2_2D_Eta,
//          NF_C_T_RR2_2D_T, NF_C_T_RR2_2D_EvalAll, NF_C_T_RR2_2D_EvalEdge);
