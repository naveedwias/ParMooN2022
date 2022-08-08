static double NF_C_T_UL1_2D_Xi[10] = { 0.0, 1.0, 0.0,
                  0.333333333333333333333333333333333,
                  0.797426985353087322398025276169754,
                  0.101286507323456338800987361915123,
                  0.101286507323456338800987361915123,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447, 
                  0.470142064105115089770441209513447 };

static double NF_C_T_UL1_2D_Eta[] = { 0.0, 0.0, 1.0,
                  0.333333333333333333333333333333333, 
                  0.101286507323456338800987361915123,
                  0.797426985353087322398025276169754, 
                  0.101286507323456338800987361915123, 
                  0.470142064105115089770441209513447,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447 }; 

static double NF_C_T_UL1_2D_W[7] =
                { 0.1125, 
                  0.0629695902724135762978419727500906,
                  0.0629695902724135762978419727500906,
                  0.0629695902724135762978419727500906,
                  0.0661970763942530903688246939165759,
                  0.0661970763942530903688246939165759,
                  0.0661970763942530903688246939165759 };

static double NF_C_T_UL1_2D_T[2] = { -1, 1 };

void NF_C_T_UL1_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = 2*(  PointValues[3]*NF_C_T_UL1_2D_W[0]
                      + PointValues[4]*NF_C_T_UL1_2D_W[1]
                      + PointValues[5]*NF_C_T_UL1_2D_W[2]
                      + PointValues[6]*NF_C_T_UL1_2D_W[3]
                      + PointValues[7]*NF_C_T_UL1_2D_W[4]
                      + PointValues[8]*NF_C_T_UL1_2D_W[5]
                      + PointValues[9]*NF_C_T_UL1_2D_W[6] );
}

void NF_C_T_UL1_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}
