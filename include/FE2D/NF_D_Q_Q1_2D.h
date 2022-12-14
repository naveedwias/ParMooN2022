static double NF_D_Q_Q1_2D_Xi[9] = {
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704 };

static double NF_D_Q_Q1_2D_Eta[9] = {
   -0.77459666924148337704, -0.77459666924148337704, -0.77459666924148337704,
    0, 0, 0,
    0.77459666924148337704, 0.77459666924148337704, 0.77459666924148337704 };

static double *NF_D_Q_Q1_2D_t = nullptr;

static double NF_D_Q_Q1_2D_Weight0[9] = {
    0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827,
    0.12345679012345679012, 0.1975308641975308642, 0.12345679012345679012,
    0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827 };

static double NF_D_Q_Q1_2D_Weight1[9] = {
   -0.17930478454663967061, 0, 0.17930478454663967061,
   -0.28688765527462347298, 0, 0.28688765527462347298,
   -0.17930478454663967061, 0, 0.17930478454663967061 };

static double NF_D_Q_Q1_2D_Weight2[9] = {
   -0.17930478454663967061, -0.28688765527462347298, -0.17930478454663967061,
    0, 0, 0,
    0.17930478454663967061, 0.28688765527462347298, 0.17930478454663967061 };

static double NF_D_Q_Q1_2D_Weight3[9] = {
    0.41666666666666666667, 0, -0.41666666666666666667,
    0, 0, 0,
   -0.41666666666666666667, 0, 0.41666666666666666667 };

void NF_D_Q_Q1_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] =  NF_D_Q_Q1_2D_Weight0[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight0[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight0[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight0[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight0[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight0[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight0[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight0[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight0[8]*PointValues[8];

  Functionals[1] =  NF_D_Q_Q1_2D_Weight1[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight1[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight1[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight1[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight1[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight1[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight1[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight1[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight1[8]*PointValues[8];

  Functionals[2] =  NF_D_Q_Q1_2D_Weight2[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight2[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight2[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight2[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight2[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight2[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight2[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight2[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight2[8]*PointValues[8];

  Functionals[3] =  NF_D_Q_Q1_2D_Weight3[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight3[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight3[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight3[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight3[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight3[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight3[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight3[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight3[8]*PointValues[8];
}

void NF_D_Q_Q1_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}
