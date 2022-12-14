static double NF_D_Q_Q2_2D_Xi[9] = {
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704 };
                                 
static double NF_D_Q_Q2_2D_Eta[9] = { 
   -0.77459666924148337704, -0.77459666924148337704, -0.77459666924148337704,
    0, 0, 0,
    0.77459666924148337704, 0.77459666924148337704, 0.77459666924148337704 };

static double NF_D_Q_Q2_2D_Weight0[9] = {
   0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827,
   0.12345679012345679012, 0.1975308641975308642, 0.12345679012345679012,
   0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827 };

static double NF_D_Q_Q2_2D_Weight1[9] = {
  -0.17930478454663967061, 0, 0.17930478454663967061,
  -0.28688765527462347298, 0, 0.28688765527462347298,
  -0.17930478454663967061, 0, 0.17930478454663967061 };

static double NF_D_Q_Q2_2D_Weight2[9] = {
   0.15432098765432098765, -0.30864197530864197531, 0.15432098765432098765,
   0.24691358024691358025, -0.49382716049382716049, 0.24691358024691358025,
   0.15432098765432098765, -0.30864197530864197531, 0.15432098765432098765 };

static double NF_D_Q_Q2_2D_Weight3[9] = {
  -0.17930478454663967061, -0.28688765527462347298, -0.17930478454663967061,
   0, 0, 0,
   0.17930478454663967061, 0.28688765527462347298, 0.17930478454663967061 };

static double NF_D_Q_Q2_2D_Weight4[9] = {
   0.41666666666666666667, 0, -0.41666666666666666667,
   0, 0, 0,
  -0.41666666666666666667, 0, 0.41666666666666666667 };

static double NF_D_Q_Q2_2D_Weight5[9] = {
  -0.35860956909327934122, 0.71721913818655868244, -0.35860956909327934122,
   0, 0, 0,
   0.35860956909327934122, -0.71721913818655868244, 0.35860956909327934122 };

static double NF_D_Q_Q2_2D_Weight6[9] = {
  0.15432098765432098765, 0.24691358024691358025, 0.15432098765432098765,
 -0.30864197530864197531, -0.49382716049382716049, -0.30864197530864197531,
  0.15432098765432098765, 0.24691358024691358025, 0.15432098765432098765 };

static double NF_D_Q_Q2_2D_Weight7[9] = {
 -0.35860956909327934122, 0, 0.35860956909327934122,
  0.71721913818655868244, 0, -0.71721913818655868244,
 -0.35860956909327934122, 0, 0.35860956909327934122 };

static double NF_D_Q_Q2_2D_Weight8[9] = {
  0.30864197530864197531, -0.61728395061728395062, 0.30864197530864197531,
 -0.61728395061728395062, 1.2345679012345679012, -0.61728395061728395062,
  0.30864197530864197531, -0.61728395061728395062, 0.30864197530864197531 };

static double *NF_D_Q_Q2_2D_t = nullptr;

void NF_D_Q_Q2_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] =  NF_D_Q_Q2_2D_Weight0[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight0[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight0[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight0[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight0[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight0[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight0[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight0[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight0[8]*PointValues[8];

  Functionals[1] =  NF_D_Q_Q2_2D_Weight1[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight1[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight1[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight1[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight1[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight1[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight1[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight1[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight1[8]*PointValues[8];

  Functionals[2] =  NF_D_Q_Q2_2D_Weight2[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight2[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight2[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight2[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight2[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight2[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight2[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight2[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight2[8]*PointValues[8];

  Functionals[3] =  NF_D_Q_Q2_2D_Weight3[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight3[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight3[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight3[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight3[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight3[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight3[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight3[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight3[8]*PointValues[8];

  Functionals[4] =  NF_D_Q_Q2_2D_Weight4[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight4[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight4[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight4[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight4[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight4[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight4[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight4[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight4[8]*PointValues[8];

  Functionals[5] =  NF_D_Q_Q2_2D_Weight5[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight5[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight5[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight5[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight5[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight5[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight5[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight5[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight5[8]*PointValues[8];

  Functionals[6] =  NF_D_Q_Q2_2D_Weight6[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight6[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight6[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight6[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight6[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight6[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight6[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight6[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight6[8]*PointValues[8];

  Functionals[7] =  NF_D_Q_Q2_2D_Weight7[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight7[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight7[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight7[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight7[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight7[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight7[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight7[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight7[8]*PointValues[8];

  Functionals[8] =  NF_D_Q_Q2_2D_Weight8[0]*PointValues[0]
                   +NF_D_Q_Q2_2D_Weight8[1]*PointValues[1]
                   +NF_D_Q_Q2_2D_Weight8[2]*PointValues[2]
                   +NF_D_Q_Q2_2D_Weight8[3]*PointValues[3]
                   +NF_D_Q_Q2_2D_Weight8[4]*PointValues[4]
                   +NF_D_Q_Q2_2D_Weight8[5]*PointValues[5]
                   +NF_D_Q_Q2_2D_Weight8[6]*PointValues[6]
                   +NF_D_Q_Q2_2D_Weight8[7]*PointValues[7]
                   +NF_D_Q_Q2_2D_Weight8[8]*PointValues[8];
}

void NF_D_Q_Q2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}
