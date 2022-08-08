static double NF_D_Q_P2_2D_Xi[9]={ 0.774596669241483, -0.774596669241483,
                                 0.774596669241483, -0.774596669241483,
                                 0.774596669241483, -0.774596669241483,
                                 0,                  0,
                                 0 };
static double NF_D_Q_P2_2D_Eta[9]={  0.774596669241483,  0.774596669241483,
                                  -0.774596669241483, -0.774596669241483,
                                   0,                  0,
                                   0.774596669241483, -0.774596669241483,
                                   0 };
// 3*xi**2-1
static double NF_D_Q_P2_2D_3Xi2_1[9]={  0.8,  0.8,
                                      0.8,  0.8,
                                      0.8,  0.8,
                                     -1.0, -1.0,
                                     -1.0 };
// xi*eta
static double NF_D_Q_P2_2D_XiEta[9]={  0.6, -0.6,
                                    -0.6,  0.6,
                                     0.0,  0.0,
                                     0.0,  0.0,
                                     0.0 };
// 3*eta**2-1
static double NF_D_Q_P2_2D_3Eta2_1[9]={  0.8,  0.8,
                                       0.8,  0.8,
                                      -1.0, -1.0,
                                       0.8, 0.8,
                                      -1.0 };

static double *NF_D_Q_P2_2D_t = nullptr;

void NF_D_Q_P2_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  static double weights[9]={ 0.308641975308642, 0.308641975308642,
                             0.308641975308642, 0.308641975308642,
                             0.493827160493827, 0.493827160493827,
                             0.493827160493827, 0.493827160493827,
                             0.790123456790123 };
  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2]
                   +weights[3]*PointValues[3]
                   +weights[4]*PointValues[4]
                   +weights[5]*PointValues[5]
                   +weights[6]*PointValues[6]
                   +weights[7]*PointValues[7]
                   +weights[8]*PointValues[8];
  Functionals[1] =  weights[0]*PointValues[0]*NF_D_Q_P2_2D_Xi[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P2_2D_Xi[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P2_2D_Xi[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P2_2D_Xi[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P2_2D_Xi[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P2_2D_Xi[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P2_2D_Xi[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P2_2D_Xi[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P2_2D_Xi[8];
  Functionals[2] =  weights[0]*PointValues[0]*NF_D_Q_P2_2D_Eta[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P2_2D_Eta[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P2_2D_Eta[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P2_2D_Eta[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P2_2D_Eta[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P2_2D_Eta[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P2_2D_Eta[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P2_2D_Eta[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P2_2D_Eta[8];
  Functionals[3] =  weights[0]*PointValues[0]*NF_D_Q_P2_2D_3Xi2_1[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P2_2D_3Xi2_1[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P2_2D_3Xi2_1[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P2_2D_3Xi2_1[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P2_2D_3Xi2_1[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P2_2D_3Xi2_1[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P2_2D_3Xi2_1[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P2_2D_3Xi2_1[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P2_2D_3Xi2_1[8];
  Functionals[4] =  weights[0]*PointValues[0]*NF_D_Q_P2_2D_XiEta[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P2_2D_XiEta[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P2_2D_XiEta[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P2_2D_XiEta[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P2_2D_XiEta[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P2_2D_XiEta[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P2_2D_XiEta[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P2_2D_XiEta[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P2_2D_XiEta[8];
  Functionals[5] =  weights[0]*PointValues[0]*NF_D_Q_P2_2D_3Eta2_1[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P2_2D_3Eta2_1[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P2_2D_3Eta2_1[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P2_2D_3Eta2_1[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P2_2D_3Eta2_1[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P2_2D_3Eta2_1[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P2_2D_3Eta2_1[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P2_2D_3Eta2_1[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P2_2D_3Eta2_1[8];

  Functionals[0] *= 0.25;
  Functionals[1] *= 0.25;
  Functionals[2] *= 0.25;
  Functionals[3] *= 0.25;
  Functionals[4] *= 0.25;
  Functionals[5] *= 0.25;
}

void NF_D_Q_P2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}
