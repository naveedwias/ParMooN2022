static double NF_C_T_UL2_2D_Xi[13] = {
                  0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 
                  0.333333333333333333333333333333333,
                  0.797426985353087322398025276169754,
                  0.101286507323456338800987361915123,
                  0.101286507323456338800987361915123,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447, 
                  0.470142064105115089770441209513447 };

static double NF_C_T_UL2_2D_Eta[] = {
                  0.0, 0.0, 0.0, 0.5, 1.0, 0.5,
                  0.333333333333333333333333333333333, 
                  0.101286507323456338800987361915123,
                  0.797426985353087322398025276169754, 
                  0.101286507323456338800987361915123, 
                  0.470142064105115089770441209513447,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447 }; 

static double NF_C_T_UL2_2D_T[3] = { -1, 0, 1 };

static double NF_C_T_UL2_2D_W0[7] = {
                     0.2250,
                     0.12593918054482715260,
                     0.12593918054482715260,
                     0.12593918054482715260,
                     0.13239415278850618074,
                     0.13239415278850618074,
                     0.13239415278850618074 };

static double NF_C_T_UL2_2D_W1[7] = {
                     0.0,
                     0.70137089077708780990,
                     -.35068544538854390492,
                     -.35068544538854390492,
                     -.43470422411042114317,
                     0.21735211205521057155,
                     0.21735211205521057155 };

static double NF_C_T_UL2_2D_W2[7] = {
                     0.0,
                     -.35068544538854390492,
                     0.70137089077708780990,
                     -.35068544538854390492,
                     0.21735211205521057155,
                     -.43470422411042114317,
                     0.21735211205521057155 };

void NF_C_T_UL2_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];

  Functionals[6] =      PointValues[ 6]*NF_C_T_UL2_2D_W0[0]
                      + PointValues[ 7]*NF_C_T_UL2_2D_W0[1]
                      + PointValues[ 8]*NF_C_T_UL2_2D_W0[2]
                      + PointValues[ 9]*NF_C_T_UL2_2D_W0[3]
                      + PointValues[10]*NF_C_T_UL2_2D_W0[4]
                      + PointValues[11]*NF_C_T_UL2_2D_W0[5]
                      + PointValues[12]*NF_C_T_UL2_2D_W0[6] ;
  Functionals[7] =      PointValues[ 6]*NF_C_T_UL2_2D_W1[0]
                      + PointValues[ 7]*NF_C_T_UL2_2D_W1[1]
                      + PointValues[ 8]*NF_C_T_UL2_2D_W1[2]
                      + PointValues[ 9]*NF_C_T_UL2_2D_W1[3]
                      + PointValues[10]*NF_C_T_UL2_2D_W1[4]
                      + PointValues[11]*NF_C_T_UL2_2D_W1[5]
                      + PointValues[12]*NF_C_T_UL2_2D_W1[6] ;
  Functionals[8] =      PointValues[ 6]*NF_C_T_UL2_2D_W2[0]
                      + PointValues[ 7]*NF_C_T_UL2_2D_W2[1]
                      + PointValues[ 8]*NF_C_T_UL2_2D_W2[2]
                      + PointValues[ 9]*NF_C_T_UL2_2D_W2[3]
                      + PointValues[10]*NF_C_T_UL2_2D_W2[4]
                      + PointValues[11]*NF_C_T_UL2_2D_W2[5]
                      + PointValues[12]*NF_C_T_UL2_2D_W2[6] ;
}

void NF_C_T_UL2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}
