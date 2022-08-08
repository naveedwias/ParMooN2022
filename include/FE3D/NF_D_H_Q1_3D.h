/* for all functionals */
#include <stdlib.h>


static double NF_D_H_Q1_3D_Xi[]   = { 
         -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530 } ;
static double NF_D_H_Q1_3D_Eta[]  = {
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530  };
static double NF_D_H_Q1_3D_Zeta[] = {
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530  };

/* face 0                               0 */
static double *NF_D_H_Q1_3D_F0_Xi = nullptr;
static double *NF_D_H_Q1_3D_F0_Eta = nullptr;
static double *NF_D_H_Q1_3D_F0_Zeta = nullptr;

/* face 1                               1 */
static double *NF_D_H_Q1_3D_F1_Xi = nullptr;
static double *NF_D_H_Q1_3D_F1_Eta = nullptr;
static double *NF_D_H_Q1_3D_F1_Zeta = nullptr;

/* face 2                               2 */
static double *NF_D_H_Q1_3D_F2_Xi = nullptr;
static double *NF_D_H_Q1_3D_F2_Eta = nullptr;
static double *NF_D_H_Q1_3D_F2_Zeta = nullptr;

/* face 3                               3 */
static double *NF_D_H_Q1_3D_F3_Xi = nullptr;
static double *NF_D_H_Q1_3D_F3_Eta = nullptr;
static double *NF_D_H_Q1_3D_F3_Zeta = nullptr;

/* face 4                               4 */
static double *NF_D_H_Q1_3D_F4_Xi = nullptr;
static double *NF_D_H_Q1_3D_F4_Eta = nullptr;
static double *NF_D_H_Q1_3D_F4_Zeta = nullptr;

/* face 5                               5 */
static double *NF_D_H_Q1_3D_F5_Xi = nullptr;
static double *NF_D_H_Q1_3D_F5_Eta = nullptr;
static double *NF_D_H_Q1_3D_F5_Zeta = nullptr;

static double *NF_D_H_Q1_3D_XiArray[6] = { 
                        NF_D_H_Q1_3D_F0_Xi,
                        NF_D_H_Q1_3D_F1_Xi,
                        NF_D_H_Q1_3D_F2_Xi,
                        NF_D_H_Q1_3D_F3_Xi,
                        NF_D_H_Q1_3D_F4_Xi,
                        NF_D_H_Q1_3D_F5_Xi };

static double *NF_D_H_Q1_3D_EtaArray[6] = { 
                        NF_D_H_Q1_3D_F0_Eta,
                        NF_D_H_Q1_3D_F1_Eta,
                        NF_D_H_Q1_3D_F2_Eta,
                        NF_D_H_Q1_3D_F3_Eta,
                        NF_D_H_Q1_3D_F4_Eta,
                        NF_D_H_Q1_3D_F5_Eta };

static double *NF_D_H_Q1_3D_ZetaArray[6] = { 
                        NF_D_H_Q1_3D_F0_Zeta,
                        NF_D_H_Q1_3D_F1_Zeta,
                        NF_D_H_Q1_3D_F2_Zeta,
                        NF_D_H_Q1_3D_F3_Zeta,
                        NF_D_H_Q1_3D_F4_Zeta,
                        NF_D_H_Q1_3D_F5_Zeta };

static double *NF_D_H_Q1_3D_T = nullptr;
static double *NF_D_H_Q1_3D_S = nullptr;

void NF_D_H_Q1_3D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  Functionals[0] =
    ( +PointValues[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_Weight[26] ) * 0.125;
  Functionals[1] = 
    ( +PointValues[0]*NF_D_H_Q1_3D_Xi[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_Q1_3D_Xi[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_Q1_3D_Xi[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_Q1_3D_Xi[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_Q1_3D_Xi[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_Q1_3D_Xi[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_Q1_3D_Xi[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_Q1_3D_Xi[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_Q1_3D_Xi[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_Q1_3D_Xi[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Xi[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Xi[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Xi[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Xi[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Xi[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Xi[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Xi[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Xi[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Xi[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Xi[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Xi[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Xi[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Xi[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Xi[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Xi[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Xi[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Xi[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[2] = 
    ( +PointValues[0] *NF_D_H_Q1_3D_Eta[0] *NF_D_H_P2_3D_Weight[0]
      +PointValues[1] *NF_D_H_Q1_3D_Eta[1] *NF_D_H_P2_3D_Weight[1]
      +PointValues[2] *NF_D_H_Q1_3D_Eta[2] *NF_D_H_P2_3D_Weight[2]
      +PointValues[3] *NF_D_H_Q1_3D_Eta[3] *NF_D_H_P2_3D_Weight[3]
      +PointValues[4] *NF_D_H_Q1_3D_Eta[4] *NF_D_H_P2_3D_Weight[4]
      +PointValues[5] *NF_D_H_Q1_3D_Eta[5] *NF_D_H_P2_3D_Weight[5]
      +PointValues[6] *NF_D_H_Q1_3D_Eta[6] *NF_D_H_P2_3D_Weight[6]
      +PointValues[7] *NF_D_H_Q1_3D_Eta[7] *NF_D_H_P2_3D_Weight[7]
      +PointValues[8] *NF_D_H_Q1_3D_Eta[8] *NF_D_H_P2_3D_Weight[8]
      +PointValues[9] *NF_D_H_Q1_3D_Eta[9] *NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Eta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Eta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Eta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Eta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Eta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Eta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Eta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Eta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Eta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Eta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Eta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Eta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Eta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Eta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Eta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Eta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Eta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[3] = 
    ( +PointValues[0]*NF_D_H_Q1_3D_Zeta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_Q1_3D_Zeta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_Q1_3D_Zeta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_Q1_3D_Zeta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_Q1_3D_Zeta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_Q1_3D_Zeta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_Q1_3D_Zeta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_Q1_3D_Zeta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_Q1_3D_Zeta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_Q1_3D_Zeta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Zeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Zeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Zeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Zeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Zeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Zeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Zeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Zeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Zeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Zeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Zeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Zeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Zeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Zeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Zeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Zeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Zeta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[4] = 
    ( +PointValues[0] *NF_D_H_Q1_3D_Xi[0] *NF_D_H_Q1_3D_Eta[0] *NF_D_H_P2_3D_Weight[0]
      +PointValues[1] *NF_D_H_Q1_3D_Xi[1] *NF_D_H_Q1_3D_Eta[1] *NF_D_H_P2_3D_Weight[1]
      +PointValues[2] *NF_D_H_Q1_3D_Xi[2] *NF_D_H_Q1_3D_Eta[2] *NF_D_H_P2_3D_Weight[2]
      +PointValues[3] *NF_D_H_Q1_3D_Xi[3] *NF_D_H_Q1_3D_Eta[3] *NF_D_H_P2_3D_Weight[3]
      +PointValues[4] *NF_D_H_Q1_3D_Xi[4] *NF_D_H_Q1_3D_Eta[4] *NF_D_H_P2_3D_Weight[4]
      +PointValues[5] *NF_D_H_Q1_3D_Xi[5] *NF_D_H_Q1_3D_Eta[5] *NF_D_H_P2_3D_Weight[5]
      +PointValues[6] *NF_D_H_Q1_3D_Xi[6] *NF_D_H_Q1_3D_Eta[6] *NF_D_H_P2_3D_Weight[6]
      +PointValues[7] *NF_D_H_Q1_3D_Xi[7] *NF_D_H_Q1_3D_Eta[7] *NF_D_H_P2_3D_Weight[7]
      +PointValues[8] *NF_D_H_Q1_3D_Xi[8] *NF_D_H_Q1_3D_Eta[8] *NF_D_H_P2_3D_Weight[8]
      +PointValues[9] *NF_D_H_Q1_3D_Xi[9] *NF_D_H_Q1_3D_Eta[9] *NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Xi[10]*NF_D_H_Q1_3D_Eta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Xi[11]*NF_D_H_Q1_3D_Eta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Xi[12]*NF_D_H_Q1_3D_Eta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Xi[13]*NF_D_H_Q1_3D_Eta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Xi[14]*NF_D_H_Q1_3D_Eta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Xi[15]*NF_D_H_Q1_3D_Eta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Xi[16]*NF_D_H_Q1_3D_Eta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Xi[17]*NF_D_H_Q1_3D_Eta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Xi[18]*NF_D_H_Q1_3D_Eta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Xi[19]*NF_D_H_Q1_3D_Eta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Xi[20]*NF_D_H_Q1_3D_Eta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Xi[21]*NF_D_H_Q1_3D_Eta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Xi[22]*NF_D_H_Q1_3D_Eta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Xi[23]*NF_D_H_Q1_3D_Eta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Xi[24]*NF_D_H_Q1_3D_Eta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Xi[25]*NF_D_H_Q1_3D_Eta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Xi[26]*NF_D_H_Q1_3D_Eta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[5] = 
    ( +PointValues[0] *NF_D_H_Q1_3D_Xi[0] *NF_D_H_Q1_3D_Zeta[0] *NF_D_H_P2_3D_Weight[0]
      +PointValues[1] *NF_D_H_Q1_3D_Xi[1] *NF_D_H_Q1_3D_Zeta[1] *NF_D_H_P2_3D_Weight[1]
      +PointValues[2] *NF_D_H_Q1_3D_Xi[2] *NF_D_H_Q1_3D_Zeta[2] *NF_D_H_P2_3D_Weight[2]
      +PointValues[3] *NF_D_H_Q1_3D_Xi[3] *NF_D_H_Q1_3D_Zeta[3] *NF_D_H_P2_3D_Weight[3]
      +PointValues[4] *NF_D_H_Q1_3D_Xi[4] *NF_D_H_Q1_3D_Zeta[4] *NF_D_H_P2_3D_Weight[4]
      +PointValues[5] *NF_D_H_Q1_3D_Xi[5] *NF_D_H_Q1_3D_Zeta[5] *NF_D_H_P2_3D_Weight[5]
      +PointValues[6] *NF_D_H_Q1_3D_Xi[6] *NF_D_H_Q1_3D_Zeta[6] *NF_D_H_P2_3D_Weight[6]
      +PointValues[7] *NF_D_H_Q1_3D_Xi[7] *NF_D_H_Q1_3D_Zeta[7] *NF_D_H_P2_3D_Weight[7]
      +PointValues[8] *NF_D_H_Q1_3D_Xi[8] *NF_D_H_Q1_3D_Zeta[8] *NF_D_H_P2_3D_Weight[8]
      +PointValues[9] *NF_D_H_Q1_3D_Xi[9] *NF_D_H_Q1_3D_Zeta[9] *NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Xi[10]*NF_D_H_Q1_3D_Zeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Xi[11]*NF_D_H_Q1_3D_Zeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Xi[12]*NF_D_H_Q1_3D_Zeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Xi[13]*NF_D_H_Q1_3D_Zeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Xi[14]*NF_D_H_Q1_3D_Zeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Xi[15]*NF_D_H_Q1_3D_Zeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Xi[16]*NF_D_H_Q1_3D_Zeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Xi[17]*NF_D_H_Q1_3D_Zeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Xi[18]*NF_D_H_Q1_3D_Zeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Xi[19]*NF_D_H_Q1_3D_Zeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Xi[20]*NF_D_H_Q1_3D_Zeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Xi[21]*NF_D_H_Q1_3D_Zeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Xi[22]*NF_D_H_Q1_3D_Zeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Xi[23]*NF_D_H_Q1_3D_Zeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Xi[24]*NF_D_H_Q1_3D_Zeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Xi[25]*NF_D_H_Q1_3D_Zeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Xi[26]*NF_D_H_Q1_3D_Zeta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[6] = 
    ( +PointValues[0] *NF_D_H_Q1_3D_Eta[0] *NF_D_H_Q1_3D_Zeta[0] *NF_D_H_P2_3D_Weight[0]
      +PointValues[1] *NF_D_H_Q1_3D_Eta[1] *NF_D_H_Q1_3D_Zeta[1] *NF_D_H_P2_3D_Weight[1]
      +PointValues[2] *NF_D_H_Q1_3D_Eta[2] *NF_D_H_Q1_3D_Zeta[2] *NF_D_H_P2_3D_Weight[2]
      +PointValues[3] *NF_D_H_Q1_3D_Eta[3] *NF_D_H_Q1_3D_Zeta[3] *NF_D_H_P2_3D_Weight[3]
      +PointValues[4] *NF_D_H_Q1_3D_Eta[4] *NF_D_H_Q1_3D_Zeta[4] *NF_D_H_P2_3D_Weight[4]
      +PointValues[5] *NF_D_H_Q1_3D_Eta[5] *NF_D_H_Q1_3D_Zeta[5] *NF_D_H_P2_3D_Weight[5]
      +PointValues[6] *NF_D_H_Q1_3D_Eta[6] *NF_D_H_Q1_3D_Zeta[6] *NF_D_H_P2_3D_Weight[6]
      +PointValues[7] *NF_D_H_Q1_3D_Eta[7] *NF_D_H_Q1_3D_Zeta[7] *NF_D_H_P2_3D_Weight[7]
      +PointValues[8] *NF_D_H_Q1_3D_Eta[8] *NF_D_H_Q1_3D_Zeta[8] *NF_D_H_P2_3D_Weight[8]
      +PointValues[9] *NF_D_H_Q1_3D_Eta[9] *NF_D_H_Q1_3D_Zeta[9] *NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Eta[10]*NF_D_H_Q1_3D_Zeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Eta[11]*NF_D_H_Q1_3D_Zeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Eta[12]*NF_D_H_Q1_3D_Zeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Eta[13]*NF_D_H_Q1_3D_Zeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Eta[14]*NF_D_H_Q1_3D_Zeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Eta[15]*NF_D_H_Q1_3D_Zeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Eta[16]*NF_D_H_Q1_3D_Zeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Eta[17]*NF_D_H_Q1_3D_Zeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Eta[18]*NF_D_H_Q1_3D_Zeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Eta[19]*NF_D_H_Q1_3D_Zeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Eta[20]*NF_D_H_Q1_3D_Zeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Eta[21]*NF_D_H_Q1_3D_Zeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Eta[22]*NF_D_H_Q1_3D_Zeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Eta[23]*NF_D_H_Q1_3D_Zeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Eta[24]*NF_D_H_Q1_3D_Zeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Eta[25]*NF_D_H_Q1_3D_Zeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Eta[26]*NF_D_H_Q1_3D_Zeta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[7] = 
    ( +PointValues[0] *NF_D_H_Q1_3D_Xi[0] *NF_D_H_Q1_3D_Eta[0] *NF_D_H_Q1_3D_Zeta[0] *NF_D_H_P2_3D_Weight[0]
      +PointValues[1] *NF_D_H_Q1_3D_Xi[1] *NF_D_H_Q1_3D_Eta[1] *NF_D_H_Q1_3D_Zeta[1] *NF_D_H_P2_3D_Weight[1]
      +PointValues[2] *NF_D_H_Q1_3D_Xi[2] *NF_D_H_Q1_3D_Eta[2] *NF_D_H_Q1_3D_Zeta[2] *NF_D_H_P2_3D_Weight[2]
      +PointValues[3] *NF_D_H_Q1_3D_Xi[3] *NF_D_H_Q1_3D_Eta[3] *NF_D_H_Q1_3D_Zeta[3] *NF_D_H_P2_3D_Weight[3]
      +PointValues[4] *NF_D_H_Q1_3D_Xi[4] *NF_D_H_Q1_3D_Eta[4] *NF_D_H_Q1_3D_Zeta[4] *NF_D_H_P2_3D_Weight[4]
      +PointValues[5] *NF_D_H_Q1_3D_Xi[5] *NF_D_H_Q1_3D_Eta[5] *NF_D_H_Q1_3D_Zeta[5] *NF_D_H_P2_3D_Weight[5]
      +PointValues[6] *NF_D_H_Q1_3D_Xi[6] *NF_D_H_Q1_3D_Eta[6] *NF_D_H_Q1_3D_Zeta[6] *NF_D_H_P2_3D_Weight[6]
      +PointValues[7] *NF_D_H_Q1_3D_Xi[7] *NF_D_H_Q1_3D_Eta[7] *NF_D_H_Q1_3D_Zeta[7] *NF_D_H_P2_3D_Weight[7]
      +PointValues[8] *NF_D_H_Q1_3D_Xi[8] *NF_D_H_Q1_3D_Eta[8] *NF_D_H_Q1_3D_Zeta[8] *NF_D_H_P2_3D_Weight[8]
      +PointValues[9] *NF_D_H_Q1_3D_Xi[9] *NF_D_H_Q1_3D_Eta[9] *NF_D_H_Q1_3D_Zeta[9] *NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_Q1_3D_Xi[10]*NF_D_H_Q1_3D_Eta[10]*NF_D_H_Q1_3D_Zeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_Q1_3D_Xi[11]*NF_D_H_Q1_3D_Eta[11]*NF_D_H_Q1_3D_Zeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_Q1_3D_Xi[12]*NF_D_H_Q1_3D_Eta[12]*NF_D_H_Q1_3D_Zeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_Q1_3D_Xi[13]*NF_D_H_Q1_3D_Eta[13]*NF_D_H_Q1_3D_Zeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_Q1_3D_Xi[14]*NF_D_H_Q1_3D_Eta[14]*NF_D_H_Q1_3D_Zeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_Q1_3D_Xi[15]*NF_D_H_Q1_3D_Eta[15]*NF_D_H_Q1_3D_Zeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_Q1_3D_Xi[16]*NF_D_H_Q1_3D_Eta[16]*NF_D_H_Q1_3D_Zeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_Q1_3D_Xi[17]*NF_D_H_Q1_3D_Eta[17]*NF_D_H_Q1_3D_Zeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_Q1_3D_Xi[18]*NF_D_H_Q1_3D_Eta[18]*NF_D_H_Q1_3D_Zeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_Q1_3D_Xi[19]*NF_D_H_Q1_3D_Eta[19]*NF_D_H_Q1_3D_Zeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_Q1_3D_Xi[20]*NF_D_H_Q1_3D_Eta[20]*NF_D_H_Q1_3D_Zeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_Q1_3D_Xi[21]*NF_D_H_Q1_3D_Eta[21]*NF_D_H_Q1_3D_Zeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_Q1_3D_Xi[22]*NF_D_H_Q1_3D_Eta[22]*NF_D_H_Q1_3D_Zeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_Q1_3D_Xi[23]*NF_D_H_Q1_3D_Eta[23]*NF_D_H_Q1_3D_Zeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_Q1_3D_Xi[24]*NF_D_H_Q1_3D_Eta[24]*NF_D_H_Q1_3D_Zeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_Q1_3D_Xi[25]*NF_D_H_Q1_3D_Eta[25]*NF_D_H_Q1_3D_Zeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_Q1_3D_Xi[26]*NF_D_H_Q1_3D_Eta[26]*NF_D_H_Q1_3D_Zeta[26]*NF_D_H_P2_3D_Weight[26] ) * 3.375;
}

void NF_D_H_Q1_3D_EvalFace(const TCollection *, const TBaseCell *, int,
                           const double *, double *)
{
}

static int NF_D_H_Q1_3D_N_AllFunctionals = 8;
static int NF_D_H_Q1_3D_N_PointsAll = 27;
static int NF_D_H_Q1_3D_N_FaceFunctionals[] = { 0, 0, 0, 0, 0, 0 };
static int NF_D_H_Q1_3D_N_PointsFace[] = { 0, 0, 0, 0, 0, 0 };
