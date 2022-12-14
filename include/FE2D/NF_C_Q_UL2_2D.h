static double NF_C_Q_UL2_2D_Xi[] = {
-1, 0, 1, 1, 1, 0, -1, -1, -.7745966692414833770358530, 0., .7745966692414833770358530,
-.7745966692414833770358530, 0., .7745966692414833770358530, -.7745966692414833770358530, 0.,
.7745966692414833770358530
};
static double NF_C_Q_UL2_2D_Eta[] = {
-1, -1, -1, 0, 1, 1, 1, 0, -.7745966692414833770358530, -.7745966692414833770358530,
-.7745966692414833770358530, 0., 0., 0., .7745966692414833770358530, .7745966692414833770358530,
.7745966692414833770358530
};
static double NF_C_Q_UL2_2D_T[] = { -1, 0, 1 };

static double NF_C_Q_UL2_2D_W8[] = {
.3086419753086419753086424, .4938271604938271604938270, .3086419753086419753086424,
.4938271604938271604938270, .7901234567901234567901219, .4938271604938271604938270,
.3086419753086419753086424, .4938271604938271604938270, .3086419753086419753086424
};

static double NF_C_Q_UL2_2D_W9[] = {
-.2390730460621862274802019, -.3825168736994979639683223, -.2390730460621862274802019, 0., 0.,
0., .2390730460621862274802019, .3825168736994979639683223, .2390730460621862274802019
};

static double NF_C_Q_UL2_2D_W10[] = {
-.2390730460621862274802019, 0., .2390730460621862274802019, -.3825168736994979639683223, 0.,
.3825168736994979639683223, -.2390730460621862274802019, 0., .2390730460621862274802019
};

void NF_C_Q_UL2_2D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];

  Functionals[8] = NF_C_Q_UL2_2D_W8[0]*PointValues[ 8]
                  +NF_C_Q_UL2_2D_W8[1]*PointValues[ 9]
                  +NF_C_Q_UL2_2D_W8[2]*PointValues[10]
                  +NF_C_Q_UL2_2D_W8[3]*PointValues[11]
                  +NF_C_Q_UL2_2D_W8[4]*PointValues[12]
                  +NF_C_Q_UL2_2D_W8[5]*PointValues[13]
                  +NF_C_Q_UL2_2D_W8[6]*PointValues[14]
                  +NF_C_Q_UL2_2D_W8[7]*PointValues[15]
                  +NF_C_Q_UL2_2D_W8[8]*PointValues[16];

  Functionals[9] = NF_C_Q_UL2_2D_W9[0]*PointValues[ 8]
                  +NF_C_Q_UL2_2D_W9[1]*PointValues[ 9]
                  +NF_C_Q_UL2_2D_W9[2]*PointValues[10]
                  +NF_C_Q_UL2_2D_W9[3]*PointValues[11]
                  +NF_C_Q_UL2_2D_W9[4]*PointValues[12]
                  +NF_C_Q_UL2_2D_W9[5]*PointValues[13]
                  +NF_C_Q_UL2_2D_W9[6]*PointValues[14]
                  +NF_C_Q_UL2_2D_W9[7]*PointValues[15]
                  +NF_C_Q_UL2_2D_W9[8]*PointValues[16];

  Functionals[10] = NF_C_Q_UL2_2D_W10[0]*PointValues[ 8]
                   +NF_C_Q_UL2_2D_W10[1]*PointValues[ 9]
                   +NF_C_Q_UL2_2D_W10[2]*PointValues[10]
                   +NF_C_Q_UL2_2D_W10[3]*PointValues[11]
                   +NF_C_Q_UL2_2D_W10[4]*PointValues[12]
                   +NF_C_Q_UL2_2D_W10[5]*PointValues[13]
                   +NF_C_Q_UL2_2D_W10[6]*PointValues[14]
                   +NF_C_Q_UL2_2D_W10[7]*PointValues[15]
                   +NF_C_Q_UL2_2D_W10[8]*PointValues[16];
}

void NF_C_Q_UL2_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                            const double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}
