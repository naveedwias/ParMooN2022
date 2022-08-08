static double NF_N_T_P1_2D_Xi[] = 
        { 0.11270166537925831149, 0.5, 
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149,
          0, 0, 0 };
static double NF_N_T_P1_2D_Eta[] = 
        { 0, 0, 0,
          0.11270166537925831149, 0.5,
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149 };
static double NF_N_T_P1_2D_T[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

void NF_N_T_P1_2D_EvalAll(const TCollection *, const TBaseCell *,
                          const double *PointValues, double *Functionals)
{
  static double weights[]={ 0.277777777777777777777778,
                            0.444444444444444444444444,
                            0.277777777777777777777778 };
  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2];
  Functionals[1] =  weights[0]*PointValues[3]
                   +weights[1]*PointValues[4]
                   +weights[2]*PointValues[5];
  Functionals[2] =  weights[0]*PointValues[6]
                   +weights[1]*PointValues[7]
                   +weights[2]*PointValues[8];
}

void NF_N_T_P1_2D_EvalEdge(const TCollection *, const TBaseCell *, int,
                           const double *PointValues, double *Functionals)
{
  static double weights[3] = { 0.5555555555555555555555555555555556,
                               0.88888888888888888888888888888888889,
                               0.5555555555555555555555555555555556 };
  Functionals[0] =(  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2])*0.5;
}
