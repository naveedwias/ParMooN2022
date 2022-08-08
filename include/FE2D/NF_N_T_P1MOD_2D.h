static double NF_N_T_P1MOD_2D_Xi[] = 
        { 0.11270166537925831149, 0.5, 
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149,
          0, 0, 0 };

static double NF_N_T_P1MOD_2D_Eta[] = 
        { 0, 0, 0,
          0.11270166537925831149, 0.5,
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149 };

static double NF_N_T_P1MOD_2D_T[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

void NF_N_T_P1MOD_2D_EvalAll(const TCollection *Coll, const TBaseCell *Cell,
                             const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  const TBaseCell *neigh;

  static double weights[]={ 0.27777777777777777777777777778,
                            0.44444444444444444444444444444,
                            0.27777777777777777777777777778 };

  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2];
  Functionals[2] =  weights[0]*PointValues[3]
                   +weights[1]*PointValues[4]
                   +weights[2]*PointValues[5];
  Functionals[4] =  weights[0]*PointValues[6]
                   +weights[1]*PointValues[7]
                   +weights[2]*PointValues[8];

  Functionals[1] = 60*( weights[0]*(NF_N_T_P1MOD_2D_Xi[0]-0.5)*PointValues[0]
                       +weights[1]*(NF_N_T_P1MOD_2D_Xi[1]-0.5)*PointValues[1]
                       +weights[2]*(NF_N_T_P1MOD_2D_Xi[2]-0.5)*PointValues[2]);
  Functionals[3] = 60*( weights[0]*(NF_N_T_P1MOD_2D_Eta[3]-0.5)*PointValues[3]
                       +weights[1]*(NF_N_T_P1MOD_2D_Eta[4]-0.5)*PointValues[4]
                       +weights[2]*(NF_N_T_P1MOD_2D_Eta[5]-0.5)*PointValues[5]);
  Functionals[5] = 60*( weights[0]*(0.5-NF_N_T_P1MOD_2D_Eta[6])*PointValues[6]
                       +weights[1]*(0.5-NF_N_T_P1MOD_2D_Eta[7])*PointValues[7]
                       +weights[2]*(0.5-NF_N_T_P1MOD_2D_Eta[8])*PointValues[8]);
  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[1] = -Functionals[1];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[3] = -Functionals[3];
    if(Cell->GetVertex(2) > Cell->GetVertex(0))
      Functionals[5] = -Functionals[5];
  }
  */

  if(Cell)
  {
    OwnNum = Coll->get_cell_index(Cell);

    neigh = Cell->GetJoint(0)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[1] = -Functionals[1];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[3] = -Functionals[3];
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 5] = -Functionals[ 5];
    } // endif neigh
  } // endif Cell
}

void NF_N_T_P1MOD_2D_EvalEdge(const TCollection *Coll, const TBaseCell *Cell, int Joint,
                              const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  const TBaseCell *neigh;

  static double weights[3]={ 0.27777777777777777777777777778,
                            0.44444444444444444444444444444,
                            0.27777777777777777777777777778 };

  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2];

  Functionals[1] = 60*( weights[0]*NF_N_T_P1MOD_2D_T[0]*PointValues[0]
                       +weights[1]*NF_N_T_P1MOD_2D_T[1]*PointValues[1]
                       +weights[2]*NF_N_T_P1MOD_2D_T[2]*PointValues[2])*0.5;

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%3))
    //   Functionals[1] = -Functionals[1];
    // /*
    neigh = Cell->GetJoint(Joint)->GetNeighbour(Cell);
    if(neigh)
    {
      OwnNum = Coll->get_cell_index(Cell);
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[1] = -Functionals[1];
    } // endif neigh
    // */
  }
}
