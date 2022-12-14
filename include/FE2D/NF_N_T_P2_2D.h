static double NF_N_T_P2_2D_Xi[16] = 
        { 
          0.112701665379258311482073460022, 0.5,
          0.887298334620741688517926539978,
          0.887298334620741688517926539978, 0.5,
          0.112701665379258311482073460022,
          0, 0, 0,
          0.333333333333333333333333333333,
          0.5, 0.5, 0.0,
          0.0, 1.0, 0.0
        };

static double NF_N_T_P2_2D_Eta[16] = 
        { 
          0, 0, 0,
          0.112701665379258311482073460022, 0.5,
          0.887298334620741688517926539978,
          0.887298334620741688517926539978, 0.5,
          0.112701665379258311482073460022,
          0.333333333333333333333333333333,
          0.0, 0.5, 0.5,
          0.0, 0.0, 1.0
        };

static double NF_N_T_P2_2D_T[3] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

static double NF_N_T_P2_2D_EdgeWeight0[3] = 
        { 0.277777777777777777777777777778,
          0.444444444444444444444444444444,
          0.277777777777777777777777777778 };

static double NF_N_T_P2_2D_EdgeWeight1[3] = 
        { -0.64549722436790281419654423329706660,
           0.0,
           0.64549722436790281419654423329706660 };

static double NF_N_T_P2_2D_CellWeight0[7] =
        { 0.45,
          0.1333333333333333333333333333333333,
          0.1333333333333333333333333333333333,
          0.1333333333333333333333333333333333,
          0.05, 0.05, 0.05 };

void NF_N_T_P2_2D_EvalAll(const TCollection *Coll, const TBaseCell *Cell,
                          const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  const TBaseCell *neigh;

  Functionals[0] = ( NF_N_T_P2_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_T_P2_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_T_P2_2D_EdgeWeight0[2]*PointValues[2]);
  Functionals[1] = ( NF_N_T_P2_2D_EdgeWeight0[0]*PointValues[3]
                    +NF_N_T_P2_2D_EdgeWeight0[1]*PointValues[4]
                    +NF_N_T_P2_2D_EdgeWeight0[2]*PointValues[5]);
  Functionals[2] = ( NF_N_T_P2_2D_EdgeWeight0[0]*PointValues[6]
                    +NF_N_T_P2_2D_EdgeWeight0[1]*PointValues[7]
                    +NF_N_T_P2_2D_EdgeWeight0[2]*PointValues[8]);

  Functionals[3] = ( NF_N_T_P2_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_T_P2_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_T_P2_2D_EdgeWeight1[2]*PointValues[2]);
  Functionals[4] = ( NF_N_T_P2_2D_EdgeWeight1[0]*PointValues[3]
                    +NF_N_T_P2_2D_EdgeWeight1[1]*PointValues[4]
                    +NF_N_T_P2_2D_EdgeWeight1[2]*PointValues[5]);
  Functionals[5] = ( NF_N_T_P2_2D_EdgeWeight1[0]*PointValues[6]
                    +NF_N_T_P2_2D_EdgeWeight1[1]*PointValues[7]
                    +NF_N_T_P2_2D_EdgeWeight1[2]*PointValues[8]);

  Functionals[6] =( NF_N_T_P2_2D_CellWeight0[0]*PointValues[ 9]
                   +NF_N_T_P2_2D_CellWeight0[1]*PointValues[10]
                   +NF_N_T_P2_2D_CellWeight0[2]*PointValues[11]
                   +NF_N_T_P2_2D_CellWeight0[3]*PointValues[12]
                   +NF_N_T_P2_2D_CellWeight0[4]*PointValues[13]
                   +NF_N_T_P2_2D_CellWeight0[5]*PointValues[14]
                   +NF_N_T_P2_2D_CellWeight0[6]*PointValues[15] );

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[3] = -Functionals[3];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[4] = -Functionals[4];
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
        Functionals[3] = -Functionals[3];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 4] = -Functionals[ 4];
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

void NF_N_T_P2_2D_EvalEdge(const TCollection *Coll, const TBaseCell *Cell, int Joint,
                           const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  const TBaseCell *neigh;

  Functionals[0] = ( NF_N_T_P2_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_T_P2_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_T_P2_2D_EdgeWeight0[2]*PointValues[2]);

  Functionals[1] = ( NF_N_T_P2_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_T_P2_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_T_P2_2D_EdgeWeight1[2]*PointValues[2]);

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
