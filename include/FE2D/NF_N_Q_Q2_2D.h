static double NF_N_Q_Q2_2D_Xi[21] = 
        { -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992, 
           1, 1, 1, 
           0.77459666924148337703585307995647992, 0, 
          -0.77459666924148337703585307995647992, 
          -1, -1, -1,
          -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992,
          -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992,
          -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992
         };

static double NF_N_Q_Q2_2D_Eta[21] = 
        { -1, -1, -1, 
          -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992,
           1, 1, 1,
           0.77459666924148337703585307995647992, 0,
          -0.77459666924148337703585307995647992,

          -0.77459666924148337703585307995647992,
          -0.77459666924148337703585307995647992,
          -0.77459666924148337703585307995647992,
           0, 0, 0,
           0.77459666924148337703585307995647992,
           0.77459666924148337703585307995647992,
           0.77459666924148337703585307995647992
 };

static double NF_N_Q_Q2_2D_T[3] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

static double NF_N_Q_Q2_2D_EdgeWeight0[3] = 
        { 0.277777777777777777777777777778,
          0.444444444444444444444444444444,
          0.277777777777777777777777777778 };

static double NF_N_Q_Q2_2D_EdgeWeight1[3] = 
        { -0.64549722436790281419654423329706660,
           0.0,
           0.64549722436790281419654423329706660 };

static double NF_N_Q_Q2_2D_CellWeight0[9] = {
           0.07716049382716049382716049, 0.1234567901234567901234568,
           0.07716049382716049382716049, 0.1234567901234567901234568,
           0.1975308641975308641975309, 0.1234567901234567901234568,
           0.07716049382716049382716049, 0.1234567901234567901234568,
           0.07716049382716049382716049 };

void NF_N_Q_Q2_2D_EvalAll(const TCollection *Coll, const TBaseCell *Cell,
                          const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_Q_Q2_2D_EdgeWeight0[0]*PointValues[ 0]
                    +NF_N_Q_Q2_2D_EdgeWeight0[1]*PointValues[ 1]
                    +NF_N_Q_Q2_2D_EdgeWeight0[2]*PointValues[ 2]);
  Functionals[1] = ( NF_N_Q_Q2_2D_EdgeWeight0[0]*PointValues[ 3]
                    +NF_N_Q_Q2_2D_EdgeWeight0[1]*PointValues[ 4]
                    +NF_N_Q_Q2_2D_EdgeWeight0[2]*PointValues[ 5]);
  Functionals[2] = ( NF_N_Q_Q2_2D_EdgeWeight0[0]*PointValues[ 6]
                    +NF_N_Q_Q2_2D_EdgeWeight0[1]*PointValues[ 7]
                    +NF_N_Q_Q2_2D_EdgeWeight0[2]*PointValues[ 8]);
  Functionals[3] = ( NF_N_Q_Q2_2D_EdgeWeight0[0]*PointValues[ 9]
                    +NF_N_Q_Q2_2D_EdgeWeight0[1]*PointValues[10]
                    +NF_N_Q_Q2_2D_EdgeWeight0[2]*PointValues[11]);

  Functionals[4] = ( NF_N_Q_Q2_2D_EdgeWeight1[0]*PointValues[ 0]
                    +NF_N_Q_Q2_2D_EdgeWeight1[1]*PointValues[ 1]
                    +NF_N_Q_Q2_2D_EdgeWeight1[2]*PointValues[ 2]);
  Functionals[5] = ( NF_N_Q_Q2_2D_EdgeWeight1[0]*PointValues[ 3]
                    +NF_N_Q_Q2_2D_EdgeWeight1[1]*PointValues[ 4]
                    +NF_N_Q_Q2_2D_EdgeWeight1[2]*PointValues[ 5]);
  Functionals[6] = ( NF_N_Q_Q2_2D_EdgeWeight1[0]*PointValues[ 6]
                    +NF_N_Q_Q2_2D_EdgeWeight1[1]*PointValues[ 7]
                    +NF_N_Q_Q2_2D_EdgeWeight1[2]*PointValues[ 8]);
  Functionals[7] = ( NF_N_Q_Q2_2D_EdgeWeight1[0]*PointValues[ 9]
                    +NF_N_Q_Q2_2D_EdgeWeight1[1]*PointValues[10]
                    +NF_N_Q_Q2_2D_EdgeWeight1[2]*PointValues[11]);

  Functionals[8] = ( NF_N_Q_Q2_2D_CellWeight0[0]*PointValues[12]
                    +NF_N_Q_Q2_2D_CellWeight0[1]*PointValues[13]
                    +NF_N_Q_Q2_2D_CellWeight0[2]*PointValues[14]
                    +NF_N_Q_Q2_2D_CellWeight0[3]*PointValues[15]
                    +NF_N_Q_Q2_2D_CellWeight0[4]*PointValues[16]
                    +NF_N_Q_Q2_2D_CellWeight0[5]*PointValues[17]
                    +NF_N_Q_Q2_2D_CellWeight0[6]*PointValues[18]
                    +NF_N_Q_Q2_2D_CellWeight0[7]*PointValues[19]
                    +NF_N_Q_Q2_2D_CellWeight0[8]*PointValues[20]);
                    
  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[4] = -Functionals[4];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[5] = -Functionals[5];
    if(Cell->GetVertex(2) > Cell->GetVertex(3))
      Functionals[6] = -Functionals[6];
    if(Cell->GetVertex(3) > Cell->GetVertex(0))
      Functionals[7] = -Functionals[7];
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
        Functionals[ 4] = -Functionals[ 4];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 5] = -Functionals[ 5];
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 6] = -Functionals[ 6];
    } // endif neigh

    neigh = Cell->GetJoint(3)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 7] = -Functionals[ 7];
    } // endif neigh
  } // endif Cell
}

void NF_N_Q_Q2_2D_EvalEdge(const TCollection *Coll, const TBaseCell *Cell, int Joint,
                           const double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_Q_Q2_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_Q_Q2_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_Q_Q2_2D_EdgeWeight0[2]*PointValues[2]);
  Functionals[1] = ( NF_N_Q_Q2_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_Q_Q2_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_Q_Q2_2D_EdgeWeight1[2]*PointValues[2]);

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%4))
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
