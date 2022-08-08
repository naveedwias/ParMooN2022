static double NF_N_T_BDM1_2D_T[3] = {-std::sqrt(3./5.), 0, std::sqrt(3./5.)};

// three point Gauss quadrature for edge dofs, no inner dofs
static const double NF_N_T_BDM1_2D_eq[3] = 
 { 0.5*NF_N_T_BDM1_2D_T[0]+0.5,  0.5*NF_N_T_BDM1_2D_T[1]+0.5,
   0.5*NF_N_T_BDM1_2D_T[2]+0.5 };
static const double NF_N_T_BDM1_2D_ew[3] = { 5./18., 8./18., 5./18. };

static double NF_N_T_BDM1_2D_Xi[9] = 
 { NF_N_T_BDM1_2D_eq[0], NF_N_T_BDM1_2D_eq[1], NF_N_T_BDM1_2D_eq[2],
   NF_N_T_BDM1_2D_eq[2], NF_N_T_BDM1_2D_eq[1], NF_N_T_BDM1_2D_eq[0],
   0, 0, 0 };
static double NF_N_T_BDM1_2D_Eta[9] = 
 { 0, 0, 0,
   NF_N_T_BDM1_2D_eq[0], NF_N_T_BDM1_2D_eq[1], NF_N_T_BDM1_2D_eq[2],
   NF_N_T_BDM1_2D_eq[2], NF_N_T_BDM1_2D_eq[1], NF_N_T_BDM1_2D_eq[0] };

void NF_N_T_BDM1_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                            const double *PointValues, double *Functionals)
{
  // short names
  const double * p = NF_N_T_BDM1_2D_T;
  const double * ew = NF_N_T_BDM1_2D_ew;
  
  // set all Functionals to zero first
  for(unsigned int i = 0; i < 6; ++i)
    Functionals[i] = 0;
  
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    // outer dofs
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[0] -= PointValues[9+i] * ew[i];
      Functionals[1] -= PointValues[9+i] * ew[i] * p[i];
      
      Functionals[2] += (PointValues[3+i] + PointValues[12+i]) * ew[i];
      Functionals[3] += (PointValues[3+i] + PointValues[12+i]) * ew[i] * p[i];
      
      Functionals[4] -= PointValues[6+i] * ew[i];
      Functionals[5] -= PointValues[6+i] * ew[i] * p[i];
    }
  }
  else // on a real cell
  {
    double x0, x1, x2, y0, y1, y2, z; // z remains zero in 2D
    Cell->GetVertex(0)->GetCoords(x0, y0, z);
    Cell->GetVertex(1)->GetCoords(x1, y1, z);
    Cell->GetVertex(2)->GetCoords(x2, y2, z);
    // length of edge, and outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[9+i]*ny)*ew[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[9+i]*ny)*ew[i]*p[i];
    }
    Functionals[0] *= Cell->GetNormalOrientation(0);
    //Functionals[1] *= Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[2] += (PointValues[3+i]*nx + PointValues[12+i]*ny)*ew[i];
      Functionals[3] += (PointValues[3+i]*nx + PointValues[12+i]*ny)*ew[i]*p[i];
    }
    Functionals[2] *= Cell->GetNormalOrientation(1);
    //Functionals[3] *= Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    for(unsigned int i = 0; i < 3; ++i)
    {
      Functionals[4] += (PointValues[6+i]*nx + PointValues[15+i]*ny)*ew[i];
      Functionals[5] += (PointValues[6+i]*nx + PointValues[15+i]*ny)*ew[i]*p[i];
    }
    Functionals[4] *= Cell->GetNormalOrientation(2);
    //Functionals[5] *= Cell->GetNormalOrientation(2);
  }
}

void NF_N_T_BDM1_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int Joint,
                             const double *PointValues, double *Functionals)
{
  // this is needed for setting boundary conditions
  Functionals[0] = 0.;
  Functionals[1] = 0.;
  for(unsigned int  i = 0; i < 3; ++i)
  {
    Functionals[0] += PointValues[i] * NF_N_T_BDM1_2D_ew[i];
    Functionals[1] += PointValues[i] * NF_N_T_BDM1_2D_ew[i]*NF_N_T_BDM1_2D_T[i];
  }
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1, y1, z); // 3=number of edges
  // length of joint
  const double l = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] *= l;
  Functionals[1] *= l;
}
