// Zeroth order Raviart-Thomas vector element, nonconforming, 2D

// the nodal functionals on the edges E_i are
// N_i(v) = int_{E_i} v.n 
// we use a 2-point Gauss quadrature on each edge (which is exact for 
// polynomials up to order 3)
static double NF_N_T_RT0_2D_Xi[6] =
{ (-std::sqrt(1./3.)+1.)/2., (std::sqrt(1./3.)+1.)/2.,
  (std::sqrt(1./3.)+1.)/2., (-std::sqrt(1./3.)+1.)/2.,
  0, 0};
static double NF_N_T_RT0_2D_Eta[6] = 
{ 0, 0 , 
  (-std::sqrt(1./3.)+1.)/2., (std::sqrt(1./3.)+1.)/2., 
  (std::sqrt(1./3.)+1.)/2., (-std::sqrt(1./3.)+1.)/2.};
static double NF_N_T_RT0_2D_T[] = {-std::sqrt(1./3.), std::sqrt(1./3.)};

void NF_N_T_RT0_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    Functionals[0] = -0.5*(PointValues[6] + PointValues[7]);
    Functionals[1] = 0.5*(PointValues[2] + PointValues[3]) 
                    +0.5*(PointValues[8] + PointValues[9]);
    Functionals[2] = -0.5*(PointValues[4] + PointValues[5]);
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
    Functionals[0] = (PointValues[0] + PointValues[1])*nx
                    +(PointValues[6] + PointValues[7])*ny;
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[1] = (PointValues[2] + PointValues[3])*nx
                    +(PointValues[8] + PointValues[9])*ny;
    Functionals[1] *= 0.5*Cell->GetNormalOrientation(1);
     
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    Functionals[2] = (PointValues[4] + PointValues[5])*nx
                    +(PointValues[10]+ PointValues[11])*ny;
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(2);
  }
}

void NF_N_T_RT0_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
{
  // this is needed for setting boundary conditions
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1, y1, z);// 3=number of edges
  // length of joint, 0.5 due to 1D-reference cell having measure 2
  double l = 0.5*std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] = (PointValues[0] + PointValues[1])*l;
}
