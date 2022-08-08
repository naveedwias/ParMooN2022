// Second order Raviart-Thomas vector element, nonconforming, 2D

// points for 1D Gauss quadrature with four points (two are symmetric)
// { -0.86114, -0.33998, 0.33998, 0.86114 }
static double NF_N_T_RT2_2D_T[4] = 
 {-std::sqrt(3./7. + (2./7.)*std::sqrt(6./5.)), -std::sqrt(3./7. - (2./7.)*std::sqrt(6./5.)),
   std::sqrt(3./7. - (2./7.)*std::sqrt(6./5.)),  std::sqrt(3./7. + (2./7.)*std::sqrt(6./5.)) };

// the points on the line [0,1]
static double NF_N_T_RT2_2D_q[4] = 
{ 0.5*NF_N_T_RT2_2D_T[0] + 0.5, 0.5*NF_N_T_RT2_2D_T[1] + 0.5,
  0.5*NF_N_T_RT2_2D_T[2] + 0.5, 0.5*NF_N_T_RT2_2D_T[3] + 0.5 };
// weights for the 1D Gauss quadrature with four points
// { 0.34785, 0.65215, 0.65215, 0.34785 }
static double NF_N_T_RT2_2D_w[4] =
{ (18. - std::sqrt(30.)) / 72., (18. + std::sqrt(30.)) / 72.,
  (18. + std::sqrt(30.)) / 72., (18. - std::sqrt(30.)) / 72. }; 

// P_2(x) with P_2 being the second Legendre polynomial and x the quad points 
// from above
static double NF_N_T_RT2_2D_p2[4] = 
  { 0.5*(3*NF_N_T_RT2_2D_T[0]*NF_N_T_RT2_2D_T[0] - 1.),
    0.5*(3*NF_N_T_RT2_2D_T[1]*NF_N_T_RT2_2D_T[1] - 1.),
    0.5*(3*NF_N_T_RT2_2D_T[2]*NF_N_T_RT2_2D_T[2] - 1.),
    0.5*(3*NF_N_T_RT2_2D_T[3]*NF_N_T_RT2_2D_T[3] - 1.) };

// four point Gauss quadrature for edge dofs, 19-point formula for inner dofs
static double NF_N_T_RT2_2D_Xi[] = 
{ NF_N_T_RT2_2D_q[0], NF_N_T_RT2_2D_q[1] , NF_N_T_RT2_2D_q[2], 
  NF_N_T_RT2_2D_q[3],
  NF_N_T_RT2_2D_q[3], NF_N_T_RT2_2D_q[2], NF_N_T_RT2_2D_q[1], 
  NF_N_T_RT2_2D_q[0],
  0, 0, 0, 0,
  0.333333333333333333333333333333333,
  0.489682519198737627783706924836192,
  0.489682519198737627783706924836192,
  0.020634961602524744432586150327616,
  0.437089591492936637269930364435354,
  0.437089591492936637269930364435354,
  0.125820817014126725460139271129292,
  0.188203535619032730240961280467335,
  0.188203535619032730240961280467335,
  0.623592928761934539518077439065330,
  0.0447295133944527098651065899662763,
  0.0447295133944527098651065899662763,
  0.910540973211094580269786820067447,
  0.741198598784498020690079873523423,
  0.0368384120547362836348175987833851,
  0.741198598784498020690079873523423,
  0.221962989160765695675102527693192,
  0.0368384120547362836348175987833851,
  0.221962989160765695675102527693192 };
static double NF_N_T_RT2_2D_Eta[] = 
{ 0, 0, 0, 0,
  NF_N_T_RT2_2D_q[0], NF_N_T_RT2_2D_q[1] , NF_N_T_RT2_2D_q[2], 
  NF_N_T_RT2_2D_q[3],
  NF_N_T_RT2_2D_q[3], NF_N_T_RT2_2D_q[2], NF_N_T_RT2_2D_q[1], 
  NF_N_T_RT2_2D_q[0],
  0.333333333333333333333333333333333,
  0.489682519198737627783706924836192,
  0.020634961602524744432586150327616,
  0.489682519198737627783706924836192,
  0.437089591492936637269930364435354,
  0.125820817014126725460139271129292,
  0.437089591492936637269930364435354,
  0.188203535619032730240961280467335,
  0.623592928761934539518077439065330,
  0.188203535619032730240961280467335,
  0.0447295133944527098651065899662763,
  0.910540973211094580269786820067447,
  0.0447295133944527098651065899662763,
  0.0368384120547362836348175987833851,
  0.741198598784498020690079873523423,
  0.221962989160765695675102527693192,
  0.741198598784498020690079873523423,
  0.221962989160765695675102527693192,
  0.0368384120547362836348175987833851 };

double NF_N_T_RT2_2D_iw[19] = 
{ 0.0485678981413994169096209912536443,
  0.0156673501135695352684274156436046,
  0.0156673501135695352684274156436046,
  0.0156673501135695352684274156436046,
  0.0389137705023871396583696781497019,
  0.0389137705023871396583696781497019,
  0.0389137705023871396583696781497019,
  0.0398238694636051265164458871320226,
  0.0398238694636051265164458871320226,
  0.0398238694636051265164458871320226,
  0.0127888378293490156308393992794999,
  0.0127888378293490156308393992794999,
  0.0127888378293490156308393992794999,
  0.0216417696886446886446886446886446,
  0.0216417696886446886446886446886446,
  0.0216417696886446886446886446886446,
  0.0216417696886446886446886446886446,
  0.0216417696886446886446886446886446,
  0.0216417696886446886446886446886446 };

void NF_N_T_RT2_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  // set all Functionals to zero at first
  for(unsigned int i = 0; i < 15; ++i)
   Functionals[i] = 0.;
  
  // short names
  const double * ew = NF_N_T_RT2_2D_w;
  const double * iw = NF_N_T_RT2_2D_iw;
  const double * p = NF_N_T_RT2_2D_T;
  const double * p2= NF_N_T_RT2_2D_p2;
  
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    // edge dofs
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[0] -= PointValues[31+i] * ew[i];
      Functionals[1] -= PointValues[31+i] * ew[i] * p[i];
      Functionals[2] -= PointValues[31+i] * ew[i] * p2[i];
      
      Functionals[3] += (PointValues[4+i] + PointValues[35+i]) * ew[i];
      Functionals[4] += (PointValues[4+i] + PointValues[35+i]) * ew[i] * p[i];
      Functionals[5] += (PointValues[4+i] + PointValues[35+i]) * ew[i] * p2[i];
      
      Functionals[6] -= PointValues[8+i] * ew[i];
      Functionals[7] -= PointValues[8+i] * ew[i] * p[i];
      Functionals[8] -= PointValues[8+i] * ew[i] * p2[i];
    }
    //inner dofs
    for(unsigned int i = 0; i < 19; ++i)
    {
      Functionals[9]  += PointValues[12+i] * iw[i];
      Functionals[10] += PointValues[43+i] * iw[i];
      
      Functionals[11] += PointValues[12+i] * iw[i] * NF_N_T_RT2_2D_Xi[12+i];
      Functionals[12] += PointValues[43+i] * iw[i] * NF_N_T_RT2_2D_Xi[12+i];
      
      Functionals[13] += PointValues[12+i] * iw[i] * NF_N_T_RT2_2D_Eta[12+i];
      Functionals[14] += PointValues[43+i] * iw[i] * NF_N_T_RT2_2D_Eta[12+i];
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
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[31+i]*ny)*ew[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[31+i]*ny)*ew[i]*p[i];
      Functionals[2] += (PointValues[i]*nx + PointValues[31+i]*ny)*ew[i]*p2[i];
    }
    Functionals[0] *= Cell->GetNormalOrientation(0);
    //Functionals[1] *= Cell->GetNormalOrientation(0);
    Functionals[2] *= Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[3] += (PointValues[4+i]*nx +PointValues[35+i]*ny)*ew[i];
      Functionals[4] += (PointValues[4+i]*nx +PointValues[35+i]*ny)*ew[i]*p[i];
      Functionals[5] += (PointValues[4+i]*nx +PointValues[35+i]*ny)*ew[i]*p2[i];
    }
    Functionals[3] *= Cell->GetNormalOrientation(1);
    //Functionals[4] *= Cell->GetNormalOrientation(1);
    Functionals[5] *= Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    for(unsigned int i = 0; i < 4; ++i)
    {
      Functionals[6] += (PointValues[8+i]*nx +PointValues[39+i]*ny)*ew[i];
      Functionals[7] += (PointValues[8+i]*nx +PointValues[39+i]*ny)*ew[i]*p[i];
      Functionals[8] += (PointValues[8+i]*nx +PointValues[39+i]*ny)*ew[i]*p2[i];
    }
    Functionals[6] *= Cell->GetNormalOrientation(2);
    //Functionals[7] *= Cell->GetNormalOrientation(2);
    Functionals[8] *= Cell->GetNormalOrientation(2);
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 2*Cell->GetMeasure();
    
    TTriaAffin referenceTransform;
    referenceTransform.SetCell(Cell);
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    // the dofs:
    // int_cell v . (1 0)^T 
    // int_cell v . (0 1)^T 
    // int_cell v . (x 0)^T
    // int_cell v . (0 x)^T
    // int_cell v . (y 0)^T
    // int_cell v . (0 y)^T
    for(unsigned int i = 0; i < 19; ++i) // loop over inner quadrature points
    {
      double x = NF_N_T_RT2_2D_Xi[12+i], y = NF_N_T_RT2_2D_Eta[12+i];
      double uxi[6] = {  1,0,  x,0,  y,0 };
      double ueta[6] = { 0,1,  0,x,  0,y };
      for(unsigned int d = 0; d < 6; ++d) // inner dofs
      {
        double uref = 0., uxiref = uxi[d], uetaref = ueta[d];
        double uorig, uxorig, uyorig;
        referenceTransform.GetOrigValues(x, y, 1, &uref, &uxiref, &uetaref, 
                                         &uorig, &uxorig, &uyorig);
        Functionals[9+d] += ( PointValues[12+i]*uxorig
                              +PointValues[43+i]*uyorig ) * iw[i];
      }
    }
    for(unsigned int d = 0; d < 6; ++d) 
    {
      Functionals[9+d] *= measure;
    }
  }
}

void NF_N_T_RT2_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
{
  // this is needed for setting boundary conditions
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1, y1, z); // 3=number of edges
  // length of joint
  double l = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); 
  Functionals[0] = ( NF_N_T_RT2_2D_w[0]*PointValues[0]
                    +NF_N_T_RT2_2D_w[1]*PointValues[1]
                    +NF_N_T_RT2_2D_w[2]*PointValues[2]
                    +NF_N_T_RT2_2D_w[3]*PointValues[3] )*l;
  Functionals[1] = ( NF_N_T_RT2_2D_w[0]*NF_N_T_RT2_2D_T[0]*PointValues[0]
                    +NF_N_T_RT2_2D_w[1]*NF_N_T_RT2_2D_T[1]*PointValues[1]
                    +NF_N_T_RT2_2D_w[2]*NF_N_T_RT2_2D_T[2]*PointValues[2]
                    +NF_N_T_RT2_2D_w[3]*NF_N_T_RT2_2D_T[3]*PointValues[3] )*l;
  Functionals[2] = ( NF_N_T_RT2_2D_w[0]*NF_N_T_RT2_2D_p2[0]*PointValues[0]
                    +NF_N_T_RT2_2D_w[1]*NF_N_T_RT2_2D_p2[1]*PointValues[1]
                    +NF_N_T_RT2_2D_w[2]*NF_N_T_RT2_2D_p2[2]*PointValues[2]
                    +NF_N_T_RT2_2D_w[3]*NF_N_T_RT2_2D_p2[3]*PointValues[3])*l;
}
