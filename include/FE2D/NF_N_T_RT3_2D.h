// Third order Raviart-Thomas vector element, nonconforming, 2D

// points for 1D Gauss quadrature with five points
static double NF_N_T_RT3_2D_T[5] = 
{ -std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3., -std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3., 0.,
   std::sqrt(5. - 2.*std::sqrt(10./7.)) / 3.,  std::sqrt(5. + 2.*std::sqrt(10./7.)) / 3. };

// the points on the line [0,1]
static double NF_N_T_RT3_2D_q[5] = 
{ 0.5*NF_N_T_RT3_2D_T[0] + 0.5, 0.5*NF_N_T_RT3_2D_T[1] + 0.5,
  0.5*NF_N_T_RT3_2D_T[2] + 0.5, 0.5*NF_N_T_RT3_2D_T[3] + 0.5,
  0.5*NF_N_T_RT3_2D_T[4] + 0.5
};
// weights for the 1D Gauss quadrature with four points
// { 0.34785, 0.65215, 0.65215, 0.34785 }
static double NF_N_T_RT3_2D_w[5] =
{ (322. - 13.*std::sqrt(70.))/1800., (322. + 13.*std::sqrt(70.))/1800., 128./450.,
  (322. + 13.*std::sqrt(70.))/1800., (322. - 13.*std::sqrt(70.))/1800. }; 

// P_2(x) with P_2 being the second Legendre polynomial and x the quad points 
// from above
static double NF_N_T_RT3_2D_p2[5] = 
 { 0.5*(3*NF_N_T_RT3_2D_T[0]*NF_N_T_RT3_2D_T[0] - 1.), 
   0.5*(3*NF_N_T_RT3_2D_T[1]*NF_N_T_RT3_2D_T[1] - 1.),
   0.5*(3*NF_N_T_RT3_2D_T[2]*NF_N_T_RT3_2D_T[2] - 1.),
   0.5*(3*NF_N_T_RT3_2D_T[3]*NF_N_T_RT3_2D_T[3] - 1.),
   0.5*(3*NF_N_T_RT3_2D_T[4]*NF_N_T_RT3_2D_T[4] - 1.) };
// P_3(x) with P_3 being the third Legendre polynomial and x the Gauss 
// quadrature points from above
static double NF_N_T_RT3_2D_p3[5] = 
{ 0.5*( 5.*NF_N_T_RT3_2D_T[0]*NF_N_T_RT3_2D_T[0]*NF_N_T_RT3_2D_T[0] 
       -3.*NF_N_T_RT3_2D_T[0]),
  0.5*( 5.*NF_N_T_RT3_2D_T[1]*NF_N_T_RT3_2D_T[1]*NF_N_T_RT3_2D_T[1] 
       -3.*NF_N_T_RT3_2D_T[1]),
  0.5*( 5.*NF_N_T_RT3_2D_T[2]*NF_N_T_RT3_2D_T[2]*NF_N_T_RT3_2D_T[2] 
       -3.*NF_N_T_RT3_2D_T[2]),
  0.5*( 5.*NF_N_T_RT3_2D_T[3]*NF_N_T_RT3_2D_T[3]*NF_N_T_RT3_2D_T[3] 
       -3.*NF_N_T_RT3_2D_T[3]),
  0.5*( 5.*NF_N_T_RT3_2D_T[4]*NF_N_T_RT3_2D_T[4]*NF_N_T_RT3_2D_T[4] 
       -3.*NF_N_T_RT3_2D_T[4]),
};

// four point Gauss quadrature for edge dofs, 19-point formula for inner dofs
static double NF_N_T_RT3_2D_Xi[] = 
{ NF_N_T_RT3_2D_q[0], NF_N_T_RT3_2D_q[1], NF_N_T_RT3_2D_q[2], 
  NF_N_T_RT3_2D_q[3], NF_N_T_RT3_2D_q[4],
  NF_N_T_RT3_2D_q[4], NF_N_T_RT3_2D_q[3], NF_N_T_RT3_2D_q[2], 
  NF_N_T_RT3_2D_q[1], NF_N_T_RT3_2D_q[0],
  0, 0, 0, 0, 0,
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
static double NF_N_T_RT3_2D_Eta[] = 
{ 0, 0, 0, 0, 0,
  NF_N_T_RT3_2D_q[0], NF_N_T_RT3_2D_q[1] , NF_N_T_RT3_2D_q[2], 
  NF_N_T_RT3_2D_q[3], NF_N_T_RT3_2D_q[4],
  NF_N_T_RT3_2D_q[4], NF_N_T_RT3_2D_q[3], NF_N_T_RT3_2D_q[2], 
  NF_N_T_RT3_2D_q[1], NF_N_T_RT3_2D_q[0],
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

double NF_N_T_RT3_2D_iw[19] = 
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

void NF_N_T_RT3_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  for(unsigned int i = 0; i < 24; ++i)
    Functionals[i] = 0.;
  
  // short names
  const double * q = NF_N_T_RT3_2D_T;
  const double * ew = NF_N_T_RT3_2D_w;
  const double * p2 = NF_N_T_RT3_2D_p2;
  const double * p3 = NF_N_T_RT3_2D_p3;
  const double * iw = NF_N_T_RT3_2D_iw;
  const double * xi = NF_N_T_RT3_2D_Xi+15; // +15 to skip the points on edges
  const double * et = NF_N_T_RT3_2D_Eta+15; // +15 to skip the points on edges
  
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    // edge dofs
    for(unsigned int i = 0; i < 5; ++i)
    {
      // first edge
      Functionals[0] -= PointValues[34+i] * ew[i];
      Functionals[1] -= PointValues[34+i] * ew[i] * q[i];
      Functionals[2] -= PointValues[34+i] * ew[i] * p2[i];
      Functionals[3] -= PointValues[34+i] * ew[i] * p3[i];
      // second edge
      Functionals[4] += (PointValues[5+i] + PointValues[39+i]) * ew[i];
      Functionals[5] += (PointValues[5+i] + PointValues[39+i]) * ew[i] * q[i];
      Functionals[6] += (PointValues[5+i] + PointValues[39+i]) * ew[i] * p2[i];
      Functionals[7] += (PointValues[5+i] + PointValues[39+i]) * ew[i] * p3[i];
      // third edge
      Functionals[8] -= PointValues[10+i] * ew[i];
      Functionals[9] -= PointValues[10+i] * ew[i] * q[i];
      Functionals[10]-= PointValues[10+i] * ew[i] * p2[i];
      Functionals[11]-= PointValues[10+i] * ew[i] * p3[i];
    }
    // inner dofs
    for(unsigned int i = 0; i < 19; ++i)
    {
      Functionals[12] += PointValues[15+i] * iw[i];
      Functionals[13] += PointValues[49+i] * iw[i];
      
      Functionals[14] += PointValues[15+i] * iw[i] * xi[i];
      Functionals[15] += PointValues[49+i] * iw[i] * xi[i];
      
      Functionals[16] += PointValues[15+i] * iw[i] * et[i];
      Functionals[17] += PointValues[49+i] * iw[i] * et[i];
      
      Functionals[18] += PointValues[15+i] * iw[i] * xi[i] * xi[i];
      Functionals[19] += PointValues[49+i] * iw[i] * xi[i] * xi[i];
      
      Functionals[20] += PointValues[15+i] * iw[i] * et[i] * et[i];
      Functionals[21] += PointValues[49+i] * iw[i] * et[i] * et[i];
      
      Functionals[22] += PointValues[15+i] * iw[i] * xi[i] * et[i];
      Functionals[23] += PointValues[49+i] * iw[i] * xi[i] * et[i];
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
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[0] += (PointValues[i]*nx + PointValues[34+i]*ny)*ew[i];
      Functionals[1] += (PointValues[i]*nx + PointValues[34+i]*ny)*ew[i] *q[i];
      Functionals[2] += (PointValues[i]*nx + PointValues[34+i]*ny)*ew[i] *p2[i];
      Functionals[3] += (PointValues[i]*nx + PointValues[34+i]*ny)*ew[i] *p3[i];
    }
    Functionals[0] *= Cell->GetNormalOrientation(0);
    //Functionals[1] *= Cell->GetNormalOrientation(0);
    Functionals[2] *= Cell->GetNormalOrientation(0);
    //Functionals[3] *= Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[4]+= (PointValues[5+i]*nx + PointValues[39+i]*ny)*ew[i];
      Functionals[5]+= (PointValues[5+i]*nx + PointValues[39+i]*ny)*ew[i]*q[i];
      Functionals[6]+= (PointValues[5+i]*nx + PointValues[39+i]*ny)*ew[i]*p2[i];
      Functionals[7]+= (PointValues[5+i]*nx + PointValues[39+i]*ny)*ew[i]*p3[i];
    }
    Functionals[4] *= Cell->GetNormalOrientation(1);
    //Functionals[5] *= Cell->GetNormalOrientation(1);
    Functionals[6] *= Cell->GetNormalOrientation(1);
    //Functionals[7] *= Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    for(unsigned int i = 0; i < 5; ++i)
    {
      Functionals[8] += (PointValues[10+i]*nx+PointValues[44+i]*ny)*ew[i];
      Functionals[9] += (PointValues[10+i]*nx+PointValues[44+i]*ny)*ew[i]*q[i];
      Functionals[10]+= (PointValues[10+i]*nx+PointValues[44+i]*ny)*ew[i]*p2[i];
      Functionals[11]+= (PointValues[10+i]*nx+PointValues[44+i]*ny)*ew[i]*p3[i];
    }
    Functionals[8] *= Cell->GetNormalOrientation(2);
    //Functionals[9] *= Cell->GetNormalOrientation(2);
    Functionals[10]*= Cell->GetNormalOrientation(2);
    //Functionals[11]*= Cell->GetNormalOrientation(2);
    
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
    // int_cell v . (xx 0)^T
    // int_cell v . (0 xx)^T
    // int_cell v . (yy 0)^T
    // int_cell v . (0 yy)^T
    // int_cell v . (xy 0)^T
    // int_cell v . (0 xy)^T
    for(unsigned int i = 0; i < 19; ++i) // loop over inner quadrature points
    {
      double x = xi[i], y = et[i];
      double uxi[12] = {  1,0,  x,0,  y,0,  x*x,0,  y*y,0,  x*y,0 };
      double ueta[12] = { 0,1,  0,x,  0,y,  0,x*x,  0,y*y,  0,x*y };
      for(unsigned int d = 0; d < 12; ++d) // inner dofs
      {
        double uref = 0., uxiref = uxi[d], uetaref = ueta[d];
        double uorig, uxorig, uyorig;
        referenceTransform.GetOrigValues(x, y, 1, &uref, &uxiref, &uetaref, 
                                         &uorig, &uxorig, &uyorig);
        Functionals[12+d] += ( PointValues[15+i]*uxorig
                              +PointValues[49+i]*uyorig ) * iw[i];
      }
    }
    for(unsigned int d = 0; d < 12; ++d) 
    {
      Functionals[12+d] *= measure;
    }
  }
}

void NF_N_T_RT3_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
{
  // this is needed for setting boundary conditions
  Functionals[0] = 0.;
  Functionals[1] = 0.;
  Functionals[2] = 0.;
  Functionals[3] = 0.;
  for(unsigned int i = 0; i < 5; ++i)
  {
    double pv = PointValues[i];
    Functionals[0] += pv * NF_N_T_RT3_2D_w[i];
    Functionals[1] += pv * NF_N_T_RT3_2D_w[i] * NF_N_T_RT3_2D_T[i];
    Functionals[2] += pv * NF_N_T_RT3_2D_w[i] * NF_N_T_RT3_2D_p2[i];
    Functionals[3] += pv * NF_N_T_RT3_2D_w[i] * NF_N_T_RT3_2D_p3[i];
  }
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1, y1, z);// 3=number of edges
  // length of joint
  double l = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] *= l;
  Functionals[1] *= l;
  Functionals[2] *= l;
  Functionals[3] *= l;
}
