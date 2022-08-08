// Second order Raviart-Thomas vector element on quads, nonconforming, 2D

// points for 1D Gauss quadrature with four points (two are symmetric)
static double NF_N_Q_RT2_2D_a = std::sqrt(3./7. + (2./7.)*std::sqrt(6./5.)); // 0.86114
static double NF_N_Q_RT2_2D_b = std::sqrt(3./7. - (2./7.)*std::sqrt(6./5.)); // 0.33998
// weights for the 1D Gauss quadrature with four points
static double NF_N_Q_RT2_2D_wa = (18. - std::sqrt(30.)) / 36.; // 0.34785
static double NF_N_Q_RT2_2D_wb = (18. + std::sqrt(30.)) / 36.; // 0.65215
// P_2(x) with P_2 being the second Legendre polynomial and x the two quad 
// points from above
static double NF_N_Q_RT2_2D_p2a = 0.5*(3*NF_N_Q_RT2_2D_a*NF_N_Q_RT2_2D_a - 1.);
static double NF_N_Q_RT2_2D_p2b = 0.5*(3*NF_N_Q_RT2_2D_b*NF_N_Q_RT2_2D_b - 1.);

static double NF_N_Q_RT2_2D_Xi[] = 
{ -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b , NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a,
   1, 1, 1, 1,
   NF_N_Q_RT2_2D_a, NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_a,
  -1, -1, -1, -1,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a
};
static double NF_N_Q_RT2_2D_Eta[] = 
{ -1, -1, -1, -1,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a,
   1, 1, 1, 1,
   NF_N_Q_RT2_2D_a, NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_a,
  -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_a,
  -NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_b, -NF_N_Q_RT2_2D_b,
   NF_N_Q_RT2_2D_b,  NF_N_Q_RT2_2D_b,  NF_N_Q_RT2_2D_b,  NF_N_Q_RT2_2D_b,
   NF_N_Q_RT2_2D_a,  NF_N_Q_RT2_2D_a,  NF_N_Q_RT2_2D_a,  NF_N_Q_RT2_2D_a
};

static double NF_N_Q_RT2_2D_T[] = 
 { -NF_N_Q_RT2_2D_a, -NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_b, NF_N_Q_RT2_2D_a };



void NF_N_Q_RT2_2D_EvalAll(const TCollection *, const TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  // short names
  const double a = NF_N_Q_RT2_2D_a;
  const double b = NF_N_Q_RT2_2D_b;
  const double wa = NF_N_Q_RT2_2D_wa;
  const double wb = NF_N_Q_RT2_2D_wb;
  const double p2a = NF_N_Q_RT2_2D_p2a;
  const double p2b = NF_N_Q_RT2_2D_p2b;
  
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -( wa*PointValues[32] + wb*PointValues[33]
                       +wb*PointValues[34] + wa*PointValues[35] );
    Functionals[1] = -(-wa*a*PointValues[32] - wb*b*PointValues[33]
                       +wb*b*PointValues[34] + wa*a*PointValues[35] );
    Functionals[2] = -( wa*p2a*PointValues[32] + wb*p2b*PointValues[33]
                       +wb*p2b*PointValues[34] + wa*p2a*PointValues[35] );
    Functionals[3] =  wa*PointValues[4] + wb*PointValues[5]
                     +wb*PointValues[6] + wa*PointValues[7];
    Functionals[4] = -wa*a*PointValues[4] - wb*b*PointValues[5]
                     +wb*b*PointValues[6] + wa*a*PointValues[7];
    Functionals[5] = +wa*p2a*PointValues[4] + wb*p2b*PointValues[5]
                     +wb*p2b*PointValues[6] + wa*p2a*PointValues[7];
    Functionals[6] =  wa*PointValues[40] + wb*PointValues[41]
                     +wb*PointValues[42] + wa*PointValues[43];
    Functionals[7] = -wa*a*PointValues[40] - wb*b*PointValues[41]
                     +wb*b*PointValues[42] + wa*a*PointValues[43];
    Functionals[8] =  wa*p2a*PointValues[40] + wb*p2b*PointValues[41]
                     +wb*p2b*PointValues[42] + wa*p2a*PointValues[43];
    Functionals[9] = -( wa*PointValues[12] + wb*PointValues[13]
                       +wb*PointValues[14] + wa*PointValues[15] );
    Functionals[10]= -(-wa*a*PointValues[12] - wb*b*PointValues[13]
                       +wb*b*PointValues[14] + wa*a*PointValues[15] );
    Functionals[11]= -( wa*p2a*PointValues[12] + wb*p2b*PointValues[13]
                       +wb*p2b*PointValues[14] + wa*p2a*PointValues[15] );
    
    // int_cell v . (1 0)^T 
    Functionals[12]=  wa*wa*PointValues[16] + wa*wb*PointValues[17]
                     +wa*wb*PointValues[18] + wa*wa*PointValues[19]
                     +wa*wb*PointValues[20] + wb*wb*PointValues[21]
                     +wb*wb*PointValues[22] + wa*wb*PointValues[23]
                     +wa*wb*PointValues[24] + wb*wb*PointValues[25]
                     +wb*wb*PointValues[26] + wa*wb*PointValues[27]
                     +wa*wa*PointValues[28] + wa*wb*PointValues[29]
                     +wa*wb*PointValues[30] + wa*wa*PointValues[31];
    // int_cell v . (0 1)^T 
    Functionals[13]=  wa*wa*PointValues[48] + wa*wb*PointValues[49]
                     +wa*wb*PointValues[50] + wa*wa*PointValues[51]
                     +wa*wb*PointValues[52] + wb*wb*PointValues[53]
                     +wb*wb*PointValues[54] + wa*wb*PointValues[55]
                     +wa*wb*PointValues[56] + wb*wb*PointValues[57]
                     +wb*wb*PointValues[58] + wa*wb*PointValues[59]
                     +wa*wa*PointValues[60] + wa*wb*PointValues[61]
                     +wa*wb*PointValues[62] + wa*wa*PointValues[63];
    // int_cell v . (x 0)^T 
    Functionals[14]= -wa*wa*a*PointValues[16] - wa*wb*b*PointValues[17]
                     +wa*wb*b*PointValues[18] + wa*wa*a*PointValues[19]
                     -wa*wb*a*PointValues[20] - wb*wb*b*PointValues[21]
                     +wb*wb*b*PointValues[22] + wa*wb*a*PointValues[23]
                     -wa*wb*a*PointValues[24] - wb*wb*b*PointValues[25]
                     +wb*wb*b*PointValues[26] + wa*wb*a*PointValues[27]
                     -wa*wa*a*PointValues[28] - wa*wb*b*PointValues[29]
                     +wa*wb*b*PointValues[30] + wa*wa*a*PointValues[31];
    // int_cell v . (0 x)^T 
    Functionals[15]= -wa*wa*a*PointValues[48] - wa*wb*b*PointValues[49]
                     +wa*wb*b*PointValues[50] + wa*wa*a*PointValues[51]
                     -wa*wb*a*PointValues[52] - wb*wb*b*PointValues[53]
                     +wb*wb*b*PointValues[54] + wa*wb*a*PointValues[55]
                     -wa*wb*a*PointValues[56] - wb*wb*b*PointValues[57]
                     +wb*wb*b*PointValues[58] + wa*wb*a*PointValues[59]
                     -wa*wa*a*PointValues[60] - wa*wb*b*PointValues[61]
                     +wa*wb*b*PointValues[62] + wa*wa*a*PointValues[63];
    // int_cell v . (y 0)^T 
    Functionals[16]= -wa*wa*a*PointValues[16] - wa*wb*a*PointValues[17]
                     -wa*wb*a*PointValues[18] - wa*wa*a*PointValues[19]
                     -wa*wb*b*PointValues[20] - wb*wb*b*PointValues[21]
                     -wb*wb*b*PointValues[22] - wa*wb*b*PointValues[23]
                     +wa*wb*b*PointValues[24] + wb*wb*b*PointValues[25]
                     +wb*wb*b*PointValues[26] + wa*wb*b*PointValues[27]
                     +wa*wa*a*PointValues[28] + wa*wb*a*PointValues[29]
                     +wa*wb*a*PointValues[30] + wa*wa*a*PointValues[31];
    // int_cell v . (0 y)^T 
    Functionals[17]= -wa*wa*a*PointValues[48] - wa*wb*a*PointValues[49]
                     -wa*wb*a*PointValues[50] - wa*wa*a*PointValues[51]
                     -wa*wb*b*PointValues[52] - wb*wb*b*PointValues[53]
                     -wb*wb*b*PointValues[54] - wa*wb*b*PointValues[55]
                     +wa*wb*b*PointValues[56] + wb*wb*b*PointValues[57]
                     +wb*wb*b*PointValues[58] + wa*wb*b*PointValues[59]
                     +wa*wa*a*PointValues[60] + wa*wb*a*PointValues[61]
                     +wa*wb*a*PointValues[62] + wa*wa*a*PointValues[63];
    // int_cell v . (0 xx)^T 
    Functionals[18]=  wa*wa*a*a*PointValues[48] + wa*wb*b*b*PointValues[49]
                     +wa*wb*b*b*PointValues[50] + wa*wa*a*a*PointValues[51]
                     +wa*wb*a*a*PointValues[52] + wb*wb*b*b*PointValues[53]
                     +wb*wb*b*b*PointValues[54] + wa*wb*a*a*PointValues[55]
                     +wa*wb*a*a*PointValues[56] + wb*wb*b*b*PointValues[57]
                     +wb*wb*b*b*PointValues[58] + wa*wb*a*a*PointValues[59]
                     +wa*wa*a*a*PointValues[60] + wa*wb*b*b*PointValues[61]
                     +wa*wb*b*b*PointValues[62] + wa*wa*a*a*PointValues[63];
    // int_cell v . (yy 0)^T 
    Functionals[19]=  wa*wa*a*a*PointValues[16] + wa*wb*a*a*PointValues[17]
                     +wa*wb*a*a*PointValues[18] + wa*wa*a*a*PointValues[19]
                     +wa*wb*b*b*PointValues[20] + wb*wb*b*b*PointValues[21]
                     +wb*wb*b*b*PointValues[22] + wa*wb*b*b*PointValues[23]
                     +wa*wb*b*b*PointValues[24] + wb*wb*b*b*PointValues[25]
                     +wb*wb*b*b*PointValues[26] + wa*wb*b*b*PointValues[27]
                     +wa*wa*a*a*PointValues[28] + wa*wb*a*a*PointValues[29]
                     +wa*wb*a*a*PointValues[30] + wa*wa*a*a*PointValues[31];
    // int_cell v . (xy 0)^T 
    Functionals[20]=  wa*wa*a*a*PointValues[16] + wa*wb*a*b*PointValues[17]
                     -wa*wb*a*b*PointValues[18] - wa*wa*a*a*PointValues[19]
                     +wa*wb*a*b*PointValues[20] + wb*wb*b*b*PointValues[21]
                     -wb*wb*b*b*PointValues[22] - wa*wb*a*b*PointValues[23]
                     -wa*wb*a*b*PointValues[24] - wb*wb*b*b*PointValues[25]
                     +wb*wb*b*b*PointValues[26] + wa*wb*a*b*PointValues[27]
                     -wa*wa*a*a*PointValues[28] - wa*wb*a*b*PointValues[29]
                     +wa*wb*b*a*PointValues[30] + wa*wa*a*a*PointValues[31];
    // int_cell v . (0 xy)^T 
    Functionals[21]=  wa*wa*a*a*PointValues[48] + wa*wb*a*b*PointValues[49]
                     -wa*wb*a*b*PointValues[50] - wa*wa*a*a*PointValues[51]
                     +wa*wb*a*b*PointValues[52] + wb*wb*b*b*PointValues[53]
                     -wb*wb*b*b*PointValues[54] - wa*wb*a*b*PointValues[55]
                     -wa*wb*a*b*PointValues[56] - wb*wb*b*b*PointValues[57]
                     +wb*wb*b*b*PointValues[58] + wa*wb*a*b*PointValues[59]
                     -wa*wa*a*a*PointValues[60] - wa*wb*b*a*PointValues[61]
                     +wa*wb*b*a*PointValues[62] + wa*wa*a*a*PointValues[63];
    // int_cell v . (0 xxy)^T 
    Functionals[22]= -wa*wa*a*a*a*PointValues[48] - wa*wb*a*b*b*PointValues[49]
                     -wa*wb*a*b*b*PointValues[50] - wa*wa*a*a*a*PointValues[51]
                     -wa*wb*a*a*b*PointValues[52] - wb*wb*b*b*b*PointValues[53]
                     -wb*wb*b*b*b*PointValues[54] - wa*wb*a*a*b*PointValues[55]
                     +wa*wb*a*b*a*PointValues[56] + wb*wb*b*b*b*PointValues[57]
                     +wb*wb*b*b*b*PointValues[58] + wa*wb*a*a*b*PointValues[59]
                     +wa*wa*a*a*a*PointValues[60] + wa*wb*a*b*b*PointValues[61]
                     +wa*wb*a*b*b*PointValues[62] + wa*wa*a*a*a*PointValues[63];
    // int_cell v . (xyy 0)^T 
    Functionals[23]= -wa*wa*a*a*a*PointValues[16] - wa*wb*a*a*b*PointValues[17]
                     +wa*wb*a*a*b*PointValues[18] + wa*wa*a*a*a*PointValues[19]
                     -wa*wb*a*b*b*PointValues[20] - wb*wb*b*b*b*PointValues[21]
                     +wb*wb*b*b*b*PointValues[22] + wa*wb*a*b*b*PointValues[23]
                     -wa*wb*a*b*b*PointValues[24] - wb*wb*b*b*b*PointValues[25]
                     +wb*wb*b*b*b*PointValues[26] + wa*wb*a*b*b*PointValues[27]
                     -wa*wa*a*a*a*PointValues[28] - wa*wb*a*b*a*PointValues[29]
                     +wa*wb*b*a*a*PointValues[30] + wa*wa*a*a*a*PointValues[31];
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_RT2_2D_EvalAll not tested for non affine ",
               "reference transformations");
    }
    double x0, x1, x2, x3, y0, y1, y2, y3, z; // z remains zero in 2D
    Cell->GetVertex(0)->GetCoords(x0, y0, z);
    Cell->GetVertex(1)->GetCoords(x1, y1, z);
    Cell->GetVertex(2)->GetCoords(x2, y2, z);
    Cell->GetVertex(3)->GetCoords(x3, y3, z);
    
    // more short names 
    double * xi = NF_N_Q_RT2_2D_Xi;
    double *eta = NF_N_Q_RT2_2D_Eta;
    
    // outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = ( wa*PointValues[0]  + wb*PointValues[1]
                      +wb*PointValues[2]  + wa*PointValues[3] )*nx
                    +( wa*PointValues[32] + wb*PointValues[33]
                      +wb*PointValues[34] + wa*PointValues[35] )*ny;
    Functionals[1] = (-wa*a*PointValues[0] - wb*b*PointValues[1]
                      +wb*b*PointValues[2] + wa*a*PointValues[3] )*nx
                    +(-wa*a*PointValues[32] - wb*b*PointValues[33]
                      +wb*b*PointValues[34] + wa*a*PointValues[35] )*ny;
    Functionals[2] = ( wa*p2a*PointValues[0]  + wb*p2b*PointValues[1]
                      +wb*p2b*PointValues[2]  + wa*p2a*PointValues[3] )*nx
                    +( wa*p2a*PointValues[32] + wb*p2b*PointValues[33]
                      +wb*p2b*PointValues[34] + wa*p2a*PointValues[35] )*ny;
    Functionals[0] *= 0.5*Cell->GetNormalOrientation(0);
    Functionals[1] *= 0.5; // Cell->GetNormalOrientation(0);
    Functionals[2] *= 0.5*Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[3] = ( wa*PointValues[4]  + wb*PointValues[5]
                      +wb*PointValues[6]  + wa*PointValues[7] )*nx
                    +( wa*PointValues[36] + wb*PointValues[37]
                      +wb*PointValues[38] + wa*PointValues[39] )*ny;
    Functionals[4] = (-wa*a*PointValues[4] - wb*b*PointValues[5]
                      +wb*b*PointValues[6] + wa*a*PointValues[7] )*nx
                    +(-wa*a*PointValues[36] - wb*b*PointValues[37]
                      +wb*b*PointValues[38] + wa*a*PointValues[39] )*ny;
    Functionals[5] = ( wa*p2a*PointValues[4]  + wb*p2b*PointValues[5]
                      +wb*p2b*PointValues[6]  + wa*p2a*PointValues[7] )*nx
                    +( wa*p2a*PointValues[36] + wb*p2b*PointValues[37]
                      +wb*p2b*PointValues[38] + wa*p2a*PointValues[39] )*ny;
    Functionals[3] *= 0.5*Cell->GetNormalOrientation(1);
    Functionals[4] *= 0.5; // Cell->GetNormalOrientation(1);
    Functionals[5] *= 0.5*Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[6] = ( wa*PointValues[8]  + wb*PointValues[9]
                      +wb*PointValues[10] + wa*PointValues[11] )*nx
                    +( wa*PointValues[40] + wb*PointValues[41]
                      +wb*PointValues[42] + wa*PointValues[43] )*ny;
    Functionals[7] = (-wa*a*PointValues[8] - wb*b*PointValues[9]
                      +wb*b*PointValues[10]+ wa*a*PointValues[11] )*nx
                    +(-wa*a*PointValues[40] - wb*b*PointValues[41]
                      +wb*b*PointValues[42] + wa*a*PointValues[43] )*ny;
    Functionals[8] = ( wa*p2a*PointValues[8]  + wb*p2b*PointValues[9]
                      +wb*p2b*PointValues[10] + wa*p2a*PointValues[11] )*nx
                    +( wa*p2a*PointValues[40] + wb*p2b*PointValues[41]
                      +wb*p2b*PointValues[42] + wa*p2a*PointValues[43] )*ny;
    Functionals[6] *= 0.5*Cell->GetNormalOrientation(2);
    Functionals[7] *= 0.5; // Cell->GetNormalOrientation(2);
    Functionals[8] *= 0.5*Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[9] = ( wa*PointValues[12]  + wb*PointValues[13]
                      +wb*PointValues[14] + wa*PointValues[15] )*nx
                    +( wa*PointValues[44] + wb*PointValues[45]
                      +wb*PointValues[46] + wa*PointValues[47] )*ny;
    Functionals[10]= (-wa*a*PointValues[12] - wb*b*PointValues[13]
                      +wb*b*PointValues[14]+ wa*a*PointValues[15] )*nx
                    +(-wa*a*PointValues[44] - wb*b*PointValues[45]
                      +wb*b*PointValues[46] + wa*a*PointValues[47] )*ny;
    Functionals[11]= ( wa*p2a*PointValues[12]  + wb*p2b*PointValues[13]
                      +wb*p2b*PointValues[14] + wa*p2a*PointValues[15] )*nx
                    +( wa*p2a*PointValues[44] + wb*p2b*PointValues[45]
                      +wb*p2b*PointValues[46] + wa*p2a*PointValues[47] )*ny;
    Functionals[9] *= 0.5*Cell->GetNormalOrientation(3);
    Functionals[10]*= 0.5; // Cell->GetNormalOrientation(3);
    Functionals[11]*= 0.5*Cell->GetNormalOrientation(3);
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 0.25*Cell->GetMeasure();
    
    // we could use TQuadAffin if the Cell is a parallelogram, but this also
    // works for general quads
    TQuadBilinear referenceTransform;
    referenceTransform.SetCell(Cell);
    
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    // 16 quadrature points for 2D quadrature
    // the weights:
    double w[16] = { wa*wa, wa*wb, wa*wb, wa*wa,
                     wa*wb, wb*wb, wb*wb, wa*wb,
                     wa*wb, wb*wb, wb*wb, wa*wb,
                     wa*wa, wa*wb, wa*wb, wa*wa };
    
    // int_cell v . (1 0)^T 
    Functionals[12] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i],  eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[12] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[12] *= measure;
    
    // int_cell v . (0 1)^T 
    Functionals[13] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = 1., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[13] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[13] *= measure;
    
    // int_cell v . (x 0)^T
    Functionals[14] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = xi[16+i], uetaref = 0.;
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref, 
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[14] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[14] *= measure;
    
    
    // int_cell v . (0 x)^T
    Functionals[15] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = xi[16+i], uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref, 
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[15] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[15] *= measure;
    
    // int_cell v . (y 0)^T
    Functionals[16] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = eta[16+i], uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref, 
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[16] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[16] *= measure;
    
    // int_cell v . (0 y)^T
     Functionals[17] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = eta[16+i], uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref, 
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[17] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[17] *= measure;
    
    // int_cell v . (0 xx)^T
    Functionals[18] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = xi[16+i]*xi[16+i];
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[18] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[18] *= measure;
    
    // int_cell v . (yy 0)^T
    Functionals[19] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = eta[16+i]*eta[16+i], uetaref = 0.;
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[19] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[19] *= measure;
    
    // int_cell v . (xy 0)^T
    Functionals[20] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = xi[16+i]*eta[16+i], uetaref = 0.;
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[20] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[20] *= measure;
    
    // int_cell v . (0 xy)^T
    Functionals[21] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = xi[16+i]*eta[16+i];
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[21] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[21] *= measure;
    
    // int_cell v . (0 xxy)^T
    Functionals[22] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = xi[16+i]*xi[16+i]*eta[16+i];
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[22] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[22] *= measure;
    
    // int_cell v . (xyy 0)^T
    Functionals[23] = 0.;
    for(unsigned int i = 0; i < 16; ++i)
    {
      double uref = 0., uxiref = xi[16+i]*eta[16+i]*eta[16+i], uetaref = 0.;
      double uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(xi[16+i], eta[16+i], 1, &uref, &uxiref,
                                       &uetaref, &uorig, &uxorig, &uyorig);
      Functionals[23] += ( PointValues[16+i]*uxorig + PointValues[48+i]*uyorig )
                         * w[i];
    }
    Functionals[23] *= measure;
  }
}

void NF_N_Q_RT2_2D_EvalEdge(const TCollection *, const TBaseCell *Cell, int,
                            const double *PointValues, double *Functionals)
{
  double x0, x1, y0, y1, z; // z is just a dummy
  Cell->GetVertex(Joint)->GetCoords(x0, y0, z);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1, z); // 4=number of edges
  // length of joint, 0.5 due to 1D-reference cell having measure 2
  double l = 0.5*std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); 
  Functionals[0] = ( NF_N_Q_RT2_2D_wa*PointValues[0]
                    +NF_N_Q_RT2_2D_wb*PointValues[1]
                    +NF_N_Q_RT2_2D_wb*PointValues[2]
                    +NF_N_Q_RT2_2D_wa*PointValues[3] )*l;
  Functionals[1] = ( NF_N_Q_RT2_2D_wa*NF_N_Q_RT2_2D_T[0]*PointValues[0]
                    +NF_N_Q_RT2_2D_wb*NF_N_Q_RT2_2D_T[1]*PointValues[1]
                    +NF_N_Q_RT2_2D_wb*NF_N_Q_RT2_2D_T[2]*PointValues[2]
                    +NF_N_Q_RT2_2D_wa*NF_N_Q_RT2_2D_T[3]*PointValues[3] )*l;
  Functionals[2] = ( NF_N_Q_RT2_2D_wa*NF_N_Q_RT2_2D_p2a*PointValues[0]
                    +NF_N_Q_RT2_2D_wb*NF_N_Q_RT2_2D_p2b*PointValues[1]
                    +NF_N_Q_RT2_2D_wb*NF_N_Q_RT2_2D_p2b*PointValues[2]
                    +NF_N_Q_RT2_2D_wa*NF_N_Q_RT2_2D_p2a*PointValues[3] )*l;
}
