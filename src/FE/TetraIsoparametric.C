#include "TetraIsoparametric.h"
#include "BoundFace.h"
#include "FEDescriptor.h"
#include "IsoBoundFace.h"
#include "GridCell.h"
#include "MooNMD_Io.h"
#include "QuadFormula.h"
#include <cmath>

// constexpr BaseFunction_type BaseFunctFromOrder[] = { 
//                 BF_C_T_P0_3D, BF_C_T_P1_3D, BF_C_T_B2_3D, BF_C_T_P3_3D };
constexpr BaseFunction_type BaseFunctFromOrder[] = { 
                BF_C_T_P0_3D, BF_C_T_P1_3D, BF_C_T_P2_3D, BF_C_T_P3_3D };

// constexpr FEDescriptor_type FEDescFromOrder[] = { 
//                 FE_C_T_P0_3D, FE_C_T_P1_3D, FE_C_T_B2_3D, FE_C_T_P3_3D };
constexpr FEDescriptor_type FEDescFromOrder[] = { 
                FE_C_T_P0_3D, FE_C_T_P1_3D, FE_C_T_P2_3D, FE_C_T_P3_3D };

/** constuctor */
TTetraIsoparametric::TTetraIsoparametric()
{
}

/** transfer from reference element to original element */
void TTetraIsoparametric::GetOrigFromRef(double xi, double eta, double zeta,
                                        double &X, double &Y, double &Z) const
{
  int j = -1;
  unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
  for(unsigned int i=0;i<N_QuadPoints;i++)
  {
    auto p = quadrature_formula->get_point(i);
    if(std::abs(xi-p.x)<1e-8 && std::abs(eta-p.y)<1e-8
       && std::abs(zeta-p.z)<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    ErrThrow("error in TTetraIsoparametric::GetOrigFromRef(1)");
    return;
  }

  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
	    
  for(int i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
    Z += ZDistance[i] * FctValues[j][i];
  }
}

/** transfer a set of point from reference to original element */
void TTetraIsoparametric::GetOrigFromRef(int N_Points, const double *xi,
                                         const double *eta, const double *zeta,
                                         double *X, double *Y, double *Z) const
{
  int i, j, k;
  double Xi, Eta, Zeta;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double AuxVector[4*MaxN_BaseFunctions3D];

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];

    dx1 = xc1;
    dx2 = xc2;
    dx3 = xc3;
    
    dy1 = yc1;
    dy2 = yc2;
    dy3 = yc3;
    
    dz1 = zc1;
    dz2 = zc2;
    dz3 = zc3;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;

    if (N_AuxPoints)
    {
      BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
      bf.GetDerivatives(MultiIndex3D::D000, Xi, Eta, Zeta, AuxVector);
      bf.GetDerivatives(MultiIndex3D::D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf.GetDerivatives(MultiIndex3D::D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf.GetDerivatives(MultiIndex3D::D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
      
      for(k=0;k<N_AuxPoints;k++)
      {
        j = IntAux[k];
        X[i] += XDistance[k] * AuxVector[j];
        Y[i] += YDistance[k] * AuxVector[j];
        Z[i] += ZDistance[k] * AuxVector[j];
        
        dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
        
        dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
        
        dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
      }
    }
  }
}

void TTetraIsoparametric::GetOrigFromRef(const TQuadFormula& qf_ref,
                                         TQuadFormula& qf_orig) const
{
  unsigned int n_points = qf_ref.GetN_QuadPoints();
  if(qf_ref.get_type() != qf_orig.get_type())
    qf_orig = qf_ref; // copy
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double AuxVector[4*MaxN_BaseFunctions3D];
  for(unsigned int i=0;i<n_points;i++)
  {
    auto p = qf_ref.get_point(i);
    dx1 = xc1;
    dx2 = xc2;
    dx3 = xc3;
    
    dy1 = yc1;
    dy2 = yc2;
    dy3 = yc3;
    
    dz1 = zc1;
    dz2 = zc2;
    dz3 = zc3;
    
    double x = xc0 + xc1*p.x + xc2*p.y + xc3*p.z;
    double y = yc0 + yc1*p.x + yc2*p.y + yc3*p.z;
    double z = zc0 + zc1*p.x + zc2*p.y + zc3*p.z;

    if (N_AuxPoints)
    {
      BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
      bf.GetDerivatives(MultiIndex3D::D000, p.x, p.y, p.z, AuxVector);
      bf.GetDerivatives(MultiIndex3D::D100, p.x, p.y, p.z, AuxVector+MaxN_BaseFunctions3D);
      bf.GetDerivatives(MultiIndex3D::D010, p.x, p.y, p.z, AuxVector+2*MaxN_BaseFunctions3D);
      bf.GetDerivatives(MultiIndex3D::D001, p.x, p.y, p.z, AuxVector+3*MaxN_BaseFunctions3D);
      
      for(int k=0;k<N_AuxPoints;k++)
      {
        int j = IntAux[k];
        x += XDistance[k] * AuxVector[j];
        y += YDistance[k] * AuxVector[j];
        z += ZDistance[k] * AuxVector[j];
        
        dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
        
        dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
        
        dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
        dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
        dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
      }
    }
    double absdet = std::abs(dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                           - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2);
    double weight = qf_ref.get_weight(i) * absdet;
    qf_orig.update_pair(i, {weight, parmoon::Point{x, y, z}});
  }
}

/** transfer from original element to reference element */
void TTetraIsoparametric::GetRefFromOrig(double X, double Y, double Z,
                                         double &xi, double &eta, double &zeta)
  const
{

  Output::print("TTetraIsoparametric::GetRefFromOrig is not implemented");
  Output::print("TTetraIsoparametric::GetRefFromOrig works only for corners of cells");
  double xt = X - xc0;
  double yt = Y - yc0;
  double zt = Z - zc0;
  double recdetaffine;
  
  recdetaffine = 1/( xc1*yc2*zc3-xc1*yc3*zc2-yc1*xc2*zc3
                    +yc1*xc3*zc2+zc1*xc2*yc3-zc1*xc3*yc2 );

  xi = (-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
  eta = ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
  zeta = (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;
}

/** calculate functions and derivatives from reference element
    to original element */
void TTetraIsoparametric::GetOrigValues(double xi, double eta, double zeta,
                                        int N_BaseFunct, const double *uref,
                                        const double *uxiref,
                                        const double *uetaref,
                                        const double *uzetaref,
                                        double *uorig, double *uxorig,
                                        double *uyorig, double *uzorig,
                                        int BaseVectDim) const
{
  if(BaseVectDim != 1)
  {
    ErrThrow("TTetraIsoparametric for vector valued functions not "
             "implemented. (Piola transform)");
  }
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  // D000
  for(int i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  int j = -1;
  unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
  for(unsigned int i=0;i<N_QuadPoints;i++)
  {
    auto p = quadrature_formula->get_point(i);
    if(std::abs(xi-p.x)<1e-8 && std::abs(eta-p.y)<1e-8
       && std::abs(zeta-p.z)<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    ErrThrow(": error in TTetraIsoparametric::GetOrigValues(3)");
  }

  dx1 = xc1;
  dx2 = xc2;
  dx3 = xc3;
              
  dy1 = yc1;
  dy2 = yc2;
  dy3 = yc3;
  
  dz1 = zc1;
  dz2 = zc2;
  dz3 = zc3;

  for(int k=0;k<N_AuxPoints;k++)
  {
    dx1 += XDistance[k] * XiDerValues[j][k];
    dx2 += XDistance[k] * EtaDerValues[j][k];
    dx3 += XDistance[k] * ZetaDerValues[j][k];

    dy1 += YDistance[k] * XiDerValues[j][k];
    dy2 += YDistance[k] * EtaDerValues[j][k];
    dy3 += YDistance[k] * ZetaDerValues[j][k];

    dz1 += ZDistance[k] * XiDerValues[j][k];
    dz2 += ZDistance[k] * EtaDerValues[j][k];
    dz3 += ZDistance[k] * ZetaDerValues[j][k];
  }
  double detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  double rec_detjk = 1/detjk;

  // cout << "derivatives: " << detjk << endl;
  // cout << "5: " << dx1 << " " << dx2 << " " << dx3 << endl;
  // cout << "5: " << dy1 << " " << dy2 << " " << dy3 << endl;
  // cout << "5: " << dz1 << " " << dz2 << " " << dz3 << endl;

  // D100, D010 and D001
  for(int i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((dy2*dz3 - dy3*dz2)*uxiref[i] +
                 (dy3*dz1 - dy1*dz3)*uetaref[i] +
                 (dy1*dz2 - dy2*dz1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((dx3*dz2 - dx2*dz3)*uxiref[i] +
                 (dx1*dz3 - dx3*dz1)*uetaref[i] +
                 (dx2*dz1 - dx1*dz2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((dx2*dy3 - dx3*dy2)*uxiref[i] +
                 (dx3*dy1 - dx1*dy3)*uetaref[i] + 
                 (dx1*dy2 - dx2*dy1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void TTetraIsoparametric::GetOrigAllDerivatives(
    double xi, double eta, double zeta, int N_BaseFunct,
    const double* refD000, const double* refD100,
    const double* refD010, const double* refD001,
    const double* refD200, const double* refD110,
    const double* refD101, const double* refD020,
    const double* refD011, const double* refD002,
    double* origD000, double* origD100, double* origD010, double* origD001,
    double* origD200, double* origD110, double* origD101, double* origD020,
    double* origD011, double* origD002, int BaseVectDim) const
{
  GetOrigValues(xi, eta, zeta, N_BaseFunct, refD000, refD100, refD010, refD001,
                origD000, origD100, origD010, origD001, BaseVectDim);
  bool all_zero_2nd_derivatives = true;
  for(int i=0;i<N_BaseFunct;i++)
  {
    if(refD200[i] != 0. || refD110[i] != 0. || refD101[i] != 0. 
       || refD020[i] != 0. || refD011[i] != 0. || refD002[i] != 0.)
    {
      all_zero_2nd_derivatives = false;
      break;
    }
  } // endfor i
  if(!all_zero_2nd_derivatives)
    Output::warn("TTetraIsoparametric", "transformation of second order "
                 "derivatives is not yet implemented.");
  // avoid compiler warnings
  (void)origD200;
  (void)origD110;
  (void)origD101;
  (void)origD020;
  (void)origD011;
  (void)origD002;
}

void TTetraIsoparametric::SetCell(const TBaseCell *cell)
{
  int i, j, k, N_AuxVertices=0, N_IsoFaces=0;
  const TJoint *joint;
  JointType type;
  const TShapeDesc *ShapeDesc;
  const int *TmpFaceVertex, *TmpLen;
  int MaxLen;
  int *JointDOF, LocIsoDOF[4];
  double xm, ym, zm, xp, yp, zp;
  const TVertex * const * Vertices;
  const TVertex * const * AuxVertices;
  double x[3], y[3], z[3], factor;
  double LinComb[3][3];
  
  N_AuxPoints = 0;

  Vertices = ((const TGridCell*)cell)->GetVertices();
  ShapeDesc = cell->GetShapeDesc();
  ShapeDesc->GetFaceVertex(TmpFaceVertex, TmpLen, MaxLen);

  cell->GetVertex(0)->GetCoords(x0, y0, z0);
  cell->GetVertex(1)->GetCoords(x1, y1, z1);
  cell->GetVertex(2)->GetCoords(x2, y2, z2);
  cell->GetVertex(3)->GetCoords(x3, y3, z3);

  xc0 = x0;
  xc1 = x1-x0;
  xc2 = x2-x0;
  xc3 = x3-x0;

  yc0 = y0;
  yc1 = y1-y0;
  yc2 = y2-y0;
  yc3 = y3-y0;

  zc0 = z0;
  zc1 = z1-z0;
  zc2 = z2-z0;
  zc3 = z3-z0;

  N_IsoFaces=0;
  for(i=0;i<4;i++) // run through all faces
  {
    joint = cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryFace || type == IsoBoundFace)
     {
       FEDescriptor fedesc(FEDescFromOrder[ApproximationOrder]);
       JointDOF = fedesc.GetJointDOF(i);

      //TBoundFace *bdfacebdface = (TBoundFace*)joint;
//       bdface->GetParameters(Param1, Param2);

//       // do not apply moving of nodes on the bottom and the
//       // top of a sandwich mesh
//       if (bdface->GetBoundComp()->GetID()>=1000)
//         continue;

//       dt = ds = 1.0/ApproximationOrder;
      // coordinates of vertices of the face
      Vertices[TmpFaceVertex[i*MaxLen+0]]->GetCoords(x[0],y[0],z[0]);
      Vertices[TmpFaceVertex[i*MaxLen+1]]->GetCoords(x[1],y[1],z[1]);
      Vertices[TmpFaceVertex[i*MaxLen+2]]->GetCoords(x[2],y[2],z[2]);

       if(type == IsoBoundFace)
         {
           AuxVertices = ((TIsoBoundFace *)joint)->GetVertices();
           N_AuxVertices = ((TIsoBoundFace *)joint)->GetN_Vertices();
	   N_AuxVertices -= N_AuxVertices > 3 ? 1 : 0; // face bubble not needed
           N_IsoFaces++;
           if(N_IsoFaces>1)
           {
             Output::print("More than one isoface are in a sigle cell ", N_IsoFaces);
             Output::print("Check  TTetraIsoparametric::SetCell");
             Output::print("x[0]: ", x[0], " y[0]: ", y[0], " z[0]: ", z[0]);
             Output::print("x[1]: ", x[1], " y[1]: ", y[1], " z[1]: ", z[1]);
             Output::print("x[2]: ", x[2], " y[2]: ", y[2], " z[2]: ", z[2]);
             exit(0);
           }
         }

      if(ApproximationOrder==2 && N_AuxVertices==3)
      {
//      additional vertices should have been already added in this face in the main program
        LocIsoDOF[0] = 1;
        LocIsoDOF[1] = 3;
        LocIsoDOF[2] = 4;
//         LocIsoDOF[3] = 6;      // local dof of face point, see C_T_B2_3D_J

        LinComb[0][0]=0.5; LinComb[1][0]=0.5; LinComb[2][0]=0.0;
        LinComb[0][1]=0.5; LinComb[1][1]=0.0; LinComb[2][1]=0.5;
        LinComb[0][2]=0.0; LinComb[1][2]=0.5; LinComb[2][2]=0.5;

        for(j=0;j<N_AuxVertices;j++)
         {
          AuxVertices[j]->GetCoords(xp, yp, zp);
          xm = 0.;  ym=0.; zm=0.;
          for(k=0;k<3;k++)    // number of vertices on the face
            {
              factor = LinComb[j][k];
              xm += factor*x[k];
              ym += factor*y[k];
              zm += factor*z[k];
            }

          if(std::abs(xp-xm) > 1e-8 || std::abs(yp-ym) > 1e-8 || std::abs(zp-zm) > 1e-8)
           {

//             cout << endl;
//             cout << "xpzp " << xp*xp + yp*yp +  zp*zp<< endl;
//             cout << "xmzm " << xm*xm + ym*ym +  zm*zm<< endl;
//             Output::print("xp: ", xp, " xm: ", xm, " X: ", xp-xm);
//             Output::print("yp: ", yp, " ym: ", ym, " Y: ", yp-ym);
//             Output::print("zp: ", zp, " zm: ", zm, " Z: ", zp-zm);
//             Output::print("R0: ", std::sqrt((xp-xm)*(xp-xm) + (yp-ym)*(yp-ym)+ (zp-zm)*(zp-zm)));
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            ZDistance[N_AuxPoints] = zp - zm;
            IntAux[N_AuxPoints] = JointDOF[LocIsoDOF[j]];
            N_AuxPoints++;
           }
         } //    for(j=0;j<=N_A
     }
    else if(ApproximationOrder>1)
     {
       Output::print(" ApproximationOrder ", ApproximationOrder, " ",
                     N_AuxVertices);
      ErrThrow("Tetra isoparametric for approximation order other than 2 has "
               "to be implemented !!! set USE_ISOPARAMETRIC: 0 ");
      }
    } // endif
   } // endfor i
//  Output::print(N_AuxPoints, " : N_AuxPoints, ", " N_QuadPoints  :", N_QuadPoints);

///@todo There was a special routine in MooNMD for isoparametric finite elements of the 
/// standard benchmark problem 'Flow around a circular cylinder in 3d' for second and
/// third order tetrahedra. This routine should be available again in future. 

//   if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY==1)
//   {
//     /***************************************************************************/
//     // this is for isoparametric tetrahedral finite elements in the 3d circular
//     // benchmark
//     // one has to change in addition in FEDatabase3D.C , ~ line 748
//    cout<< "Check tetraisoparametric.C !!!!!!!!!!!!!!!!!!!!! "<<endl;
// 
//     exit(0);
//     // no boundary face
//     if(N_AuxPoints == 0)
//     {
//       // for all vertices    
//       for(i=0;i<4;i++)
//       {
//         // get coordinates of vertex i
//         cell->GetVertex(i)->GetCoords(xp, yp, zp);
//         // point on cylinder
//         if(std::abs((xp-0.5)*(xp-0.5)+(yp-0.2)*(yp-0.2)-0.0025)<1e-10)
//           //if(std::abs(xp*xp+yp*yp+zp*zp-1)<1e-10)
//         {
//           // point on boundary
//           // look to all other vertices
//           for(j=i+1;j<4;j++)
//           {
//             // get coordinates of second vertex
//             cell->GetVertex(j)->GetCoords(xm, ym, zm);
//             // if second vertex on cylinder
//             if(std::abs((xm-0.5)*(xm-0.5)+(ym-0.2)*(ym-0.2)-0.0025)<1e-10)
//               //if(std::abs(xm*xm+ym*ym+zm*zm-1)<1e-10)
//             {
//               // point on boundary
//               // find number of dof for P_2
//               bf = FEDatabase::GetBaseFunct3D(
//                 BaseFunctFromOrder[ApproximationOrder]);
//               if (ApproximationOrder==2)
//               {
//                 // switch number of first vertex
//                 switch(i)
//                 {
//                   case 0:
//                     // switch number of second vertex
//                     switch(j)
//                     {
//                       case 1: m = 1; break;
//                       case 2: m = 3; break;
//                       case 3: m = 6; break;
//                     } // endswitch(j)
//                     break;
//                     
//                   case 1:
//                     switch(j)
//                     {
//                       case 2: m = 4; break;
//                       case 3: m = 7; break;
//                     } // endswitch(j)
//                     break;
//                     
//                   case 2:
//                     m = 8;
//                     break;
//                 } // endswitch(i)
//                 
//                 // centre of edge
//                 a = 0.5*(xp+xm);
//                 b = 0.5*(yp+ym);
//                 c = 0.5*(zp+zm);
//                 // projection to the boundary
//                 T = std::sqrt((a-0.5)*(a-0.5)+(b-0.2)*(b-0.2));
//                 project_a = 0.5 + 0.05*(a-0.5)/T;
//                 project_b = 0.2 + 0.05*(b-0.2)/T;    
//                 project_c = c;
//                 // T = std::sqrt(a*a+b*b+c*c);
//                 // project_a = a/T;
//                 // project_b = b/T;    
//                 //project_c = c/T;    
//                 if(std::abs(a-project_a) > 1e-8 || std::abs(b-project_b) > 1e-8  || std::abs(c-project_c) > 1e-8)
//                 {
//                   // differences
//                   XDistance[N_AuxPoints] = project_a - a;
//                   YDistance[N_AuxPoints] = project_b - b;
//                   ZDistance[N_AuxPoints] = project_c - c;
//                   // dof which is moved    
//                   IntAux[N_AuxPoints] = m;
//                   N_AuxPoints++;
//                 }
//               }
//               // P_3
//               if (ApproximationOrder==3)
//               {
//                 switch(i)
//                 {
//                   case 0:
//                     switch(j)
//                     {
//                       case 1: m = 1; n = 2; break;
//                       case 2: m = 4; n = 7; break;
//                       case 3: m = 10; n = 16;  break;
//                     } // endswitch(j)
//                     break;
//                     
//                   case 1:
//                     switch(j)
//                     {
//                       case 2: m = 6; n = 8; break;
//                       case 3: m = 12; n = 17; break;
//                     } // endswitch(j)
//                     break;
//                     
//                   case 2:
//                     m = 15; n = 18;
//                     break;
//                 } // endswitch(i)
//                 
//                 //  first dof
//                 a = (2*xp+xm)/3;
//                 b = (2*yp+ym)/3;
//                 // projection to the boundary
//                 T = std::sqrt((a-0.5)*(a-0.5)+(b-0.2)*(b-0.2));
//                 project_a = 0.5 + 0.05*(a-0.5)/T;
//                 project_b = 0.2 + 0.05*(b-0.2)/T;         
//                 // Output::print("0: a ", a, " project a ", project_a, " b ", b,
//                 //       " project b ", project_b, " ",
//                 //       std::sqrt((project_a-0.5)*(project_a-0.5)+(project_b-0.2)*(project_b-0.2)));
//                 // differences
//                 XDistance[N_AuxPoints] = project_a - a;
//                 YDistance[N_AuxPoints] = project_b - b;
//                 ZDistance[N_AuxPoints] = 0;
//                 IntAux[N_AuxPoints] = m;
//                 N_AuxPoints++;
//                 //  second dof
//                 a = (xp+2*xm)/3;
//                 b = (yp+2*ym)/3;
//                 // projection to the boundary
//                 project_a = 0.5 + 0.05*(a-0.5)/T;
//                 project_b = 0.2 + 0.05*(b-0.2)/T;         
//                 //Output::print("1: a ", a, " project a ", project_a, " b ", b,
//                 //       " project b ", project_b, " ",
//                 //      std::sqrt((project_a-0.5)*(project_a-0.5)+(project_b-0.2)*(project_b-0.2)));
//                 // differences
//                 XDistance[N_AuxPoints] = project_a - a;
//                 YDistance[N_AuxPoints] = project_b - b;
//                 ZDistance[N_AuxPoints] = 0;
//                 IntAux[N_AuxPoints] = n;
//                 N_AuxPoints++;
//               }
//               
//             } // endif
//           } // endfor j
//         } // endif
//       } // endfor i
//     } // endif
//   }
 

  /*if (N_AuxPoints)
  {
    Output::print("auxa ", N_AuxPoints);
    Output::print(x0, " ", y0, " ", z0);
    Output::print(x1, " ", y1, " ", z1);
    Output::print(x2, " ", y2, " ", z2);
    Output::print(x3, " ", y3, " ", z3);
    }*/

//  Output::print(N_AuxPoints, " : N_AuxPoints, ",  " N_QuadPoints  >", N_QuadPoints);
  if(N_AuxPoints)
  {
    BaseFunctions bf(BaseFunctFromOrder[ApproximationOrder]);
    unsigned int N_QuadPoints = quadrature_formula->GetN_QuadPoints();
    for(unsigned int i=0;i<N_QuadPoints;i++)
    {
      auto p = quadrature_formula->get_point(i);
      bf.GetDerivatives(MultiIndex3D::D000, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]];

      bf.GetDerivatives(MultiIndex3D::D100, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]];

      bf.GetDerivatives(MultiIndex3D::D010, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]];

      bf.GetDerivatives(MultiIndex3D::D001, p.x, p.y, p.z, DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        ZetaDerValues[i][j] = DoubleAux[IntAux[j]];
    } // endfor i
  } // endif
}

void TTetraIsoparametric::SetQuadFormula(QuadratureFormula_type formula)
{
  QuadFormula_type = formula;
  quadrature_formula = std::unique_ptr<TQuadFormula>(new TQuadFormula(formula));
}

void TTetraIsoparametric::PiolaMapOrigFromRef(double, double, double, int,
                                              const double *, double *) const
{
  ErrThrow("TTetraIsoparametric::PiolaMapOrigFromRef not implemented");
}
