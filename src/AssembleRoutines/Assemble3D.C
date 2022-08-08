
// =======================================================================
// %W% %G%
// 
// Class:       TAssemble3D
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <MooNMD_Io.h>
#include <Assemble3D.h>
#include <Enumerations_fe.h>
#include <Matrix3D.h>
#include <AuxParam3D.h>
#include <Joint.h>
#include <BoundFace.h>
#include "FEDatabase.h"
#include <SquareMatrix3D.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <Database.h>
#include <Convolution.h>
#include <Edge.h>
#include "MainUtilities.h"
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "QuadratureFormulaDatabase.h"



void Assemble3D(int n_fespaces, const TFESpace3D** fespaces, int n_sqmatrices,
                TSquareMatrix3D** sqmatrices, int n_matrices,
                TMatrix3D** matrices, int n_rhs, double** rhs,
                const TFESpace3D** ferhs, BoundCondFunct3D** BoundaryConditions,
                BoundValueFunct3D** BoundaryValues, LocalAssembling3D& la,
                bool assemble_dirichlet_rows)
{
//   if(n_rhs != la.get_n_rhs())
//   {
//     ErrThrow("the number of right-hand-sides in Assemble2D does not match that "
//              "in the LocalAssembling2D object, ", n_rhs, " != ",
//              la.get_n_rhs());
//   }
  
  std::vector<const TFESpace3D*> further_fe_spaces;
  for(size_t i = 0, n = la.n_fe_functions(); i < n; ++i)
  {
    // space used in local assembling object
    const TFESpace3D* space = la.get_fe_function(i)->GetFESpace3D().get();
    // check if this space is not one of the `fespaces`
    auto first = fespaces, last = fespaces + n_fespaces;
    while(first != last)
    {
      if(*first == space) break;
      ++first;
    }
    if(first == last)
      further_fe_spaces.push_back(space);
  }
  
  std::vector< const FiniteElement*> LocalUsedElements(n_fespaces + further_fe_spaces.size(),
                                                       nullptr);
          int N_AllMatrices = n_sqmatrices+n_matrices;
          int i,j,k,l,l1,l2,l3,n,m,ij,N_Vertex;
          int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
          int N_Test, N_Ansatz, N_Joints;
          const TQuadFormula *qf2;
          TJoint *joint;
          //  TIsoBoundEdge *isoboundedge;
          const double *xi, *eta, *zeta;
          const double *t, *s;
          double xf, yf, zf;
          double *Param[MaxN_QuadPoints_3D];
          double *local_rhs =nullptr;
          double *righthand =nullptr;
          double **Matrices =nullptr, *aux =nullptr;
          double **Matrix=nullptr;
          double ***LocMatrices=nullptr, **LocRhs=nullptr;
          std::vector<const BaseFunctions*> LocBF(n_fespaces + further_fe_spaces.size());
          int ActiveBound, DirichletBound, end;
          const int *TestDOF, *AnsatzDOF;
          double *Entries;
          const int *ColInd, *RowPtr;
          double *RHS, *MatrixRow;
          double t0, t1, t2;
          BoundCond Cond0;
          BoundCondFunct3D *BoundaryCondition;
          BoundValueFunct3D *BoundaryValue;
          ReferenceTransformation_type reftrans;
          double PointValues[MaxN_PointsForNodal3D];
          double FunctionalValues[MaxN_BaseFunctions3D];
          int *EdgeDOF, N_EdgeDOF;
          double **JointValues, *JointValue;
          // static bool *SecondDer = nullptr;
          bool *SecondDer = nullptr;
          double LinComb[4];

          const int *TmpFV, *TmpLen, *TmpFE, *ETmpLen;
          int MaxLen, EMaxLen;
          double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
          double nx, ny, nz;

          double time1, time2, time_all, time_total;

#ifdef _MPI
  const int *EdgeVertex;
  int N_FaceEdges;
  int N_Edges, N_VertInCell, dof;
  TEdge *edge;
  TVertex *Vert;
#endif

          bool InnerBoundary, OuterBoundary;

          time_total = GetTime();
          time_all = 0;

        // ########################################################################
        // store information in local arrays
        // ########################################################################
          if(n_rhs)
          {
            LocRhs = new double* [n_rhs];
            righthand = new double [n_rhs*MaxN_BaseFunctions3D];
            for(i=0;i<n_rhs;i++)
              LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;

          } // endif n_rhs

          N_Parameters = la.GetN_Parameters(); //ask local assembling object for number of parameters

          if(N_Parameters)
          {
            aux = new double [MaxN_QuadPoints_3D*N_Parameters];
            for(j=0;j<MaxN_QuadPoints_3D;j++)
              Param[j] = aux + j*N_Parameters;
          }

          if(N_AllMatrices)
          {
            aux = new double
                    [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
            Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
            for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
              Matrices[j] = aux+j*MaxN_BaseFunctions3D;

            LocMatrices = new double** [N_AllMatrices];
            for(i=0;i<N_AllMatrices;i++)
              LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
          } // endif N_AllMatrices

          SecondDer = la.GetNeeds2ndDerivatives();

          if(SecondDer == nullptr)
          {
            SecondDer = new bool[n_fespaces];
            for(i=0;i<n_fespaces;i++)
              SecondDer[i] = false;
          }


          // all spaces use same Coll
          auto Coll = fespaces[0]->GetCollection();
          N_Cells = Coll->GetN_Cells();
        //   cout << "N_Cells: " << N_Cells << endl;

          Coll->mark_all_cells();
          
          TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
          TQuadFormula qf_orig(qf_ref);

        // ########################################################################
        // loop over all cells
        // ########################################################################
          for(int i=0;i<N_Cells;i++)
          {
            auto cell = Coll->GetCell(i);

            // ####################################################################
            // find local used elements on this cell
            // ####################################################################
            for(int j=0;j<n_fespaces;j++)
            {
              auto& element = fespaces[j]->get_fe(i);
              LocalUsedElements[j] = &element;
              LocBF[j] = element.GetBaseFunct();
            }
            for(auto j = 0u; j < further_fe_spaces.size(); ++j)
            {
              auto& element = further_fe_spaces[j]->get_fe(i);
              LocalUsedElements[j+n_fespaces] = &element;
              LocBF[j+n_fespaces] = element.GetBaseFunct();
            }

            // ####################################################################
            // calculate values on original element
            // ####################################################################

            //Output::print("CELL ", i);
            reftrans = FEDatabase::GetOrig(LocalUsedElements, Coll, cell,
                                           SecondDer, qf_ref, qf_orig);

            // TODO: local cell vertices etc. should really be communicated more tidily

            N_Vertex = cell->GetN_Vertices();

            for (ij = 0; ij < N_Vertex; ij++)
            {
              TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
              TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
              TDatabase::ParamDB->INTERNAL_VERTEX_Z[ij] = cell->GetVertex(ij)->GetZ();
            }

            if (N_Vertex == 4)
            {
              TDatabase::ParamDB->INTERNAL_VERTEX_X[4] = -4711;
            }

            TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;

            // use local assembling object to assemble a few matrices and
            // right-hand sides at once
            la.GetLocalForms(qf_orig, LocBF, cell, i, N_AllMatrices, n_rhs,
                             LocMatrices, LocRhs);

            //Output::print("local form done ", i);
            time1 = GetTime();
            // ####################################################################
            // add local matrices to global matrices (ansatz == test)
            // ####################################################################
            for(j=0;j<n_sqmatrices;j++)
            {
              if(sqmatrices[j] == nullptr)
                continue;
             // find space for this bilinear form
              auto fespace = sqmatrices[j]->GetFESpace3D();
              auto element = fespace->get_fe(i);
              N_ = element.GetN_DOF();

              Matrix = Matrices+j*MaxN_BaseFunctions3D;
              Entries = sqmatrices[j]->GetEntries();
              RowPtr = sqmatrices[j]->get_row_ptr();
              ColInd = sqmatrices[j]->get_vector_columns();

              ActiveBound = fespace->get_n_active();
              DirichletBound = fespace->get_n_active_non_hanging();
              const int * DOF = fespace->GetGlobalDOF(i);

              // add local matrix to global
              for(m=0;m<N_;m++)
              {
                l=DOF[m];
                MatrixRow = Matrix[m];
                // cout << "DOF: " << l << endl;
                if(l<DirichletBound || assemble_dirichlet_rows)
                {
                  // node l is inner or Neumann node
                  // for all dof
                  for(k=0;k<N_;k++)
                  {
                    // for all columns
                    l2 = 0;
                    end=RowPtr[l+1];
                    for(n=RowPtr[l];n<end;n++)
                    {
                      // if column of current dof found
                      if(DOF[k] == ColInd[n])
                      {
                        // cout << m << "   " << k << endl << n << endl;
                        Entries[n] += MatrixRow[k];
                        l2 = 1;
                        break;
                      } // endif
                    } // endfor n
                    if(l2 == 0)
                      cout << "not found it's me" << endl;
                  } // endfor k
                } // endif l
                else
                {
                  // Dirichlet node
                  sqmatrices[j]->set(l, l, 1.);
                }
              } // endfor m
            } // endfor j


            // ####################################################################
            // add local matrices to global matrices (ansatz != test)
            // ####################################################################
            for(j=0;j<n_matrices;j++)
            {
              if(matrices[j] == nullptr)
                continue;
              auto test_space = matrices[j]->GetTestSpace3D();
              auto ansatz_space = matrices[j]->GetAnsatzSpace3D();
              auto TestElement = test_space->get_fe(i);
              auto AnsatzElement = ansatz_space->get_fe(i);

              // cout << "non square matrix: " << j << endl;
              // cout << "TestElement: " << TestElement.GetID() << endl;
              // cout << "AnsatzElement: " << AnsatzElement.GetID() << endl;

              N_Test = TestElement.GetN_DOF();
              N_Ansatz = AnsatzElement.GetN_DOF();

              Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions3D;

              Entries = matrices[j]->GetEntries();
              RowPtr = matrices[j]->get_row_ptr();
              ColInd = matrices[j]->get_vector_columns();

              TestDOF = test_space->GetGlobalDOF(i);
              AnsatzDOF = ansatz_space->GetGlobalDOF(i);

              // add local matrix to global
              for(m=0;m<N_Test;m++)
              {
                l=TestDOF[m];
                MatrixRow = Matrix[m];
                // cout << "DOF: " << l << endl;
                for(k=0;k<N_Ansatz;k++)
                {
                  end = RowPtr[l+1];
                  for(n=RowPtr[l];n<end;n++)
                  {
                    if(AnsatzDOF[k] == ColInd[n])
                    {
                      // cout << m << "   " << k << endl << n << endl;
                      Entries[n] += MatrixRow[k];
                      break;
                    } // endif
                  } // endfor n
                } // endfor k
              } // endfor m
            } // endfor j
            time2 = GetTime();
            time_all += time2-time1;

            // ####################################################################
            // add local right-hand sides to global right-hand side
            // ####################################################################
            for(j=0;j<n_rhs;j++)
            {
              auto fespace = ferhs[j];
              ActiveBound = fespace->get_n_active();
              auto element = fespace->get_fe(i);

              N_ = element.GetN_DOF();

              local_rhs = righthand+j*MaxN_BaseFunctions3D;
              RHS = rhs[j];
              if(RHS == nullptr)
              {
                continue;
              }
              // find space for this linear form

              ActiveBound = fespace->get_n_active();
              DirichletBound = fespace->get_n_active_non_hanging();
              const int * DOF = fespace->GetGlobalDOF(i);

              // add local right-hand side to the global one
              for(m=0;m<N_;m++)
              {
                l=DOF[m];
                // cout << "DOF: " << l << endl;
                if(l<DirichletBound)
                {
                  // node l is inner or Neumann node
                  RHS[l] += local_rhs[m];
                } // endif l
              } // endfor m

              BoundaryCondition = BoundaryConditions[j];
              BoundaryValue = BoundaryValues[j];
              auto nf = element.GetNodalFunctional();

              auto FEDesc_Obj = element.GetFEDesc();
              N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

              // setting Dirichlet boundary condition
              N_Joints = cell->GetN_Faces();
              for(m=0;m<N_Joints;m++)
              {
                joint = cell->GetJoint(m);
                InnerBoundary = false;
                OuterBoundary = false;

                if(joint->GetType() == BoundaryFace ||
                   joint->GetType() == IsoBoundFace)
                  OuterBoundary = true;

                /*
                // check whether neighbour does not belong to Coll
                neigh = joint->GetNeighbour(cell);
                if(neigh)
                {
                  // check for neighbour's clipboard
                  if(neigh->GetClipBoard() == -1)
                    InnerBoundary = true;
                }
                */

                if(InnerBoundary || OuterBoundary)
                {
                  cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
                  cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
                  t0 = 1.0/TmpLen[m];

                  xf = 0; yf = 0; zf = 0;
                  double X[4], Y[4], Z[4];
                  for (l1=0;l1<TmpLen[m];l1++)
                  {
                    cell->GetVertex(TmpFV[m*MaxLen+l1])
                                ->GetCoords(X[l1], Y[l1], Z[l1]);
                    //LinComb[l1] = t0;
                    xf += t0*X[l1];
                    yf += t0*Y[l1];
                    zf += t0*Z[l1];
                  }

                  auto boundface              = (TBoundFace*)joint;
                  const TBoundComp* BoundComp = boundface->GetBoundComp();
                  int comp                    = BoundComp->get_physical_id();

                  /*if(OuterBoundary)
                  {
                    boundface = (TBoundFace *)joint;
                    BoundComp = boundface->GetBoundComp();

                    boundface->GetParameters(Param1, Param2);
                    comp=BoundComp->GetID();

                    BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                           X, Y, Z, Param1, Param2,
                                           xf, yf, zf, t0, t1);
                                           }*/
                  // the face gets the b.c. which is valid at its center
//                   BoundaryCondition(xf, yf, zf, Cond0);
                  BoundaryCondition(comp, xf, yf, zf, Cond0);

            /**
             @brief change the Cond0 to NEUMANN according to
             the markers prescribed in the input file
             @attention this part is still on-going work
             **/
            /*TBoundFace *boundface = (TBoundFace *)joint;
            int face_marker = boundface->GetBoundComp()->get_physical_id();
            for (int ibd=0; ibd< TDatabase::ParamDB->n_neumann_boundary; ibd++)
            {
                if (face_marker==TDatabase::ParamDB->neumann_boundary_id[ibd])
                    Cond0 = NEUMANN;
            }

            for (int ibd=0; ibd< TDatabase::ParamDB->n_nitsche_boundary; ibd++)
            {
                if (face_marker==TDatabase::ParamDB->nitsche_boundary_id[ibd])
                    Cond0 = NEUMANN;
            }

            */

                  switch(Cond0)
                  {
                    case DIRICHLET:
                      if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
                      {
                        nf->GetPointsForFace(N_Points, t, s);

                        for(l1=0;l1<N_Points;l1++)
                        {
                          switch(TmpLen[m])
                          {
                            case 4:
                              LinComb[0] = (1-t[l1])*(1-s[l1]);
                              LinComb[1] = t[l1]*(1-s[l1]);
                              LinComb[2] = t[l1]*s[l1];
                              LinComb[3] = (1-t[l1])*s[l1];
                              xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                                  +LinComb[2]*X[2] + LinComb[3]*X[3];
                              yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                                  +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                              zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                                  +LinComb[2]*Z[2] + LinComb[3]*Z[3];
                            break;
                            case 3:
                              LinComb[0] = 1-t[l1]-s[l1];
                              LinComb[1] = t[l1];
                              LinComb[2] = s[l1];
                              xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                                  +LinComb[2]*X[2];
                              yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                                  +LinComb[2]*Y[2];
                              zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                                  +LinComb[2]*Z[2];
                            break;
                          }
                          /*if(OuterBoundary)
                            BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                                   X, Y, Z, Param1, Param2,
                                                   xf, yf, zf, t0, t1);
                          */
                          //cout << l1 << ": " << t0 << " " << t1 << " <> ";
                          //cout << xf << " " << yf << " " << zf << endl;

                          BoundaryValue(comp, xf, yf, zf, PointValues[l1]);
                          // cout << " PV: " << PointValues[l1] << endl;
                        }
                      }
                      else
                      {
                        // no isoparametric
                        nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
                        std::vector<double> XYZ(3*N_Points);

                        FEDatabase::GetOrigFromRef(reftrans, N_Points,
                                                   xi, eta, zeta, XYZ.data(),
                                                   XYZ.data()+N_Points,
                                                   XYZ.data()+2*N_Points);

                        for(l1=0;l1<N_Points;l1++)
                        {
                          BoundaryValue(comp, XYZ[l1], XYZ[l1+N_Points],
                                        XYZ[l1+2*N_Points], PointValues[l1]);
                          // PointValues[l1] = 0;
                          // cout << " PV: " << PointValues[l1] << endl;
                        }
                      }

                      nf->GetFaceFunctionals(Coll, cell, m, PointValues,
                                             FunctionalValues);
                      EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
                      for(l=0;l<N_EdgeDOF;l++)
                      {
                        RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                        // cout << " DOF: " << DOF[EdgeDOF[l]];
                        // cout << " FV: " << FunctionalValues[l] << endl;
                      }
                    break;

                    case NEUMANN:
                    case ROBIN:  // lhs of robin will be assembled in main program
                    {
                      // cout << "Neumann condition in Assemble3D" << endl;
                      const BaseFunctions* bf = element.GetBaseFunct();
                      l = bf->GetPolynomialDegree();
                      const Shapes* face_types;
                      cell->GetShapeDesc()->GetFaceType(face_types);
                      qf2 = QuadratureFormulaDatabase::qf_from_degree(
                          2*l, face_types[m]);
                      N_Points = qf2->GetN_QuadPoints();
                      // values of base functions in all quadrature points on face
                      JointValues = FEDatabase::GetJointDerivatives3D(
                        *bf, *qf2, m, MultiIndex3D::D000);
                      bf->ChangeBF(Coll, cell, N_Points, JointValues);

                      switch(TmpLen[m])
                      {
                        case 3:
                          xc1 = X[1] - X[0];
                          xc2 = X[2] - X[0];

                          yc1 = Y[1] - Y[0];
                          yc2 = Y[2] - Y[0];

                          zc1 = Z[1] - Z[0];
                          zc2 = Z[2] - Z[0];

                          // normal vector
                          nx = yc1*zc2 - zc1*yc2;
                          ny = zc1*xc2 - xc1*zc2;
                          nz = xc1*yc2 - yc1*xc2;
                          // determinant of reference trafo
                          t2 = std::sqrt(nx*nx + ny*ny + nz*nz);

                          for(l=0;l<N_Points;l++)
                          {
                            JointValue = JointValues[l];
                            auto p = qf2->get_point(l);
                            t0 = p.x;
                            t1 = p.y;
                            // cout << "t: " << t0 << " " << t1 << endl;
                            LinComb[0] = 1-t0-t1;
                            LinComb[1] = t0;
                            LinComb[2] = t1;

                            xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                                +LinComb[2]*X[2];
                            yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                                +LinComb[2]*Y[2];
                            zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                                +LinComb[2]*Z[2];

                            /*if(OuterBoundary)
                              BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                                     X, Y, Z, Param1, Param2,
                                                     xf, yf, zf, t0, t1);
                            */
                            
                            BoundaryValue(comp, xf, yf, zf, t0);
                            // cout << t1 << endl;
                            t0 *= qf2->get_weight(l)*t2;
                            for(k=0;k<N_;k++)
                              if((l3 = DOF[k])<ActiveBound)
                              {
                                RHS[l3] += t0*JointValue[k];
                              }
                          } // endfor l
                        break;

                        case 4:
                          xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
                          xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
                          xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

                          yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
                          yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
                          yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

                          zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
                          zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
                          zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

                          for(l=0;l<N_Points;l++)
                          {
                            JointValue = JointValues[l];
                            auto p = qf2->get_point(l);
                            t0 = 0.5*(p.x+1);
                            t1 = 0.5*(p.y+1);
                            // cout << "t: " << t0 << " " << t1 << endl;
                            LinComb[0] = (1-t0)*(1-t1);
                            LinComb[1] = t0*(1-t1);
                            LinComb[2] = t0*t1;
                            LinComb[3] = (1-t0)*t1;

                            xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                                +LinComb[2]*X[2] + LinComb[3]*X[3];
                            yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                                +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                            zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                                +LinComb[2]*Z[2] + LinComb[3]*Z[3];

                            /*if(OuterBoundary)
                              BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                                     X, Y, Z, Param1, Param2,
                                                     xf, yf, zf, t0, t1);
                            */
                            // cout << xf << " " << yf << " " << zf << endl;
                            BoundaryValue(comp, xf, yf, zf, t0);
                            // cout << "PV: " << t0 << endl;
                            nx = (yc1 + p.y * yc3)*(zc2 + p.x * zc3)
                                -(zc1 + p.y * zc3)*(yc2 + p.x * yc3);
                            ny = (zc1 + p.y * zc3)*(xc2 + p.x * xc3)
                                -(xc1 + p.y * xc3)*(zc2 + p.x * zc3);
                            nz = (xc1 + p.y * xc3)*(yc2 + p.x * yc3)
                                -(yc1 + p.y * yc3)*(xc2 + p.x * xc3);
                            t1 = nx*nx+ny*ny+nz*nz;
                            // cout << t1 << endl;
                            t0 *= qf2->get_weight(l)*std::sqrt(t1);
                            for(k=0;k<N_;k++)
                              if((l3 = DOF[k])<ActiveBound)
                                RHS[l3] += t0*JointValue[k];
                          } // endfor l
                        break;
                      }
                      break;
                    }
                    default:

                     break;


                  } // endswitch Cond0

        #ifdef _MPI
                /** edges on this face are already set, so no need of checking in BD edges on this face */
                N_FaceEdges = ETmpLen[m];
                for(l1=0;l1<N_FaceEdges;l1++)
                 {
                  edge = cell->GetEdge(TmpFE[m*EMaxLen+l1]);
                  edge->SetClipBoard(i);
                 } // fo

                 /** vertices on this face are already set, so no need of checking in BD vert on this face, 30.06.12, sashi */
                 for (l1=0;l1<TmpLen[m];l1++)
                  cell->GetVertex(TmpFV[m*MaxLen+l1])->SetClipBoard(i);
        #endif
                } // endif
              } // endfor m

        #ifdef _MPI
            if( (N_EdgeDOF = FEDesc_Obj->GetN_EdgeDOF()) > 0) //conforming FE
             {
              cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex);
              N_Edges=cell->GetN_Edges();

              
              for(m=0;m<N_Edges;m++)
               {
                edge = cell->GetEdge(m);

                if( (edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D) && edge->GetClipBoard()==-1)
                 {
                  // BD edge is not yet set
                  edge->SetClipBoard(i);
                  EdgeDOF = FEDesc_Obj->GetEdgeDOF(m);
                  int comp = edge->get_physical_id();
                  for(l=0;l<N_EdgeDOF;l++)
                   {
                    dof = DOF[EdgeDOF[l]];
                    fespace->GetDOFPosition(dof, xf, yf, zf);
//                     BoundaryCondition(xf, yf, zf, Cond0);
                    BoundaryCondition(comp, xf, yf, zf, Cond0);

                    if(Cond0==DIRICHLET) // nothing to do for non-Dirichlet
                     {
                     BoundaryValue(comp, xf, yf, zf, t0);
                     RHS[dof] = t0;
                     }
                  }
                  /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */
                  cell->GetVertex(EdgeVertex[2*m])->SetClipBoard(i);
                  cell->GetVertex(EdgeVertex[2*m+1])->SetClipBoard(i);
                 }
               }//   endfor m
             }
            /** correct the Vert Dirichlet   */
            if( FEDesc_Obj->GetN_VertDOF() > 0) //conforming FE
             {
              N_VertInCell = cell->GetN_Vertices();

              for(m=0;m<N_VertInCell;m++)
               {
                Vert = cell->GetVertex(m);

                if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
                 {
                  // BD vert is not yet set
                  Vert->SetClipBoard(i);
                  Vert->GetCoords(xf, yf, zf);
//                   BoundaryCondition(xf, yf, zf, Cond0);
                  int comp = Vert->get_physical_id();
                  BoundaryCondition(comp, xf, yf, zf, Cond0);

                  if(Cond0==DIRICHLET) // nothing to do for nin-Dirichlet
                   {
                    dof = DOF[FEDesc_Obj->GetVertDOF(m)];
                    BoundaryValue(comp, xf, yf, zf, t0);
                    RHS[dof] = t0;
                   }

                 }
               }//   endfor m
              } //   if( FEDesc_Obj->GetN_VertDOF() > 0)
        #endif

            } // endfor j
            // cout << "end i:" << i << endl;
          } // endfor i

          // ####################################################################
          // modify matrix according to coupling
          // ####################################################################
          for(j=0;j<n_sqmatrices;j++)
          {
            if(sqmatrices[j] == nullptr)
              continue;
            if(!assemble_dirichlet_rows)
              sqmatrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);

          } // endfor j

          for(j=0;j<n_rhs;j++)
          {
            auto fespace = ferhs[j];
            N_Hanging = fespace->get_n_hanging();
            // there are no hanging nodes
            if (N_Hanging == 0)
              continue;

            RHS = rhs[j];
            if(RHS == nullptr)
              continue;

            auto hanging_nodes = fespace->get_sorted_hanging_nodes();
            int hanging_bound = fespace->get_n_active_non_hanging();

            for(auto hn_dof_pair : hanging_nodes)
            {
              int N_ = hn_dof_pair.first->GetN_Nodes();
              auto Coupling = hn_dof_pair.first->GetCoeff();
              auto DOF = hn_dof_pair.first->GetDOF();

              for(int k=0;k<N_;k++)
              {
                int l = DOF[k];
                if(l<hanging_bound)
                {
                  RHS[l] += Coupling[k] * RHS[hn_dof_pair.second];
                }
              }
            }
          } // endfor j

          // ####################################################################
          // write coupling into matrix
          // ####################################################################
          if(!assemble_dirichlet_rows)
          {
            for(j=0;j<n_sqmatrices;j++)
            {
              if(sqmatrices[j] == nullptr)
                continue;
              sqmatrices[j]->correct_hanging_rows();
            } // endfor j
          }

          if(n_rhs)
          {
            delete [] righthand;
            delete [] LocRhs;
          }

          if(N_Parameters)
          {
            delete [] Param[0];
          }

          if(N_AllMatrices)
          {
            delete [] LocMatrices;
            delete [] Matrices[0];
            delete [] Matrices;
          }

        /*
        #ifdef _MPI
          int rank;
          MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

          MPI_Barrier(TDatabase::ParamDB->Comm);

          if(rank==1)
        #endif
          {

          // ####################################################################
          // print the whole matrix -- SECOND
          // ####################################################################
          for(k=0;k<n_sqmatrices;k++)
          {
            cout << "sqmatrix: " << k << endl;
            RowPtr = sqmatrices[k]->get_row_ptr();
            Entries = sqmatrices[k]->GetEntries();
            ColInd = sqmatrices[k]->get_vector_columns();
            N_Rows = sqmatrices[k]->get_n_rows();

            auto fespace = sqmatrices[k]->GetFESpace3D();
            ActiveBound = fespace->GetN_ActiveDegrees();
            for(i=ActiveBound;i<N_Rows;i++)
            {
              end=RowPtr[i+1];
              for(j=RowPtr[i];j<end;j++)
              {
                // cout << j << endl;
                 Output::print("a", k, "(", i+1, ",", ColInd[j]+1, " )=  ", 
                               Entries[j], ";");
              }
            }
            cout << endl;
          } // endfor k

          for(k=0;k<n_matrices;k++)
          {
            cout << endl;
            cout << "matrix: " << k << endl;
            RowPtr = matrices[k]->get_row_ptr();
            Entries = matrices[k]->GetEntries();
            ColInd = matrices[k]->get_vector_columns();
            N_Rows = matrices[k]->get_n_rows();
            for(i=0;i<N_Rows;i++)
            {
              end=RowPtr[i+1];
              for(j=RowPtr[i];j<end;j++)
              {
                // cout << j << endl;
                 Output::print("b", k, "(", i+1, ",", ColInd[j]+1, " )=  "
                               Entries[j], ";");
                      }
            }
            cout << endl;
          } // endfor k

        //   for(k=0;k<n_rhs;k++)
        //   {
        //     cout << "rhs: " << k << endl;
        //     N_Rows = ferhs[k]->get_n_dof();
        //     RHS=rhs[k];
        //     for(i=0;i<N_Rows;i++)
        //       cout << setw(5) << i << setw(20) << RHS[i] << endl;
        //   }
        }
         */
          time_total = GetTime() - time_total;
          //Output::print("LocalToGlobal: ", time_all, " ", time_total, " ",
          //              time_all/time_total);
}
// =======================================================================
//
// Assemble3DSlipBC
//
// some manipulations in matrices and the rhs are necessary
//
// =======================================================================
void Assemble3DSlipBC(int n_fespaces, const TFESpace3D **fespaces,
                      int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                      int n_matrices, TMatrix3D **matrices,
                      int n_rhs, double **rhs, const TFESpace3D **ferhs,
                      BoundCondFunct3D **BoundaryConditions,
                      BoundValueFunct3D **)
{
  double hK;
  const int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,l,l1,m;
  int N_Cells, N_Points, N_;
  int N_Joints;
  const TQuadFormula *qf2;
  double xf, yf, zf;
  double *righthand=nullptr;
  double **Matrices, *aux;
  double ***LocMatrices=nullptr, **LocRhs=nullptr;
  double *AuxArray[MaxN_QuadPoints_3D];
  const int *DOF;
  int ActiveBound;
  double t0;
  BoundCond Cond0;
  BoundCondFunct3D *BoundaryCondition;
  int N_EdgePoints;
  double **JointValues, *JointValue;
  bool *SecondDer;
  double hE;

  const int *TmpFV, *TmpLen;
  int MaxLen;
  double nn, nx, ny, nz, t11, t12, t13, t21, t22, t23, eps=1e-12;
  double x10, y10, z10, x20, y20, z20;
  const double *EdgePoints1, *EdgePoints2;
  double *Entries1=nullptr,*Entries2=nullptr,*Entries3=nullptr;
  double *Entries4=nullptr,*Entries5=nullptr,*Entries6=nullptr,*Entries7=nullptr;
  const int *ColInd1=nullptr, *RowPtr1=nullptr,*ColInd2=nullptr, *RowPtr2=nullptr, *ColInd3=nullptr, *RowPtr3=nullptr;
  const int *ColInd4=nullptr, *RowPtr4=nullptr,*RowPtr5=nullptr,*RowPtr6=nullptr,*RowPtr7=nullptr;
  double *RHS;
  double penetration_penalty, friction_parameter;
  double friction_constant= TDatabase::ParamDB->FRICTION_CONSTANT;
  double friction_power = TDatabase::ParamDB->FRICTION_POWER;
  double penetration_constant = TDatabase::ParamDB->PENETRATION_CONSTANT;
  double penetration_power = TDatabase::ParamDB->PENETRATION_POWER;
  double integral[3], val, det;
  int ii, jj, dof_ii, dof_jj, ll, found;
  int *EdgeDOF, N_EdgeDOF;

// ########################################################################
// store information in local arrays
// ########################################################################

  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions3D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;
  } // endif n_rhs

//  N_Parameters = Parameters->GetN_Parameters();
//  if(N_Parameters)
//  {
//    aux = new double [MaxN_QuadPoints_3D*N_Parameters];
//    for(j=0;j<MaxN_QuadPoints_3D;j++)
//      Param[j] = aux + j*N_Parameters;
//  }

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  if(N_AllMatrices)
  {
    aux = new double
            [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions3D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
  } // endif N_AllMatrices

 SecondDer = new bool[n_fespaces];
 SecondDer[0] = false;
// ########################################################################
// loop over all cells
// ########################################################################
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector< const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  // cout << "Anzahl Zellen: " << N_Cells << endl;
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    //cout << "Zellenr.: " << i << endl;
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      LocalUsedElements[j] = &fespaces[j]->get_fe(i);
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################

    // SecondDer not set !!!
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);

    // ####################################################################
    // manipulate global matrices
    // manipulate global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      // get fe space
      auto fespace = ferhs[j];
      // get current finite element
      auto element = fespace->get_fe(i);
      // get number of basis function of current finite element
      N_ = element.GetN_DOF();

      RHS = rhs[j];
    
      // find bounds in fe space
      ActiveBound = fespace->get_n_active();

      // dof of the rhs nodes connected to this cell
      DOF = fespace->GetGlobalDOF(i);

      // only for faces on the boundary
      BoundaryCondition = BoundaryConditions[j];
      auto nf = element.GetNodalFunctional();
      nf->GetPointsForFace(N_EdgePoints, EdgePoints1, EdgePoints2);

      N_EdgeDOF = element.GetFEDesc()->GetN_JointDOF();

      // find slip type bc
      N_Joints = cell->GetN_Faces();
      for(m=0;m<N_Joints;m++)
      {
        auto joint = cell->GetJoint(m);
        if (joint->GetType() == BoundaryFace)
           //|| (joint->GetType() == IsoBoundaryFace))
        {
            // !!! everything is commented which is not necessary for computing
            // !!! (xf, yf, zf)
            /*if(joint->GetType() == BoundaryFace)
          {
            boundface = (TBoundFace *)joint;
            BoundComp = boundface->GetBoundComp();
            boundface->GetParameters(Param1, Param2);
            }*/
          /* else
          {
            isoboundface = (TIsoBoundFace *)joint;
            BoundComp = isoboundface->GetBoundComp();
            isoboundface->GetParameters(Param1, Param2);
            }*/
          // get id of the boundary component
          //comp=BoundComp->GetID();

          // compute point on the boundary face
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

          t0 = 1.0/TmpLen[m];
          xf = 0; yf = 0; zf = 0;
          double X[4], Y[4], Z[4];
          for (l1=0;l1<TmpLen[m];l1++)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])
                        ->GetCoords(X[l1], Y[l1], Z[l1]);
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }
          /*
          Output::print(xf, " ", yf, " ", zf);

          t0 = 1.0/TmpLen[m];
          
          for (l1=0;l1<TmpLen[m];l1++)
          {
            // get coordinates of the vertices belonging to the 
            // boundary face
            cell->GetVertex(TmpFV[m*MaxLen+l1])
                        ->GetCoords(X[l1], Y[l1], Z[l1]);
            LinComb[l1] = t0;
          }
          // get coordinates (xf, yf, zf) and parameters (t0, t1) of the 
          // barycentre of  the boundary face
          BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                 X, Y, Z, Param1, Param2,
                                 xf, yf, zf, t0, t1);
          Output::print(xf, " ", yf, " ", zf);
          */
          // get boundary condition
          // it is implicitly assumed that only one bc is 
          // assigned to the face
//           BoundaryCondition(xf, yf, zf, Cond0);
          auto boundface              = (TBoundFace*)joint;
          const TBoundComp* BoundComp = boundface->GetBoundComp();
          int comp                    = BoundComp->get_physical_id();
          BoundaryCondition(comp, xf, yf, zf, Cond0);


          switch(Cond0)
          {
            case DIRICHLET:
              break;
              
            case NEUMANN:
              break;
              
            case SLIP:
              ErrThrow("SLIP Condition not implemented here");
              break;
              
            case SLIP_FRICTION_PENETRATION_RESISTANCE:
            {
              // get polynomial degree of finite element in current cell
              const BaseFunctions* bf = element.GetBaseFunct();
              l = bf->GetPolynomialDegree();
              const Shapes* face_types;
              cell->GetShapeDesc()->GetFaceType(face_types);
              // define quadrature formulas
              qf2 = QuadratureFormulaDatabase::qf_from_degree(
                  2*l, face_types[m]);
              // get quad points and weights
              N_Points = qf2->GetN_QuadPoints();
              // get values of test functions in all quadrature points
              // on joint m 
              JointValues = FEDatabase::GetJointDerivatives3D(
                *bf, *qf2, m, MultiIndex3D::D000);

              // compute unit normal vector of the face
              // the face is assumed to be a plane spanned by the first 
              // three vertices
              x10 = X[1] - X[0];
              y10 = Y[1] - Y[0];
              z10 = Z[1] - Z[0];
              x20 = X[2] - X[0];
              y20 = Y[2] - Y[0];
              z20 = Z[2] - Z[0];
              nx = y10*z20-z10*y20;
              ny = z10*x20-z20*x10;
              nz = x10*y20-x20*y10;
              hE= std::sqrt(nx*nx+ny*ny+nz*nz);
              nx /= hE;
              ny /= hE;
              nz /= hE;
              //Output::print(X[0], " ", Y[0], " ", Z[0]); 
              //Output::print(X[1], " ", Y[1], " ", Z[1]); 
              //Output::print(X[2], " ", Y[2], " ", Z[2]);
              // compute two tangential vectors of the face
              if ( (std::abs(nx)>=0.5) || (std::abs(ny)>=0.5))
              {
                nn = std::sqrt(nx*nx+ny*ny);
                t11 = ny/nn;
                t12 = -nx/nn;
                t13 = 0;
                t21 = -t12*nz;
                t22 = t11*nz;
                t23 = t12*nx-t11*ny;
              }
              else
              {
                nn = std::sqrt(ny*ny+nz*nz);
                t11 = 0;
                t12 = -nz/nn;
                t13 = ny/nn;
                t21 = t13*ny-t12*nz;
                t22 = - t13*nx;
                t23 = t12*nx;
              }
              //Output::print("nx ", nx, " ny ", ny, " nz ", nz);
              //Output::print(" t11 ", t11, " t12 ", t12, " t13 ", t13);
              //Output::print(" t21 ", t21, " t22 ", t22, " t23 ", t23);
 
              switch(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY)
              {
                case 0:
                  // standard case, same value on all walls
                  // penalty value for weak imposition of no penetration bc
                  penetration_penalty = penetration_constant*std::pow(hK,penetration_power);              
                  // parameter for friction
                  //BoundaryValue(comp, (t0+t1)/2.0,friction_parameter);
                  friction_parameter = friction_constant * std::pow(hK,friction_power);
                  //Output::print(penetration_penalty, " ", friction_parameter);
                  break;
                case 1:
                  // ChannelStepSlip3D.h
                  penetration_penalty = penetration_constant*std::pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * std::pow(hK,friction_power);
                  // free slip else
                  if ((std::abs(Z[0])<1e-6) && (std::abs(Z[1])<1e-6) && (std::abs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Z[0]-10)<1e-6) && (std::abs(Z[1]-10)<1e-6) && (std::abs(Z[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  break;
                case 2:
                  // ChannelStepSlipFreeSlip.h
                  penetration_penalty = penetration_constant*std::pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * std::pow(hK,friction_power);
                  // free slip else
                  if ((std::abs(Z[0])<1e-6) && (std::abs(Z[1])<1e-6) && (std::abs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Z[0]-10)<1e-6) && (std::abs(Z[1]-10)<1e-6) && (std::abs(Z[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Y[0]-10)<1e-6) && (std::abs(Y[1]-10)<1e-6) && (std::abs(Y[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  break;
                case 3:
                  // CircularCylinder22000.h
                  penetration_penalty = penetration_constant*std::pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * std::pow(hK,friction_power);
                  //Output::print(penetration_penalty, " : ", friction_parameter);
                  // free slip left and right
                  if ((std::abs(Y[0])<1e-6) && (std::abs(Y[1])<1e-6) && (std::abs(Y[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Y[0]-1.4)<1e-6) && (std::abs(Y[1]-1.4)<1e-6) && (std::abs(Y[2]-1.4)<1e-6))
                    friction_parameter = 0.0;
                  // free slip top and bottom
                  if ((std::abs(Z[0])<1e-6) && (std::abs(Z[1])<1e-6) && (std::abs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Z[0]-0.4)<1e-6) && (std::abs(Z[1]-0.4)<1e-6) && (std::abs(Z[2]-0.4)<1e-6))
                    friction_parameter = 0.0;
                  //if ( friction_parameter != 0.0)
                  // Output::print(penetration_penalty, " : ", friction_parameter, 
                  //               " ", X[0], " ", Y[0], " ", Z[0]);
                  break;
                case 4:
                  // WallMountedCube.h
                  penetration_penalty = penetration_constant*std::pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * std::pow(hK,friction_power);
                  //Output::print(penetration_penalty, " : ", friction_parameter);
                  // free slip left and right
                  if ((std::abs(Y[0])<1e-6) && (std::abs(Y[1])<1e-6) && (std::abs(Y[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((std::abs(Y[0]-0.7)<1e-6) && (std::abs(Y[1]-0.7)<1e-6) && (std::abs(Y[2]-0.7)<1e-6))
                    friction_parameter = 0.0;
                  //if ( friction_parameter != 0.0)
                  //Output::print(penetration_penalty, " : ", friction_parameter,
                  //              " ", X[0], " ", Y[0], " ", Z[0]);
                  break;
                case 5:
                  // windtunnel_fine.h
                  // upper wall, symmetrie wall -> free slip, no penetration
                  if ((std::abs(Z[0])<1e-6) && (std::abs(Z[1])<1e-6) && (std::abs(Z[2])<1e-6))
                  {
                    penetration_penalty = penetration_constant*std::pow(hK,penetration_power);
                  }
                  else
                    penetration_penalty = 0.0;
                  friction_parameter = 0.0;
                  break;
                case 6:
                  // windtunnel_m2.h
                  // free slip, no penetration
                  penetration_penalty = 1e12;
                  friction_parameter = 0.0;
                  break;
                default:
                  ErrThrow("INTERNAL_SLIP_WITH_FRICTION_IDENTITY not implemented !!!");
                  break;
                   
              }  

              EdgeDOF = element.GetFEDesc()->GetJointDOF(m);

              // compute additional matrix entries
              // for all velo dof in the mesh cell
              // ii - test function
              for (ii=0;ii<N_;ii++)
              {
                // look for 'ii'-th row in all matrices
                dof_ii = DOF[ii];
                // Dirichlet node
                if (dof_ii>=ActiveBound)
                  continue;
                
                // !!!!!! assumed that A_11 - A_33 are in sqmatrices[0] - [8]
                // ordered as in main program 
                // first velocity component -> matrices A_11, A_12, A13 (and M_11)
                if (j==0)
                {
                  // A11
                  Entries1 = sqmatrices[0]->GetEntries();
                  RowPtr1 = sqmatrices[0]->get_row_ptr();
                  ColInd1 = sqmatrices[0]->get_vector_columns();
                  
                  if (n_sqmatrices>3)
                  {
                    // A12
                    Entries2 = sqmatrices[3]->GetEntries();
                    RowPtr2 = sqmatrices[3]->get_row_ptr();
                    ColInd2 = sqmatrices[3]->get_vector_columns();
                     // A13
                    Entries3 = sqmatrices[4]->GetEntries();
                    RowPtr3 = sqmatrices[4]->get_row_ptr();
                    ColInd3 = sqmatrices[4]->get_vector_columns();
                  }
                  
                  // time dependent problem and NSTYPE 4
                  // entries 4 = M11, entries5 = M12, entries6=M13
                  if (n_sqmatrices == 18)
                  {
                    Entries4 = sqmatrices[9]->GetEntries();
                    RowPtr4  = sqmatrices[9]->get_row_ptr();
                    ColInd4  = sqmatrices[9]->get_vector_columns();
                    Entries5 = sqmatrices[12]->GetEntries();
                    RowPtr5  = sqmatrices[12]->get_row_ptr();
                    Entries6 = sqmatrices[13]->GetEntries();
                    RowPtr6  = sqmatrices[13]->get_row_ptr();
                  }

                  if (n_matrices == 3)
                  {
                    Entries7 = matrices[0]->GetEntries();
                    RowPtr7  = matrices[0]->get_row_ptr();
                  }
                }
                // second velocity component -> matrices A_21, A_22, A23
                if (j==1)
                {
                  if (n_sqmatrices>3)
                  {
                    // A21 
                    Entries1 = sqmatrices[5]->GetEntries();
                    RowPtr1 = sqmatrices[5]->get_row_ptr();
                    ColInd1 = sqmatrices[5]->get_vector_columns();
                    // A23 
                    Entries3 = sqmatrices[6]->GetEntries();
                    RowPtr3 = sqmatrices[6]->get_row_ptr();
                    ColInd3 = sqmatrices[6]->get_vector_columns();
                  }
                  // A22
                  Entries2 = sqmatrices[1]->GetEntries();
                  RowPtr2 = sqmatrices[1]->get_row_ptr();
                  ColInd2 = sqmatrices[1]->get_vector_columns();
                  
                  // time dependent problem and NSTYPE 4
                  // entries 4 = M22, entries5 = M21, entries6=M23
                  if (n_sqmatrices == 18)
                  {
                    Entries4 = sqmatrices[10]->GetEntries();
                    RowPtr4  = sqmatrices[10]->get_row_ptr();
                    ColInd4  = sqmatrices[10]->get_vector_columns();
                    Entries5 = sqmatrices[14]->GetEntries();
                    RowPtr5  = sqmatrices[14]->get_row_ptr();
                    Entries6 = sqmatrices[15]->GetEntries();
                    RowPtr6  = sqmatrices[15]->get_row_ptr();
                  }

                  if (n_matrices == 3)
                  {
                    Entries7 = matrices[1]->GetEntries();
                    RowPtr7  = matrices[1]->get_row_ptr();
                  }
                }
                // third velocity component -> matrices A_31, A_32, A33
                if (j==2)
                {
                  if (n_sqmatrices>3)
                  {
                    // A31 
                    Entries1 = sqmatrices[7]->GetEntries();
                    RowPtr1 = sqmatrices[7]->get_row_ptr();
                    ColInd1 = sqmatrices[7]->get_vector_columns();
                    // A32 
                    Entries2 = sqmatrices[8]->GetEntries();
                    RowPtr2 = sqmatrices[8]->get_row_ptr();
                    ColInd2 = sqmatrices[8]->get_vector_columns();
                  }
                  // A33
                  Entries3 = sqmatrices[2]->GetEntries();
                  RowPtr3 = sqmatrices[2]->get_row_ptr();
                  ColInd3 = sqmatrices[2]->get_vector_columns();
                  
                  // time dependent problem and NSTYPE 4
                  // entries 4 = M33, entries5 = M31, entries6=M32
                  if (n_sqmatrices == 18)
                  {
                    Entries4 = sqmatrices[11]->GetEntries();
                    RowPtr4  = sqmatrices[11]->get_row_ptr();
                    ColInd4  = sqmatrices[11]->get_vector_columns();
                    Entries5 = sqmatrices[16]->GetEntries();
                    RowPtr5  = sqmatrices[16]->get_row_ptr();
                    Entries6 = sqmatrices[17]->GetEntries();
                    RowPtr6  = sqmatrices[17]->get_row_ptr();
                  }

                  if (n_matrices == 3)
                  {
                    Entries7 = matrices[2]->GetEntries();
                    RowPtr7  = matrices[2]->get_row_ptr();
                  }
                }
                //Output::print("ii ", dof_ii);
                
                // for all dof in the mesh cell
                // jj - ansatz function
                for (jj=0;jj<N_;jj++)
                {
                  dof_jj = DOF[jj];
                  // Output::print("jj ", dof_jj);
                  // initialize the boundary integrals
                  for (l=0;l<3;l++)
                    integral[l] = 0;
                  
                  // compute boundary integrals
                  // first component of the velocity
                  if (j==0)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = qf2->get_weight(l)*hE;
                      // (A_11)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nx;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t11;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t21;
                      integral[0] += val*det;
                      // (A_12)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nx;
                      val+= friction_parameter*JointValue[jj]*t11*JointValue[ii]*t12;
                      val+= friction_parameter*JointValue[jj]*t21*JointValue[ii]*t22;
                      integral[1] += val*det;
                      // (A_13)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*nx;
                      val+= friction_parameter*JointValue[jj]*t11*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t21*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // surface not parallel to plane Oyz or penetration
                    if ( TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM || //if weak..
                        (!TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM && //..or strong with conditions
                            ((std::abs(ny)>eps) || (std::abs(nz)>eps) || (penetration_penalty<1e3))) )
                    {
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                    {
                      if (ColInd1[ll] == dof_jj)
                      {
                        //Output::print("a11 ", integral[0]);
                        Entries1[ll] += integral[0];
                        //Output::print("a11 ", integral[0], " ", Entries1[ll]);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      ErrThrow("ERROR A_11 ");
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                      {
                        if (ColInd2[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a12 ", integral[1]);
                          Entries2[ll] += integral[1];
                          //Output::print("a12 ", integral[1], " ", Entries2[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_12 ");
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                      {
                        if (ColInd3[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a13 ", integral[2]);
                          Entries3[ll] += integral[2];
                          //Output::print("a13 ", integral[2], " ", Entries3[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_13 ");
                      }
                    }
                    }
                    else // surface parallel to plane Oyz, i.e. normal to x axis and no penetration
                    {
                      found = 0;
                      for (ll=0; ll<N_EdgeDOF;ll++)
                      {
                        if (dof_ii==DOF[EdgeDOF[ll]])
                        {
                          found = 1;
                          break;
                        }
                      }
                      if (!found)
                        continue;
                      // Output::print("y");
                      // update first matrix, set diagonal entry to 1
                      // all other entries to zero
                      for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                      {
                        if (ColInd1[ll] == dof_ii)
                          Entries1[ll] = 1;
                        else
                          Entries1[ll] = 0;
                      }

                      // update second and third matrix, set all entries to zero
                      if (n_sqmatrices>3)
                      {
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                          Entries2[ll] = 0;

                        // update third matrix
                        for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          Entries3[ll] = 0;
                      }

                      if (n_sqmatrices==18)
                      { // M_11, set off diagonal to zero
                        for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                        {
                          if (ColInd4[ll] != dof_ii)
                            Entries4[ll] = 0;
                          else
                            Entries4[ll] = 1;
                        }
                        // M_12, set row to zero
                        for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                          Entries5[ll] = 0;
                        // M_13, set row to zero
                        for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                          Entries6[ll] = 0;
                      }

                      if (n_matrices==3)
                        for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                          Entries7[ll] = 0;

                      // set rhs to zero
                      RHS[dof_ii] = 0;
                    }
                  } // end first component (j==0)
                  
                  // compute boundary integrals
                  // second component of the velocity
                  if (j==1)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = qf2->get_weight(l)*hE;
                      // (A_21)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*ny;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t12;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t22;
                      integral[0] += val*det;
                      // (A_22)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*ny;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t12;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t22;
                      integral[1] += val*det;
                      // (A_23)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*ny;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // surface not parallel to plane Oxz or penetration
                    if ( TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM || //if weak..
                        (!TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM && //..or strong with conditions
                         ((std::abs(nx)>eps) || (std::abs(nz)>eps) || (penetration_penalty<1e3))) )
                    {
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                    {
                      if (ColInd2[ll] == dof_jj)
                      {
                          //Output::print("a22 ", integral[1]);
                        Entries2[ll] += integral[1];
                        // if (std::abs(ny)!>1e-6)
                        //Output::print("a22 ", integral[1], " ", Entries2[ll]);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      ErrThrow("ERROR A_22 ");
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                      {
                        if (ColInd1[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a21 ", integral[0]);
                          Entries1[ll] += integral[0];
                          //Output::print("a21 ", Entries1[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_21 ");
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                      {
                        if (ColInd3[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a23 ", integral[2]);
                          Entries3[ll] += integral[2];
                          //Output::print("a23 ", Entries3[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_23 ");
                      }
                    }
                    }
                    else // surface parallel to plane Oxz, i.e. normal to y axis and no penetration
                    {
                      found = 0;
                      for (ll=0; ll<N_EdgeDOF;ll++)
                      {
                        if (dof_ii==DOF[EdgeDOF[ll]])
                        {
                          found = 1;
                          break;
                        }
                      }
                      if (!found)
                        continue;
                      //Output::print("y ", yf, " ", n_sqmatrices);
                      // update first matrix, set all entries to zero
                      if (n_sqmatrices>3)
                      {
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                          Entries1[ll] = 0;

                        // update second matrix, set diagonal entry to 1
                        // all other entries to zero
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                        {
                          if (ColInd2[ll] == dof_ii)
                            Entries2[ll] = 1;
                          else
                            Entries2[ll] = 0;
                        }

                        // update third matrix, set all entries to zero
                        for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          Entries3[ll] = 0;
                      }

                      if (n_sqmatrices==18)
                      { // M_22, set off diagonal to zero
                        for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                        {
                          if (ColInd4[ll] != dof_ii)
                            Entries4[ll] = 0;
                          else
                            Entries4[ll] = 1;
                        }
                        // M_21, set row to zero
                        for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                          Entries5[ll] = 0;
                        // M_23, set row to zero
                        for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                          Entries6[ll] = 0;
                      }

                      if (n_matrices==3)
                        for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                          Entries7[ll] = 0;

                      // set rhs to zero
                      RHS[dof_ii] = 0;
                    }
                  } // end second component (j==1)


                  // compute boundary integrals
                  // third component of the velocity
                  if (j==2)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = qf2->get_weight(l)*hE;
                      // (A_31)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nz;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t13;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t23;
                      integral[0] += val*det;
                      // (A_23)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nz;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t23;
                      integral[1] += val*det;
                      // (A_33)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*nz;
                      val+= friction_parameter*JointValue[jj]*t13*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t23*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // surface not parallel to plane Oxy or penetration
                    if ( TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM || //if weak..
                       (!TDatabase::ParamDB->INTERNAL_SLIP_WEAK_FORM && //..or strong with conditions
                           ((std::abs(nx)>eps) || (std::abs(ny)>eps) || (penetration_penalty<1e3))) )
                    {
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                    {
                      if (ColInd3[ll] == dof_jj)
                      {
                        Entries3[ll] += integral[2];
                        //Output::print("a33 ", Entries3[ll]);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      ErrThrow("ERROR A_33 ");
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                      {
                        if (ColInd1[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a31 ", integral[0]);
                          Entries1[ll] += integral[0];
                          //Output::print("a31 ", Entries1[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_31 ");
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                      {
                        if (ColInd2[ll] == dof_jj)
                        {
                          found = 1;
                          //Output::print("a32 ", integral[1]);
                          Entries2[ll] += integral[1];
                          //Output::print("a32 " Entries2[ll]);
                          break;
                        }
                      }
                      if (!found)
                      {
                        ErrThrow("ERROR A_32 ");
                      }
                    }
                    }
                    else // surface parallel to plane Oxy, i.e. normal to z axis and no penetration
                    {
                      found = 0;
                      for (ll=0; ll<N_EdgeDOF;ll++)
                      {
                        if (dof_ii==DOF[EdgeDOF[ll]])
                        {
                          found = 1;
                         break;
                        }
                      }
                      if (!found)
                        continue;

                      if (n_sqmatrices>3)
                      {
                        // update first matrix, set all entries to zero
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                          Entries1[ll] = 0;

                        // update second matrix, set all entries to zero
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                          Entries2[ll] = 0;
                      }

                      // update third matrix, set diagonal entry to 1
                      // all other entries to zero
                      for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                      {
                        if (ColInd3[ll] == dof_ii)
                          Entries3[ll] = 1;
                        else
                          Entries3[ll] = 0;
                      }

                      if (n_sqmatrices==18)
                      { // M_33, set off diagonal to zero
                        for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                        {
                          if (ColInd4[ll] != dof_ii)
                            Entries4[ll] = 0;
                          else
                            Entries4[ll] = 1;
                        }
                        // M_31, set row to zero
                        for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                          Entries5[ll] = 0;

                        // M_32, set row to zero
                        for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                          Entries6[ll] = 0;
                      }

                      if (n_matrices==3)
                      {
                        for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                          Entries7[ll] = 0;
                      }

                      // set rhs to zero
                      RHS[dof_ii] = 0;
                    }
                  } // end third component (j==2)
                  
                }// endfor ansatz functions (jj)
              }// enfor test functions (ii)
              break;
            }
            default:
           
           break;              
              
         } // end switch (Cond0)
        } // endif (boundary joint)
      }  // endfor m (N_Joints)
    } // endfor j (n_rhs)
  }  // endfor i (N_Cells)

  if(n_rhs)
  {
    delete [] righthand;
    delete [] LocRhs;
  }

//  if(N_Parameters)
//  {
//    delete [] Param[0];
//  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];
  delete [] SecondDer;

/*
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  for(k=0;k<n_sqmatrices;k++)
  {
    cout << endl;
    cout << "sqmatrix: " << k << endl;
    RowPtr = sqmatrices[k]->get_row_ptr();
    Entries = sqmatrices[k]->GetEntries();
    ColInd = sqmatrices[k]->get_vector_columns();
    N_Rows = sqmatrices[k]->get_n_rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
        cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
        cout << setw(10) << Entries[j] << endl;
      }
    }
    cout << endl;
  } // endfor k
  
  for(k=0;k<n_matrices;k++)
  {
    cout << endl;
    cout << "matrix: " << k << endl;
    RowPtr = matrices[k]->get_row_ptr();
    Entries = matrices[k]->GetEntries();
    ColInd = matrices[k]->get_vector_columns();
    N_Rows = matrices[k]->get_n_rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
        cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
        cout << setw(10) << Entries[j] << endl;
      }
    }
    cout << endl;
  } // endfor k

  for(k=0;k<n_rhs;k++)
  {
    cout << "rhs: " << k << endl;
    N_Rows = ferhs[k]->get_n_dof();
    RHS=rhs[k];
    for(i=0;i<N_Rows;i++)
      cout << setw(5) << i << setw(20) << RHS[i] << endl;
  }
*/

} // end of Assemble

/*************************************+******************************/
// 
// Modification of Matrices for slip with friction bc for 
// better condition
//
/*******************************************************************/
void ModifyMatrixSlipBC(TSquareMatrix3D **sqmatrices, TMatrix3D **matrices,
    int N_U, double *rhs)
{
    int i, j, k, j0, j1, index, row_off, col_off;
    const int *RowPtr, *ColInd;
    double *setzero, bound, *Entries, maximal = -1;; 

    return;

    bound = TDatabase::ParamDB->PENETRATION_CONSTANT*1e-4;

    if (std::abs(TDatabase::ParamDB->PENETRATION_POWER+2) >1e6)
      Output::print("ModifyMatrixSlipBC : recommended to set PENETRATION_POWER to -2 !!!");

    setzero = new double[3*N_U];
    memset(setzero, 0, sizeof(double)*3*N_U);
    
    // find largest values of the matrices A_1k
    for (k=0;k<3;k++)
    {
        ColInd = sqmatrices[3*k]->get_vector_columns();
        RowPtr = sqmatrices[3*k]->get_row_ptr();
        Entries = sqmatrices[3*k]->GetEntries();
        for (i=0;i<N_U;i++)
        {
            // i-th row of sqmatrix
            j0 = RowPtr[i];
            j1 = RowPtr[i+1];
            for(j=j0;j<j1;j++)
            {
                // column
                index = ColInd[j];
                if (std::abs(Entries[j])>bound)
                    setzero[index] = 1;
                if (std::abs(Entries[j])>maximal)
                    maximal = std::abs(Entries[j]);         
            }
        }
    }

    // find largest values of the matrices A_2k
    for (k=0;k<3;k++)
    {
        ColInd = sqmatrices[3*k+1]->get_vector_columns();
        RowPtr = sqmatrices[3*k+1]->get_row_ptr();
        Entries = sqmatrices[3*k+1]->GetEntries();
        for (i=0;i<N_U;i++)
        {
            // i-th row of sqmatrix
            j0 = RowPtr[i];
            j1 = RowPtr[i+1];
            for(j=j0;j<j1;j++)
            {
                // column
                index = ColInd[j];
                if (std::abs(Entries[j])>bound)
                    setzero[N_U+index] = 1;
                if (std::abs(Entries[j])>maximal)
                    maximal = std::abs(Entries[j]);         
            }
        }
    }
    // find largest values of the matrices A_3k
    for (k=0;k<3;k++)
    {
        ColInd = sqmatrices[3*k+2]->get_vector_columns();
        RowPtr = sqmatrices[3*k+2]->get_row_ptr();
        Entries = sqmatrices[3*k+2]->GetEntries();
        for (i=0;i<N_U;i++)
        {
            // i-th row of sqmatrix
            j0 = RowPtr[i];
            j1 = RowPtr[i+1];
            for(j=j0;j<j1;j++)
            {
                // column
                index = ColInd[j];
                if (std::abs(Entries[j])>bound)
                    setzero[2*N_U+index] = 1;
                if (std::abs(Entries[j])>maximal)
                    maximal = std::abs(Entries[j]);         
            }
        }
    }
    Output::print("maximal matrix entry ", maximal);

    // set homogeneous Dirichlet bc for the marked nodes
    // loop over the matrix blocks A_ij
    for (k=0;k<9;k++)
    {
        ColInd = sqmatrices[k]->get_vector_columns();
        RowPtr = sqmatrices[k]->get_row_ptr();
        Entries = sqmatrices[k]->GetEntries();
        switch(k)
        {
            case 0:
                row_off = 0;
                col_off = 0;
                break;
            case 1:
                row_off = 0;
                col_off = N_U;
                break;
            case 2:
                row_off = 0;
                col_off = 2*N_U;
                break;
            case 3:
                row_off = N_U;
                col_off = 0;
                break;
            case 4:
                row_off = N_U;
                col_off = N_U;
                break;
            case 5:
                row_off = N_U;
                col_off = 2*N_U;
                break;
            case 6:
                row_off = 2*N_U;
                col_off = 0;
                break;
            case 7:
                row_off = 2*N_U;
                col_off = N_U;
                break;
            case 8:
                row_off = 2*N_U;
                col_off = 2*N_U;
                break;
        }
        // loop over the dof
        for (i=0;i<N_U;i++)
        {
            // i-th row of sqmatrix
            j0 = RowPtr[i];
            j1 = RowPtr[i+1];
            // find diagonal
            for(j=j0;j<j1;j++)
            {
                index = ColInd[j];
                if (setzero[i+row_off])
                {
                    if ((index==i)&&(row_off==col_off))
                        Entries[j] = 1.0;
                    else 
                        Entries[j] = 0.0;
                    continue;
                }
                if (setzero[index+col_off])
                    Entries[j] = 0.0;
            }
        }
    }
   
    // loop over the matrix blocks B_k
    for (k=0;k<3;k++)
    {
        ColInd = matrices[k]->get_vector_columns();
        RowPtr = matrices[k]->get_row_ptr();
        Entries = matrices[k]->GetEntries();
        switch(k)
        {
            case 0:
                row_off = 0;
                break;
            case 1:
                row_off = N_U;
                break;
            case 2:
                row_off = 2*N_U;
                break;
        }
        for (i=0;i<N_U;i++)
        {
            // i-th row of sqmatrix
            j0 = RowPtr[i];
            j1 = RowPtr[i+1];
            for(j=j0;j<j1;j++)
            {
                if (setzero[i+row_off])
                {
                    Entries[j] = 0.0;
                }
            }
        }
    }
    // rhs
    for (i=0;i<3*N_U;i++)
    {
        if (setzero[i])
            rhs[i] = 0;
    }

    delete  [] setzero;
}

/*
  Assemble3D_mixed:
    Assemble for vector finite elements (Raviart-Thomas)
    Need the global orientation of normal at each inner edge/face

    implementation: Alfonso (07.09.2010)
*/
void Assemble3D_mixed(int n_fespaces, const TFESpace3D **fespaces,
int n_sqmatrices, TSquareMatrix3D **sqmatrices,
int n_matrices, TMatrix3D **matrices,
int n_rhs, double **rhs, const TFESpace3D **ferhs,
LocalAssembling3D& la,
BoundCondFunct3D **BoundaryConditions,
BoundValueFunct3D **BoundaryValues)
{
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int N_Points, N_;
  int N_Test, N_Ansatz, N_Joints;
  const TQuadFormula *qf2;
  TJoint *joint;
  const double *xi, *eta, *zeta;
  double xf, yf, zf;
  double *local_rhs;
  double *righthand =nullptr;
  double **Matrices =nullptr, *aux;
  double **Matrix;
  double ***LocMatrices =nullptr, **LocRhs =nullptr;
  std::vector<const BaseFunctions*> LocBF(n_fespaces);
  const int *DOF;
  int ActiveBound, begin, end, middle;
  const int *TestDOF, *AnsatzDOF;
  double *Entries;
  const int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double t0, t1, t2;
  BoundCond Cond0;
  BoundCondFunct3D *BoundaryCondition;
  BoundValueFunct3D *BoundaryValue;
  //TOutput3D *Output;
  ReferenceTransformation_type reftrans;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_BaseFunctions3D];
  int *EdgeDOF, N_EdgeDOF;
  double **JointValues, *JointValue;

  // static bool *SecondDer = nullptr;
  bool *SecondDer;
  double LinComb[4];

  const int *TmpFV, *TmpLen;
  int MaxLen;
  double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
  double nx, ny, nz;
  
  double time1, time2;
  
  bool OuterBoundary;

  double time_total = GetTime();
  double time_all = 0;

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions3D];
    for(int i = 0; i < n_rhs; ++i)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;
  }                                               // endif n_rhs
  
  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
    for(int j = 0; j < N_AllMatrices*MaxN_BaseFunctions3D; ++j)
      Matrices[j] = aux+j*MaxN_BaseFunctions3D;

    LocMatrices = new double** [N_AllMatrices];
    for(int i = 0; i < N_AllMatrices; ++i)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
  }                                               // endif N_AllMatrices
  SecondDer = la.GetNeeds2ndDerivatives();
  

  // all spaces use same Coll
  auto Coll = fespaces[0]->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  for(int i = 0; i < N_Cells; ++i)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector< const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  

  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(int i = 0; i < N_Cells; ++i)
  {
    auto cell = Coll->GetCell(i);
    cell->SetNormalOrientation();
    
    double hK = 0.;
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
      // cases 4 and 5 are specially treated for special problems
      case 0:                                     // diameter
      case 4:                                     // diameter
      case 5:                                     // piecewise constant array
        hK = cell->GetDiameter();
        break;
      // case 3 is specially treated for special problem
      case 1:                                     // with reference map
      case 3:                                     // measure
        hK = cell->GetMeasure();
        hK = std::pow(hK,1.0/3.0);
        break;
      case 2:                                     // shortest edge
        hK = cell->GetShortestEdge();
        break;
      default:                                   // diameter
        hK = cell->GetDiameter();
        Output::print("CELL_MEASURE ", TDatabase::ParamDB->CELL_MEASURE,
                      " not available, set CELL_MEASURE: 0 !!!");
        TDatabase::ParamDB->CELL_MEASURE = 0;
        break;
    }
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(int j = 0; j < n_fespaces; ++j)
    {
      auto& element = fespaces[j]->get_fe(i);
      LocalUsedElements[j] = &element;
      LocBF[j] = element.GetBaseFunct();
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    //Output::print("CELL ", i);
    reftrans = FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer,
                                   qf_ref, qf_orig);
    
    // assemble a few matrices and right-hand sides at once
    la.GetLocalForms(qf_orig, LocBF, cell, i, N_AllMatrices, n_rhs, LocMatrices,
                     LocRhs);

    time1 = GetTime();
    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################

    for(int j = 0; j < n_sqmatrices; ++j)
    {
      // find space for this bilinear form
      auto fespace = sqmatrices[j]->GetFESpace3D();
      auto element = fespace->get_fe(i);
      N_ = element.GetN_DOF();
      Matrix = Matrices+j*MaxN_BaseFunctions3D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr = sqmatrices[j]->get_row_ptr();
      ColInd = sqmatrices[j]->get_vector_columns();

      ActiveBound = fespace->get_n_active();
      DOF = fespace->GetGlobalDOF(i);
      // add local matrix to global
      for(int m = 0; m < N_; ++m)
      {
        int l=DOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          // for all dof
          for(int k = 0; k < N_; ++k)
          {

            begin = RowPtr[l];
            end = RowPtr[l+1]-1;
            int l1 = DOF[k];
            if(ColInd[begin] == l1)
            {
              Entries[begin] += MatrixRow[k];
            }
            else
            {
              begin++;
              middle = (begin+end)/2;
              while(middle>begin)
              {
                if(l1 == ColInd[middle]) break;
                if(l1 < ColInd[middle]) end = (begin+end)/2;
                if(l1 > ColInd[middle]) begin = (begin+end)/2;
                middle = (begin+end)/2;
              }
              if(ColInd[middle] == l1)
              {
                Entries[middle] += MatrixRow[k];
              }
              else
              {
                if(ColInd[middle+1] == l1)
                {
                  Entries[middle+1] += MatrixRow[k];
                }
              }
            }
          }                                       // endfor k
        }                                         // endif l
        else
        {
          // Dirichlet node
          int n=RowPtr[l];
          if(ColInd[n]==l)
          {
            Entries[n]=1.0;
          }
        }
      }                                           // endfor m
    }                                             // endfor j
    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(int j = 0; j < n_matrices; ++j)
    {
      auto TestElement = matrices[j]->GetTestSpace3D()->get_fe(i);
      auto AnsatzElement = matrices[j]->GetAnsatzSpace3D()->get_fe(i);

      N_Test = TestElement.GetN_DOF();
      N_Ansatz = AnsatzElement.GetN_DOF();

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions3D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->get_row_ptr();
      ColInd = matrices[j]->get_vector_columns();

      TestDOF = matrices[j]->GetTestSpace3D()->GetGlobalDOF(i);
      AnsatzDOF = matrices[j]->GetAnsatzSpace3D()->GetGlobalDOF(i);

      int ActiveBound = matrices[j]->GetTestSpace3D()->get_n_active();
      
      // add local matrix to global
      for(int m = 0; m < N_Test; ++m)
      {
        int l=TestDOF[m];
        if(l >= ActiveBound)
          continue;
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        for(int k = 0; k < N_Ansatz; ++k)
        {
          int l1 = AnsatzDOF[k];
          begin = RowPtr[l];
          end = RowPtr[l+1]-1;
          middle = (begin+end)/2;
          while(middle>begin)
          {
            if(l1 == ColInd[middle]) break;
            if(l1 < ColInd[middle]) end = (begin+end)/2;
            if(l1 > ColInd[middle]) begin = (begin+end)/2;
            middle = (begin+end)/2;
          }
          if(ColInd[middle] == l1)
          {
            Entries[middle] += MatrixRow[k];
          }
          else
          {
            if(ColInd[middle+1] == l1)
            {
              Entries[middle+1] += MatrixRow[k];
            }
          }
        }                                         // endfor k
      }                                           // endfor m
    }                                             // endfor j
    time2 = GetTime();
    time_all += time2-time1;
    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(int j = 0; j < n_rhs; ++j)
    {
      //Output::print("rhs ", j);
      auto fespace = ferhs[j];
      ActiveBound = fespace->get_n_active();
      auto element = fespace->get_fe(i);

      N_ = element.GetN_DOF();

      local_rhs = righthand+j*MaxN_BaseFunctions3D;
      RHS = rhs[j];
      // find space for this linear form

      ActiveBound = fespace->get_n_active();
      DOF = fespace->GetGlobalDOF(i);

      // add local right-hand side to the global one
      for(int m = 0; m < N_; ++m)
      {
        int l=DOF[m];
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
        }                                         // endif l
      }                                           // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      auto nf = element.GetNodalFunctional();
      auto FEDesc_Obj = element.GetFEDesc();
      
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Faces();
      for(int m = 0; m < N_Joints; ++m)
      {
        joint = cell->GetJoint(m);
        OuterBoundary = false;

        if(joint->GetType() == BoundaryFace ||
          joint->GetType() == IsoBoundFace)
          OuterBoundary = true;

       
        if(OuterBoundary)
        {
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

          t0 = 1.0/TmpLen[m];

          // coordines of center are computed, this is where boundary 
          // conditions are taken
          double X[4], Y[4], Z[4];
          xf = 0; yf = 0; zf = 0;
          for (int l1 = 0; l1 < TmpLen[m]; ++l1)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }

          // the face gets the b.c. which is valid at its center
//           BoundaryCondition(xf, yf, zf, Cond0);
          auto boundface              = (TBoundFace*)joint;
          const TBoundComp* BoundComp = boundface->GetBoundComp();
          int comp                    = BoundComp->get_physical_id();
          BoundaryCondition(comp, xf, yf, zf, Cond0);


          switch(Cond0)
          {
            case DIRICHLET:
            {
              nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
              std::vector<double> XYZ(3*N_Points);
              FEDatabase::GetOrigFromRef(reftrans, N_Points, xi, eta, zeta,
                                         XYZ.data(), XYZ.data()+N_Points,
                                         XYZ.data()+2*N_Points);
              
              for(int l1 = 0; l1 < N_Points; ++l1)
              {
                BoundaryValue(comp, XYZ[l1], XYZ[l1+N_Points],
                              XYZ[l1+2*N_Points], PointValues[l1]);
              }

              nf->GetFaceFunctionals(Coll, cell, m, PointValues,
                                     FunctionalValues);
              EdgeDOF = FEDesc_Obj->GetJointDOF(m);
              N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              
              for(int l = 0; l < N_EdgeDOF; ++l)
              {
                RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              }
              break;
            }
            case NEUMANN:
            {
              // Output::print("Neumann condition in Assemble3D");
              const BaseFunctions* bf = element.GetBaseFunct();
              int l = bf->GetPolynomialDegree();
              const Shapes* face_types;
              cell->GetShapeDesc()->GetFaceType(face_types);
              qf2 = QuadratureFormulaDatabase::qf_from_degree(
                  2*l, face_types[m]);
              N_Points = qf2->GetN_QuadPoints();
              // values of base functions in all quadrature points on face
              JointValues = FEDatabase::GetJointDerivatives3D(
                *bf, *qf2, m, MultiIndex3D::D000);
              bf->ChangeBF(Coll, cell, N_Points, JointValues);

              switch(TmpLen[m])
              {
                case 3:
                  xc1 = X[1] - X[0];
                  xc2 = X[2] - X[0];

                  yc1 = Y[1] - Y[0];
                  yc2 = Y[2] - Y[0];

                  zc1 = Z[1] - Z[0];
                  zc2 = Z[2] - Z[0];

                  // normal vector
                  nx = yc1*zc2 - zc1*yc2;
                  ny = zc1*xc2 - xc1*zc2;
                  nz = xc1*yc2 - yc1*xc2;
                  // determinant of reference trafo
                  t2 = std::sqrt(nx*nx + ny*ny + nz*nz);

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    auto p = qf2->get_point(l);
                    t0 = p.x;
                    t1 = p.y;
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = 1-t0-t1;
                    LinComb[1] = t0;
                    LinComb[2] = t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                      +LinComb[2]*X[2];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                      +LinComb[2]*Y[2];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                      +LinComb[2]*Z[2];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
                    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(comp, xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    // cout << t1 << endl;
                    t0 *= qf2->get_weight(l)*t2;
                    for(int k = 0; k < N_; ++k)
                    {
                      int l3 = DOF[k];
                      if(l3 < ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                    }
                  }                               // endfor l
                  break;

                case 4:
                  xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
                  xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
                  xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

                  yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
                  yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
                  yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

                  zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
                  zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
                  zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    auto p = qf2->get_point(l);
                    t0 = 0.5*(p.x+1);
                    t1 = 0.5*(p.y+1);
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = (1-t0)*(1-t1);
                    LinComb[1] = t0*(1-t1);
                    LinComb[2] = t0*t1;
                    LinComb[3] = (1-t0)*t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                      +LinComb[2]*X[2] + LinComb[3]*X[3];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                      +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                      +LinComb[2]*Z[2] + LinComb[3]*Z[3];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
                    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(comp, xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    nx = (yc1+p.y*yc3)*(zc2+p.x*zc3)
                      -(zc1+p.y*zc3)*(yc2+p.x*yc3);
                    ny = (zc1+p.y*zc3)*(xc2+p.x*xc3)
                      -(xc1+p.y*xc3)*(zc2+p.x*zc3);
                    nz = (xc1+p.y*xc3)*(yc2+p.x*yc3)
                      -(yc1+p.y*yc3)*(xc2+p.x*xc3);
                    t1 = nx*nx+ny*ny+nz*nz;
                    // cout << t1 << endl;
                    t0 *= qf2->get_weight(l)*std::sqrt(t1);
                    for(int k = 0; k < N_; ++k)
                    {
                      int l3 = DOF[k];
                      if(l3 < ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                    }
                  }                               // endfor l
                  break;
              }
              break;
            }
         default:
           break;
          }                                       // endswitch Cond0=
        }                                         // endif
      }                                           // endfor m
    }                                             // endfor j
    // cout << "end i:" << i << endl;

  }                                               // endfor i (loop over all cells)
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  
  // ####################################################################
  // write coupling into matrix
  // ####################################################################

  if(n_rhs)
  {
    delete [] righthand;
    delete [] LocRhs;
  }


  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  time_total = GetTime() - time_total;
}
