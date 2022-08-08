// =======================================================================
// @(#)Assemble2D.C        1.16 04/13/00
//
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#include <Assemble2D.h>
#include <Enumerations_fe.h>
#include <Matrix2D.h>
#include <AuxParam2D.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include "FEDatabase.h"
#include "HangingNode.h"
#include <NodalFunctional.h>
#include <SquareMatrix2D.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <Convolution.h>
#include "BaseCell.h"
#include "BoundaryAssembling2D.h"
#include "QuadratureFormulaDatabase.h"

//#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <vector>

// =======================================================================
//
// Assemble2DSlipBC
//
// some manipulations in matrices and the rhs are necessary
//
// =======================================================================

void Assemble2DSlipBC(int n_fespaces, const TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, const TFESpace2D **ferhs,
BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **,
TFEFunction2D *u1, TFEFunction2D *u2)
{
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,l,m, ii,jj,ll;
  int N_Cells, N_;
  int N_Joints;
  std::vector<const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  const TQuadFormula *qf1;
  double *righthand = nullptr;
  double **Matrices = nullptr;
  double *aux = nullptr;
  double ***LocMatrices = nullptr;
  double **LocRhs = nullptr;
  double *AuxArray[MaxN_QuadPoints_2D];
  const int *DOF;
  int ActiveBound;

  double *Entries1 = nullptr;
  double *Entries2 = nullptr;
  double *Entries3 = nullptr;
  double *Entries4 = nullptr;
  double *Entries5 = nullptr;

  const int *ColInd1 = nullptr;
  const int* RowPtr1 = nullptr;
  const int* ColInd2 = nullptr;
  const int* RowPtr2 = nullptr;
  const int* ColInd3 = nullptr;
  const int* RowPtr3 = nullptr;
  const int* RowPtr4 = nullptr;
  const int* RowPtr5 = nullptr;
  double *RHS;
  const TBoundComp *BoundComp;
  double t0, t1, s,integral[2];
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double x0 = 0;
  double x1 = 0;
  double y0 = 0;
  double y1 = 0;
  double hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty;
  double friction_parameter = 0;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double friction_constant= TDatabase::ParamDB->FRICTION_CONSTANT;
  double friction_power = TDatabase::ParamDB->FRICTION_POWER;
  double penetration_constant = TDatabase::ParamDB->PENETRATION_CONSTANT;
  double penetration_power = TDatabase::ParamDB->PENETRATION_POWER;
  int friction_type = TDatabase::ParamDB->FRICTION_TYPE;
  double RE_NR, tangential_velo, U0, denominator;
  bool *SecondDer;
#ifdef __3D__
  double z0, z1;
#endif


  // ########################################################################
  // store information in local arrays
  // ########################################################################

  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_2D*40];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*40;

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  SecondDer = new bool[n_fespaces];
  SecondDer[0] = false;
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  // ########################################################################
  // loop over all cells
  // ########################################################################
  auto Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

    // double hK = cell->GetDiameter();
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
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);

    // ####################################################################
    // manipulate global matrices
    // manipulate global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      auto fespace = ferhs[j];
      auto ele = fespace->get_fe(i);
      N_ = ele.GetN_DOF();

      RHS = rhs[j];

      // find bounds in fe space
      ActiveBound = fespace->get_n_active();

      // dof of the rhs nodes connected to this cell
      DOF = fespace->GetGlobalDOF(i);

      // only for edges on the boundary
      BoundaryCondition = BoundaryConditions[j];      
      // auto nf = ele.GetNodalFunctional2D();
      // nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      auto FEDesc_Obj = ele.GetFEDesc();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

      // find slip type bc
      N_Joints = cell->GetN_Edges();
      for(m=0;m<N_Joints;m++)
      {
        auto joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge)
          {
            auto boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            auto isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            switch(Cond0)
            {
              case DIRICHLET:
                break;

              case NEUMANN:
                break;

              case SLIP:
                ErrThrow("Slip boundary condition not implemented here");
                break;

              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // edge is assumed to be straight line
                // get polynomial degree of fe
                l = ele.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                qf1 = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                N_LinePoints = qf1->GetN_QuadPoints();
                // get values of test functions in all quadrature points
                // on joint m
                JointValues=FEDatabase::GetJointDerivatives2D(
                    *ele.GetBaseFunct(), *qf1, m, MultiIndex2D::D00);
                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                // compute normal vector to this boundary (normalized)
                nx = (y1-y0)/hE;
                ny = (x0-x1)/hE;
                // tangential normal vector to this boundary (normalized)
                tx = (x1-x0)/hE;
                ty = (y1-y0)/hE;
                // Output::print("A ", x0, " ", y0, " B ", x1, " ", y1);
                // double delta = CharacteristicFilterWidth(hK);
#ifdef  __CHANNELSTEPSLIP__
                // upper boundary - free slip
                if (comp==6)
                  friction_constant = 0;
#endif

                // penalty value for weak imposition of no penetration bc
                penetration_penalty = penetration_constant*std::pow(hE,penetration_power);;

                // parameter for friction
                //BoundaryValue(comp, (t0+t1)/2.0,friction_parameter);
                switch(friction_type)
                {
                  case 1:
                  default:
                    // linear friction
                    friction_parameter = friction_constant * std::pow(hE,friction_power);
                    //Output::print(" friction parameter ", friction_parameter);
                    break;
                  case 2:
                    // nonlinear type 1
                    // beta = std::sqrt(Re)(w\cdot tau)/(U_0/2-(w\cdot tau)
                    // centre of the boundary edge
                    x = (x1+x0)/2.0;
                    y = (y1+y0)/2.0;
                    // compute velocity in (x,y);
                    u1->FindGradientLocal(cell,i,x,y,u1_values);
                    u2->FindGradientLocal(cell,i,x,y,u2_values);
                    // compute tangential velocity
                    tangential_velo = u1_values[0]*tx + u2_values[0]*ty;
                    // get Reynolds number
                    RE_NR = TDatabase::ParamDB->RE_NR;
                    U0 = TDatabase::ParamDB->FRICTION_U0;
                    denominator = U0/2-tangential_velo;
                    if (std::abs(denominator)<1e-8)
                    {
                      ErrThrow("nonlinear slip bc type 1, denominator zero !!!");
                    }
                    friction_parameter = std::sqrt(RE_NR) * tangential_velo/denominator;
                    friction_parameter = std::abs(friction_parameter);
                    Output::print("x ", x, " y ", y);
                    Output::print(" friction parameter ", friction_parameter);
                    break;
                }

                hE = hE/2;

                EdgeDOF = FEDesc_Obj->GetJointDOF(m);

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

                  // !!!!!! assumed that A_11 - A_22 are in sqmatrices[0] - [3]
                  // COMMENT UPDATE (2018-08-10) A_11 - A_22 are in sqmatrices[0] - [4]
                  // first velocity component -> matrices A_11 and A_12 (and M_11)
                  if (j==0)
                  {
                    Entries1 = sqmatrices[0]->GetEntries();
                    RowPtr1 = sqmatrices[0]->get_row_ptr();
                    ColInd1 = sqmatrices[0]->get_vector_columns();

                    if (n_sqmatrices>2)  //nstype 3 and 4
                    {
                      Entries2 = sqmatrices[2]->GetEntries();
                      RowPtr2 = sqmatrices[2]->get_row_ptr();
                      ColInd2 = sqmatrices[2]->get_vector_columns();
                    }

                    // time dependent problem and NSTYPE 4
                    // entries 3 = M11, entries4 = M12
                    if (n_sqmatrices==8)
                    {
                      Entries3 = sqmatrices[4]->GetEntries();
                      RowPtr3 = sqmatrices[4]->get_row_ptr();
                      ColInd3 = sqmatrices[4]->get_vector_columns();
                      Entries4 = sqmatrices[6]->GetEntries();
                      RowPtr4 = sqmatrices[6]->get_row_ptr();
                    }

                    if (n_matrices==2)
                    {
                      Entries5 = matrices[0]->GetEntries();
                      RowPtr5 = matrices[0]->get_row_ptr();                      
                    }
                  }
                  // second velocity component -> matrices A_21 and A_22
                  if (j==1)
                  {
                    if (n_sqmatrices>2)
                    { // entries1 = A21
                      Entries1 = sqmatrices[3]->GetEntries();
                      RowPtr1 = sqmatrices[3]->get_row_ptr();
                      ColInd1 = sqmatrices[3]->get_vector_columns();
                    }
                    // entries 2 = A22
                    Entries2 = sqmatrices[1]->GetEntries();
                    RowPtr2 = sqmatrices[1]->get_row_ptr();
                    ColInd2 = sqmatrices[1]->get_vector_columns();

                    // time dependent problem and NSTYPE 4
                    // entries 3 = M22, entries4 = M21
                    if (n_sqmatrices==8)
                    {
                      Entries3 = sqmatrices[5]->GetEntries();
                      RowPtr3 = sqmatrices[5]->get_row_ptr();
                      ColInd3 = sqmatrices[5]->get_vector_columns();
                      Entries4 = sqmatrices[7]->GetEntries();
                      RowPtr4 = sqmatrices[7]->get_row_ptr();
                    }
                    if (n_matrices==2)
                    {
                      Entries5 = matrices[1]->GetEntries();
                      RowPtr5 = matrices[1]->get_row_ptr();                      
                    }
                  }

                  // for all dof in the mesh cell
                  // jj - ansatz function
                  for (jj=0;jj<N_;jj++)
                  {
                    dof_jj = DOF[jj];
                    // initialize the boundary integrals
                    for (l=0;l<2;l++)
                      integral[l] = 0;

                    // compute boundary integrals
                    // first component of the velocity
                    if (j==0)
                    {
                      for(l=0;l<N_LinePoints;l++)
                      {
                        // values of test functions in this quadrature point
                        JointValue = JointValues[l];
                        // get quadrature point on the boundary
                        x = x0 + 0.5*(x1-x0)*(qf1->get_point(l).x+1);
                        y = y0 + 0.5*(y1-y0)*(qf1->get_point(l).x+1);

                        //t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);

                        // weight times determinant of reference trafo
                        s = hE * qf1->get_weight(l);
                        // (A_11)_{ii,jj}
                        //Output::print("before ", integral[0], " ", integral[1]);
                        val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nx;
                        val += friction_parameter*JointValue[jj]*tx*JointValue[ii]*tx;
                        integral[0] += val*s;
                        // (A_12)_{ii,jj}
                        val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nx;
                        val+= friction_parameter*JointValue[jj]*ty*JointValue[ii]*tx;
                        integral[1] += val*s;
                      }
                      
                      // edge not parallel to y axis or penetration
                      if ((std::abs(ny)>eps)||(penetration_penalty<1e3))
                      {
                        // update first matrix
                        found = 0;
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                        {
                          if (ColInd1[ll] == dof_jj)
                          {
                            Entries1[ll] += integral[0];
                            found = 1;
                            break;
                          }
                        }
                        if (!found)
                        {
                          ErrThrow("ERROR A_11 ");
                        }
                        // update second matrix
                        if (n_sqmatrices>2)
                        {
                          found = 0;
                          for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                          {
                            if (ColInd2[ll] == dof_jj)
                            {
                              found = 1;
                              Entries2[ll] += integral[1];
                              break;
                            }
                          }
                          if (!found)
                          {
                            ErrThrow("ERROR A_12 ");
                          }
                        }
                      }
                      else                        // edge parallel to y-axis and no panetration
                      {
                        found = 0;
                        for (ll=0;ll<N_EdgeDOF; ll++)
                        {
                          if (dof_ii==DOF[EdgeDOF[ll]])
                          {
                            found =1;
                            break;
                          }
                        }
                        if (!found)
                          continue;
                        // update first matrix, set diagonal entry to 1
                        // all other entries to zero
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                        {
                          if (ColInd1[ll] == dof_ii)
                            Entries1[ll] = 1;
                          else
                            Entries1[ll] = 0;
                        }
                        // update second matrix, set all entries to zero
                        if (n_sqmatrices>2)
                        {
                          for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                            Entries2[ll] = 0;
                        }

                        if (n_sqmatrices==8)
                        {                         // M_11, set off diagonal to zero
                          for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          {
                            if (ColInd3[ll] != dof_ii)
                              Entries3[ll] = 0;
                          }
                          // M_12, set row to zero
                          for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                            Entries4[ll] = 0;
                        }

                        if (n_matrices==2)
                          for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                            Entries5[ll] = 0;

                        // set rhs to zero
                        RHS[dof_ii] = 0;
                      }
                    }                             // end first component (j==0)

                    // second component
                    if (j==1)
                    {
                      for(l=0;l<N_LinePoints;l++)
                      {
                        // values of test functions in this quadrature point
                        JointValue = JointValues[l];
                        // get quadrature point on the boundary
                        x = x0 + 0.5*(x1-x0)*(qf1->get_point(l).x+1);
                        y = y0 + 0.5*(y1-y0)*(qf1->get_point(l).x+1);

                        //t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                        // get velocity in this quadrature point

                        // weight times determinant of reference trafo
                        s = hE * qf1->get_weight(l);
                        // (A_21)_{ii,jj}
                        val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*ny;
                        val += friction_parameter*JointValue[jj]*ty*JointValue[ii]*tx;
                        integral[0] += s*val;
                        // (A_22)_{ii,jj}
                        val = penetration_penalty*JointValue[jj]*ny*JointValue[ii]*ny;
                        val += friction_parameter*JointValue[jj]*ty*JointValue[ii]*ty;
                        integral[1] += s*val;
                      }
                      //if ((integral[0]==0)&&( integral[1]==0)) continue;
                      //Output::print("2 ", integral[0], " ", integral[1]);

                      // edge not parallel to x-axis or pentration
                      if ((std::abs(nx)>eps)|| (penetration_penalty < 1e3))
                      {
                        if (n_sqmatrices>2)
                        {
                          // update first matrix
                          found = 0;
                          for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                          {
                            if (ColInd1[ll] == dof_jj)
                            {
                              Entries1[ll] += integral[0];
                              found =1 ;
                              break;
                            }
                          }
                          if (!found)
                          {
                            ErrThrow("ERROR A_21 ");
                          }
                        }

                        // update second matrix
                        found = 0;
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                        {
                          if (ColInd2[ll] == dof_jj)
                          {
                            Entries2[ll] += integral[1];
                            found =1 ;
                            break;
                          }
                        }
                        if (!found)
                        {
                          ErrThrow("ERROR A_22 ");
                        }
                      }
                      else                        // edge parallel to x-axis and no penetration
                      {
                        found = 0;
                        for (ll=0;ll<N_EdgeDOF; ll++)
                        {
                          if (dof_ii==DOF[EdgeDOF[ll]])
                          {
                            found =1;
                            break;
                          }
                        }
                        if (!found)
                          continue;

                        // update first matrix, set all entries to zero
                        if (n_sqmatrices>2)
                        {
                          for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                            Entries1[ll] = 0;
                        }
                        // update second matrix, set diagonal entry to 1,
                        // all other entries to 0
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                        {
                          if (ColInd2[ll] == dof_ii)
                          {
                            Entries2[ll] = 1;
                          }
                          else
                            Entries2[ll] = 0;
                        }

                        // set rhs to zero
                        RHS[dof_ii] = 0;

                        // update mass matrix
                        if (n_sqmatrices==8)
                        {
                          for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          {                       //M_22
                            if (ColInd3[ll] != dof_ii)
                              Entries3[ll] = 0;
                          }
                          // M_21
                          for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                            Entries4[ll] = 0;
                        }

                        if (n_matrices==2)
                          for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                            Entries5[ll] = 0;

                      }
                    }                             // end first component (j==1)
                  }                               // end inner loop over dof (jj)
                }                                 // end outer loop over dof (ii)

                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                break;                            // end slip with friction and penetration with resistance bc

                  default:
                    ErrThrow("This boundary condition is not handled here.");
            }                                     // endswitch Cond0
          }                                       // endif (Cond0==Cond1)
          else
          {
            ErrThrow("different boundary condition on one edge are not allowed!");
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)

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

  delete [] AuxArray[0];
  delete [] SecondDer;

  /*
    int N_Rows;
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

}    
 // end of Assemble


/*
  Assemble2D_VectFE:
    Assemble for vector finite elements (Raviart-Thomas, Brezzi-Douglas-Marini)
    Need the global orientation of normal at each inner edge/face

    implementation: Alfonso (07.09.2010)
*/

#ifdef __2D__
void Assemble2D_VectFE(int n_fespaces, const TFESpace2D** fespaces,
                       int n_sqmatrices, TSquareMatrix2D** sqmatrices,
                       int n_matrices, TMatrix2D** matrices, int n_rhs,
                       double** rhs, const TFESpace2D** ferhs,
                       LocalAssembling2D& la,
                       BoundCondFunct2D** BoundaryConditions, 
                       BoundValueFunct2D * const * const BoundaryValues)
{
  int N_AllMatrices = n_sqmatrices+n_matrices;
  std::vector<double> AbsDetjk(MaxN_QuadPoints_2D, 1.);
  double *righthand;
  
  
  // check if hanging nodes exist, I don't know how this works, we quit the
  // program in that case.
  for (int iSpace=0; iSpace<n_fespaces; iSpace++)
  {
    if(fespaces[iSpace]->get_n_hanging() != 0)
    {
      ErrThrow("Assemble2D_VectFE: hanging entries not supported. Exiting");
    }
  }
  bool *SecondDer = la.GetNeeds2ndDerivatives();
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  
  // ########################################################################
  // loop over all cells
  // ########################################################################
  auto Coll = fespaces[0]->GetCollection();// all spaces use same Coll
  int N_Cells = Coll->GetN_Cells(); // number of cells in this collection
  // set cell indices for all cells
  for(int icell=0; icell<N_Cells; icell++)
    Coll->GetCell(icell)->SetCellIndex(icell);
  for(int icell=0; icell<N_Cells; icell++)
  {
    double **LocRhs = nullptr;
    if(n_rhs)
    {
      LocRhs = new double* [n_rhs];
      righthand = new double [n_rhs*MaxN_BaseFunctions2D];
      memset(righthand,0,sizeof(double)*n_rhs*MaxN_BaseFunctions2D);
      for(int i=0;i<n_rhs;i++)
        LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
    }
    TBaseCell *cell = Coll->GetCell(icell); // current cell
    int N_Edges  = cell->GetN_Edges(); // number of edges of this cell
    
    cell->SetNormalOrientation();
    // ########################################################################
    // find local used elements on this cell
    // ########################################################################
    std::vector<const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
    std::vector<const BaseFunctions*> LocBF(n_fespaces);
    for(int iSpace=0; iSpace<n_fespaces; iSpace++)
    {
      const FiniteElement& fe = fespaces[iSpace]->get_fe(icell);
      LocalUsedElements[iSpace] = &fe;
      LocBF[iSpace] = fe.GetBaseFunct();
    }
    
    // ########################################################################
    // calculate values on original element
    // ########################################################################
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    
    // ########################################################################
    // assemble local matrices and right hand sides
    // ########################################################################
    // maximum number of used basis functions
    //int max_n_BF = *max_element(LocN_BF.begin(),LocN_BF.end()); 
    // for every matrix we allocate a local matrix with corresponding number of 
    // rows and columns
    double ***LocMatrices = nullptr;
    if(N_AllMatrices)
    {
      LocMatrices = new double**[N_AllMatrices];
      for(int i=0;i<N_AllMatrices;i++)
      {
        int n_rows = LocBF[la.rowSpaceOfMat(i)]->GetDimension();//number of rows
        int n_cols = LocBF[la.colSpaceOfMat(i)]->GetDimension();//number of columns
        LocMatrices[i] = new double*[n_rows];
        for(int j=0; j<n_rows; j++)
        {
          LocMatrices[i][j] = new double[n_cols]; 
          memset(LocMatrices[i][j], 0, sizeof(double)*n_cols);
        }
      }
    }                                               // endif N_AllMatrices
    
    la.GetLocalForms(qf_orig, LocBF, cell, icell, N_AllMatrices, n_rhs,
                     LocMatrices, LocRhs);
    
    // ########################################################################
    // add local matrices to global matrices (ansatz == test)
    // ########################################################################
    for(int iSqMat=0;iSqMat<n_sqmatrices;iSqMat++)
    {
      // fe space for this square matrix
      const TFESpace2D *fespace = fespaces[la.rowSpaceOfMat(iSqMat)];
      // the number of local basis functions (= size of local matrix)
      int N_BaseFunctions = LocBF[la.rowSpaceOfMat(iSqMat)]->GetDimension();
      
      double **Matrix = LocMatrices[iSqMat];
      int ActiveBound = fespace->get_n_active();
      const int *DOF = fespace->GetGlobalDOF(icell);
      
      // add local matrix to global
      for(int irow = 0; irow < N_BaseFunctions; irow++)
      {
        int RowDOF = DOF[irow];
        if(RowDOF<ActiveBound)
        { // active degree of freedom
          for(int icolumn=0;icolumn<N_BaseFunctions;icolumn++)
          {
            int columnDOF=DOF[icolumn];
            sqmatrices[iSqMat]->add(RowDOF,columnDOF,Matrix[irow][icolumn]);
          }
        }
        else
        { // nonactive degree of freedom (Dirichlet)
          sqmatrices[iSqMat]->set(RowDOF,RowDOF,1.0); // 1 on diagonal
        }
      }                                           // endfor m
    }                                             // endfor j
    
    // ########################################################################
    // add local matrices to global matrices (ansatz != test)
    // ########################################################################
    for(int iMat=0;iMat<n_matrices;iMat++)
    {
      auto testSpace = matrices[iMat]->GetTestSpace2D();
      auto ansatzSpace = matrices[iMat]->GetAnsatzSpace2D();
      auto test_fe = testSpace->get_fe(icell);
      auto ansatz_fe = ansatzSpace->get_fe(icell);
      
      // number of test and ansatz functions
      int N_Test = test_fe.GetN_DOF();
      int N_Ansatz = ansatz_fe.GetN_DOF();
      
      double **Matrix = LocMatrices[iMat+n_sqmatrices];
      
      const int *TestDOF = testSpace->GetGlobalDOF(icell);
      const int *AnsatzDOF = ansatzSpace->GetGlobalDOF(icell);
      
      int ActiveBound = testSpace->get_n_active();
      
      // add local matrix to global
      for(int irow = 0; irow < N_Test; irow++)
      {
        int rowDOF = TestDOF[irow];
        if(rowDOF<ActiveBound)
        {
          for(int icolumn=0; icolumn<N_Ansatz; icolumn++)
          {
            int columnDOF = AnsatzDOF[icolumn];
            matrices[iMat]->add(rowDOF,columnDOF,Matrix[irow][icolumn]);
          }
        }
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)
    
    if(N_AllMatrices)
    {
      for(int i=0;i<N_AllMatrices;i++)
      {
        int n_rows = LocBF[la.rowSpaceOfMat(i)]->GetDimension();//number of rows
        for(int j=0; j<n_rows; j++)
        {
          delete [] LocMatrices[i][j];
        }
        delete [] LocMatrices[i];
      }
      delete [] LocMatrices;
    }
    // ########################################################################
    // add local right-hand sides to global right-hand side
    // ########################################################################
    for(int irhs=0;irhs<n_rhs;irhs++)
    {
      const TFESpace2D *fespace = ferhs[irhs];
      const FiniteElement& fe = fespace->get_fe(icell);

      int N_BaseFunctions = fe.GetN_DOF();

      double *local_rhs = LocRhs[irhs];
      double *RHS = rhs[irhs];
      int ActiveBound = fespace->get_n_active();
      
      // dof of the rhs nodes connected to this cell
      const int *DOF = fespace->GetGlobalDOF(icell);

      // add local right-hand side to the global one
      for(int irow = 0; irow < N_BaseFunctions; irow++)
      { 
        int rowDOF = DOF[irow];
        if(rowDOF<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[rowDOF] += local_rhs[irow];
        }                                         // endif l
      }                                           // endfor m
      //////////////////////////////////////////////////////////////////////
      // take care of boundary conditions:      
      BoundCondFunct2D *BoundaryCondition = BoundaryConditions[irhs];
      BoundValueFunct2D * const BoundaryValue = BoundaryValues[irhs];
      
      for(int ijoint=0; ijoint<N_Edges; ijoint++)
      {
        const TJoint *joint = cell->GetJoint(ijoint);
        if(joint->GetType() == BoundaryEdge)
        {
          const TBoundEdge *boundedge = (const TBoundEdge *)joint;
          double t0,t1;
          boundedge->GetParameters(t0, t1);
          // get id of the boundary component
          int comp = boundedge->GetBoundComp()->GetID();
          // get type of the boundary condition in the middle of the edge
          BoundCond Cond;
          BoundaryCondition(comp, (t0+t1)/2, Cond);
          switch(Cond)
          {
            case DIRICHLET:
            {
              const FEDescriptor *FEDesc_Obj = fe.GetFEDesc();
              int N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              // if DG
              if (N_EdgeDOF==0)
                break;
              auto nf = fe.GetNodalFunctional();
              // number of points used for computation of nodal functionals
              int N_EdgePoints; 
              // points used for computation of nodal functionals
              const double *EdgePoints; 
              nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

              // read boundary values for each point
              std::vector<double> PointValues(N_EdgePoints);
              for(int iedgePoint=0;iedgePoint<N_EdgePoints;iedgePoint++)
              {
                double s = EdgePoints[iedgePoint];
                s = 0.5*(t0*(1-s) + t1*(1+s)); // map s from [-1,1] to [t0,t1]
                BoundaryValue(comp, s, PointValues[iedgePoint]);
              }
              // compute values for each dof on the boundary edge with the 
              // nodal functionals
              std::vector<double> FunctionalValues(N_EdgeDOF);
              nf->GetEdgeFunctionals(Coll, cell, ijoint, PointValues.data(),
                                     FunctionalValues.data());
              int *EdgeDOF = FEDesc_Obj->GetJointDOF(ijoint);
              // save boundary values of each dof on the boundary
              // edge in the rhs
              for(int l=0;l<N_EdgeDOF;l++)
              {
                RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              }
              break;
            }
            case NEUMANN:
            {
              // Basis functions
              auto BaseFunct = fe.GetBaseFunct();
              // get polynomial degree of fe
              int polynomialDegree = BaseFunct->GetPolynomialDegree();
              // get a suitable line quadrature formula
              auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                  2*polynomialDegree, BFRefElements::BFUnitLine);
              unsigned int N_LinePoints = qf1->GetN_QuadPoints();
              double **JointValues=FEDatabase::GetJointDerivatives2D(
                *BaseFunct, *qf1, ijoint, MultiIndex2D::D00);
              // get vertices of boundary edge
              double x0, x1, y0, y1;
              cell->GetVertex(ijoint)->GetCoords(x0, y0);
              cell->GetVertex((ijoint+1) % N_Edges)->GetCoords(x1, y1);
              // compute (half of the) length of the boundary edge
              double hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
              double nx = (y1-y0)/(2*hE);
              double ny = (x0-x1)/(2*hE);
              // for Piola transform
              double *JointValuesTransformed = 
                    new double [N_BaseFunctions * BaseFunct->GetBaseVectDim()];
              // compute boundary integral
              for(unsigned int ilinePoint=0; ilinePoint<N_LinePoints; ilinePoint++)
              {
                double zeta = qf1->get_point(ilinePoint).x;
                // values of test functions in this quadrature point
                double *JointValue = JointValues[ilinePoint];
                FEDatabase::GetOrigValues(
                    fe.GetRefTransID(), zeta, BaseFunct, ijoint, Coll, cell,
                    JointValue, nullptr,nullptr, JointValuesTransformed, nullptr,nullptr); 
                // get quadrature point on the boundary
                double t = t0 + 0.5*(t1-t0)*(zeta+1);
                double s;
                // get value in this quadrature point (in s)
                BoundaryValue(comp, t, s);
                // multiply value with weights from quadrature formula
                // and determinant from integral transformation to the
                // unit edge (-1,1)
                s *= hE * qf1->get_weight(ilinePoint);
                // in case of the pressure right hand side, the values in 
                // JointValuesTransformed[k+N_] are not valid. Therefore we 
                // need HOMOGENEOUS Neumann boundary conditions for the 
                // pressure function
                if(s==0.0)
                  continue;
                
                // update rhs for all test functions
                for(int k=0;k<N_BaseFunctions;k++)
                {
                  int rowDOF = DOF[k];
                  if(rowDOF < ActiveBound)
                  {
                    RHS[rowDOF] -= s*(JointValuesTransformed[k]*nx 
                                 +JointValuesTransformed[k+N_BaseFunctions]*ny);
                  }
                }
              }
              
              delete [] JointValuesTransformed;
              break;
            }
            case ROBIN:
              /*// get polynomial degree of fe
              l = fe.GetBaseFunct2D()->GetPolynomialDegree();
              // get a suitable line quadrature formula
              qf1 = QuadratureFormulaDatabase::qf_from_degree(
                  2*l, BFRefElements::BFUnitLine);
              N_LinePoints = qf1->GetN_QuadPoints();
              ele.GetBaseFunct2D()->MakeRefElementData(qf1->get_type());
              JointValues=FEDatabase::GetJointDerivatives2D(
                ele.GetBaseFunct2D_ID(), qf1->get_type(), m, D00);
              // get vertices of boundary edge
              cell->GetVertex(m)->GetCoords(x0, y0);
              cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
              // compute (half of the) length of the boundary edge
              hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
              // compute boundary integral
              for(l=0;l<N_LinePoints;l++)
              {
                // values of test functions in this quadrature point
                JointValue = JointValues[l];
                // get quadrature point on the boundary
                t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                // get value in this quadrature point (in s)
                BoundaryValue(comp, t, s);
                // multiply value with weights from quadrature formula
                // and determinant from integral transformation to the
                // unit edge (-1,1)
                s *= hE * qf1->get_weight(l);
                // update rhs for all test functions
                for(k=0;k<N_BaseFunctions;k++)
                  if((l3 = DOF[k])<ActiveBound)
                    RHS[l3] += s*JointValue[k];
              }*/
              ErrThrow("Robin boundary conditions not yet supported. Exiting");
              break;
            case PERIODIC:
              break;
            default:
              ErrThrow("Unknown boundary condition !");
              break;
          }                                     // endswitch Cond0
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
    if(n_rhs)
    {
      delete [] LocRhs[0]; 
      delete [] LocRhs;
    }
  } // endfor i (N_Cells)
} // end of Assemble_VectFE
#endif // __2D__



// =======================================================================
//
// Assemble2D_CIP
//
// assembling for continuous interior penalty discretization
//
// =======================================================================
#ifdef __2D__
void Assemble2D_CIP(const CoeffFct2D& Coeff, int n_fespaces,
                    const TFESpace2D **fespaces, int n_sqmatrices,
                    TSquareMatrix2D **sqmatrices, int n_matrices,
                    TMatrix2D ** /*matrices*/, int n_rhs, double **rhs,
                    TFESpace2D **ferhs, BoundCondFunct2D **BoundaryConditions,
                    BoundValueFunct2D **BoundaryValues,
                    TAuxParam2D *Parameters)
{
  const int MaxN_BaseFunctions2D_Ersatz =100;

  double w,integrant = 0.,tau_par,sigma_par;
  int N_AllMatrices = n_sqmatrices+n_matrices,out;
    int i,j,k,l,l3,n,m,r,q,dummy,N_UsedElements,ii,jj,ll,weak;
  int N_Cells, N_Points, N_Parameters, N_Edges, N_;
  unsigned int N_Points1D;
  int N_Joints, ref_n = 0;
  int Used[N_FEs2D];
  FE_type *UsedElements, CurrentElement;
  BaseFunction_type BaseFunctCell;
  BFRefElements bf2Drefelements;
  const double *X, *Y;
  double *Param[MaxN_QuadPoints_2D];
  double *righthand = nullptr;
  double **Matrices = nullptr, *aux, *aux2, *aux4;
  double ***LocMatrices = nullptr, **LocRhs = nullptr;
  double *Coeffs[MaxN_QuadPoints_2D];
  const int *DOF;
  int ActiveBound, end;
  double *Entries,*Entries1;
  const int *ColInd, *RowPtr;
  const int *ColInd1, *RowPtr1;
  double *RHS;
  const TBoundComp *BoundComp;
  double t0, t1, t, s,integral;
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  unsigned int N_LinePoints;
  double x0, x1, y0, y1, hE, nx, ny, eps=1e-12;
  bool *SecondDer;

  double *Coefficients1D[MaxN_QuadPoints_2D];
  
  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  ReferenceTransformation_type RefTrans;
  double*** value_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** xderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** yderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double *value_basefunct_ori[6];
  double *xderiv_basefunct_ori[6];
  double *yderiv_basefunct_ori[6];
  double x_pos_ref[6];
  double y_pos_ref[6];
  double x_pos[6];
  double y_pos[6];
  double *value_basefunct_ori_neigh[6];
  double *xderiv_basefunct_ori_neigh[6];
  double *yderiv_basefunct_ori_neigh[6];
  double x_pos_neigh[6];
  double y_pos_neigh[6];

  int neigh_edge;
  int N_Neigh;
  BaseFunction_type BaseFunctNeigh;
  ReferenceTransformation_type RefTransNeigh;
  const int *DOF_neigh;
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

#ifdef __3D__
  double z0, z1;
#endif

  out=2;
  
  if(out > 1)
    Output::print("CIP ", TDatabase::TimeDB->CURRENTTIME);

  // ########################################################################
  // store information in local arrays
  // ########################################################################

  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D_Ersatz*MaxN_BaseFunctions2D_Ersatz];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D_Ersatz];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D_Ersatz;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D_Ersatz;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D_Ersatz;
  }                                               // endif N_AllMatrices

  SecondDer = new bool[n_fespaces];
  SecondDer[0] = false;

  for (i=0;i<N_BaseFuncts2D;i++)
  {
    value_basefunct_ref1D[i] = new double* [6];
    xderiv_basefunct_ref1D[i] = new double* [6];
    yderiv_basefunct_ref1D[i] = new double* [6];
    for (j=0;j<6;j++)
    {

      value_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      xderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      yderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];

      memset( value_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( xderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( yderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );

    }
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
    xietaval_ref1D[i] = new double** [4];
    xideriv_ref1D[i] = new double** [4];
    etaderiv_ref1D[i] = new double** [4];
    for (j=0;j<4;j++)
    {
      xietaval_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      xideriv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      etaderiv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      for (n=0;n<MaxN_QuadPoints_1D;n++)
      {
        xietaval_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        xideriv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        etaderiv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];

        memset( xietaval_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( xideriv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( etaderiv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      }
    }
  }

  memset(Used, 0, N_FEs2D*sizeof(int));

  for(i=0;i<n_fespaces;i++)
  {
    auto fespace = fespaces[i];                        /* fe space */
    n = fespace->GetN_UsedElements();             /* # used finite elements */
    auto UsedElements = fespace->GetUsedElements(); /* used finite elements */
    for(j=0;j<n;j++)                              /* for all finite elements */
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
    }                                             // enfor j
  }                                               // endfor i

  N_UsedElements = 0;                             /* compute number of used elements */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE_type[N_UsedElements];        /* store used finite elements */
  j=0;                                            /* in array */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
  {
    UsedElements[j] = ( FE_type )i;
    j++;
  }                                               // endif

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################
  if (out==2)
      Output::print("N_UsedElements:", N_UsedElements);
  // Output::print("N_BaseFuncts2D:", N_BaseFuncts2D);
  // Output::print("MaxN_QuadPoints_1D:", MaxN_QuadPoints_1D);
  // Output::print("MaxN_BaseFunctions2D_Ersatz:", MaxN_BaseFunctions2D_Ersatz);

  for(n=0;n<N_UsedElements;n++)                   // for used finite elements
  {
    CurrentElement = UsedElements[n];
    const FiniteElement element(CurrentElement);
    l = element.GetBaseFunct()->GetPolynomialDegree();
    auto qf1D = QuadratureFormulaDatabase::qf_from_degree(
        2*l, BFRefElements::BFUnitLine);
    N_Points1D = qf1D->GetN_QuadPoints();
    BaseFunctCell = element.GetBaseFunct_ID();
                                                  // get base functions
    auto bf = element.GetBaseFunct();
    bf2Drefelements = bf->GetRefElement();
    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
      // quadrilateral cell
      case BFRefElements::BFUnitSquare:         // edge 0
        bf->GetDerivatives(MultiIndex2D::D00, -1, 1, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D10, -1, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, -1, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = -1;
        y_pos_ref[0] = 1;
        bf->GetDerivatives(MultiIndex2D::D00, 1, -1, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(MultiIndex2D::D10, 1, -1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 1, -1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] =-1;
        bf->GetDerivatives(MultiIndex2D::D00, 1, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(MultiIndex2D::D10, 1, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 1, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 1;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(MultiIndex2D::D00, -1, -1, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(MultiIndex2D::D10, -1, -1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, -1, -1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = -1;
        y_pos_ref[3] = -1;

        ref_n=4;

        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          xi1D[BaseFunctCell][0][j] = zeta;
          eta1D[BaseFunctCell][0][j] = -1;
          bf->GetDerivatives(MultiIndex2D::D00, zeta, -1, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(MultiIndex2D::D10, zeta, -1, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(MultiIndex2D::D01, zeta, -1, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          xi1D[BaseFunctCell][1][j] = 1;
          eta1D[BaseFunctCell][1][j] = zeta;
          bf->GetDerivatives(MultiIndex2D::D00, 1, zeta, xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(MultiIndex2D::D10, 1, zeta, xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(MultiIndex2D::D01, 1, zeta, etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          xi1D[BaseFunctCell][2][j] = -zeta;
          eta1D[BaseFunctCell][2][j] = 1;
          bf->GetDerivatives(MultiIndex2D::D00, -zeta, 1, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(MultiIndex2D::D10, -zeta, 1, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(MultiIndex2D::D01, -zeta, 1, etaderiv_ref1D[BaseFunctCell][2][j]);
        }                                         // edge 3
        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          xi1D[BaseFunctCell][3][j] = -1;
          eta1D[BaseFunctCell][3][j] = -zeta;
          bf->GetDerivatives(MultiIndex2D::D00, -1, -zeta, xietaval_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(MultiIndex2D::D10, -1, -zeta, xideriv_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(MultiIndex2D::D01, -1, -zeta, etaderiv_ref1D[BaseFunctCell][3][j]);
        }
        break;

      case BFRefElements::BFUnitTriangle:            // triangular cell

        bf->GetDerivatives(MultiIndex2D::D00, 0, 0, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D10, 0, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 0, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = 0;
        y_pos_ref[0] = 0;
        bf->GetDerivatives(MultiIndex2D::D00, 1, 0, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(MultiIndex2D::D10, 1, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 1, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] = 0;
        bf->GetDerivatives(MultiIndex2D::D00, 0, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(MultiIndex2D::D10, 0, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 0, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 0;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(MultiIndex2D::D00, 0.5, 0, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(MultiIndex2D::D10, 0.5, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 0.5, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = 0.5;
        y_pos_ref[3] = 0;
        bf->GetDerivatives(MultiIndex2D::D00, 0.5, 0.5, value_basefunct_ref1D[BaseFunctCell][4]);
        bf->GetDerivatives(MultiIndex2D::D10, 0.5, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 0.5, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[4] = 0.5;
        y_pos_ref[4] = 0.5;
        bf->GetDerivatives(MultiIndex2D::D00, 0, 0.5, value_basefunct_ref1D[BaseFunctCell][5]);
        bf->GetDerivatives(MultiIndex2D::D10, 0, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(MultiIndex2D::D01, 0, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[5] = 0;
        y_pos_ref[5] = 0.5;
        
        ref_n = 6;

        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature poin
        {
          double zeta = qf1D->get_point(j).x;
          Output::print(zeta, " ", BaseFunctCell);
          xi1D[BaseFunctCell][0][j] = (zeta+1)/2;
          eta1D[BaseFunctCell][0][j] = 0;
          bf->GetDerivatives(MultiIndex2D::D00, (zeta+1)/2, 0, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(MultiIndex2D::D10, (zeta+1)/2, 0, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(MultiIndex2D::D01, (zeta+1)/2, 0, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          Output::print(zeta);
          xi1D[BaseFunctCell][1][j] = (-zeta+1)/2;
          eta1D[BaseFunctCell][1][j] = (zeta+1)/2;
          bf->GetDerivatives(MultiIndex2D::D00, (-zeta+1)/2, (zeta+1)/2, xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(MultiIndex2D::D10, (-zeta+1)/2, (zeta+1)/2, xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(MultiIndex2D::D01, (-zeta+1)/2, (zeta+1)/2, etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (unsigned int j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          double zeta = qf1D->get_point(j).x;
          Output::print(zeta);
          xi1D[BaseFunctCell][2][j] = 0;
          eta1D[BaseFunctCell][2][j] = (-zeta +1)/2;
          bf->GetDerivatives(MultiIndex2D::D00, 0, (-zeta+1)/2, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(MultiIndex2D::D10, 0, (-zeta+1)/2, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(MultiIndex2D::D01, 0, (-zeta+1)/2, etaderiv_ref1D[BaseFunctCell][2][j]);
        }
        break;
      case BFRefElements::BFUnitTetrahedron :        // edge 0
      case BFRefElements::BFUnitHexahedron :        // edge 0
      case BFRefElements::BFUnitLine:
        {
          ErrThrow("This assemble routine is for 2D not for 3D");
        }
    }
  }                                               // endfor n
  if (out==2)
      Output::print("basefunct");

  for(l=0;l<ref_n;l++)
  {
    value_basefunct_ori[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    value_basefunct_ori_neigh[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];

  }

  for(m=0;m<4;m++)                                // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[m] = new double[N_Points1D];              // coordinates of edge i
    Y1D[m] = new double[N_Points1D];
                                                  // determinant of affine mapping
    for (unsigned int j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
    }

  }                                               // endfor m

  for (unsigned int j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################

  if (out==2)
      Output::print("parameters");
  N_Parameters = Parameters->GetN_Parameters();   // get number of parameters of equation
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coeffs[j] = aux2 + j*20;

  aux4 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coefficients1D[j] = aux4 + j*20;

  // ########################################################################
  // prepare loop over cells
  // ########################################################################

  // all spaces use same Coll
  auto Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    auto cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  TQuadFormula qf_ref_neigh(qf_ref);
  TQuadFormula qf_orig_neigh(qf_ref);
  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    auto cell = Coll->GetCell(i);                      // next cell

    // if (out==2) Output::print("cell ", i);
    for(n=0;n<n_sqmatrices;n++)
    {
      // calculate all needed derivatives of this FE function
      auto fespace = sqmatrices[n]->GetFESpace2D();
      auto element = fespace->get_fe(i); // finite element on cell
      CurrentElement = element.GetID();

      BaseFunctCell = element.GetBaseFunct_ID(); // basis functions
      N_ = element.GetN_DOF();           // # basis functions
      DOF = fespace->GetGlobalDOF(i);  // dof of current mesh cell

      SecondDer[0] = false;
      std::vector<const FiniteElement*> used_elements(1, &element);
      RefTrans = FEDatabase::GetOrig(used_elements, Coll, cell, SecondDer,
                                     qf_ref, qf_orig);
      N_Points = qf_orig.GetN_QuadPoints();
      X = qf_orig.get_xi();
      Y = qf_orig.get_eta();
      if(N_Parameters>0)                          // get parameters of equ.
        Parameters->GetParameters(qf_orig, i, Param);

                                                  // get coefficients of pde
      if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);
      // prepare 1D quadrature formula
      l = element.GetBaseFunct()->GetPolynomialDegree();
      if(out>2){ Output::print("Polynomial degree on cell: ", l);}
      auto qf1D = QuadratureFormulaDatabase::qf_from_degree(
          2*l, BFRefElements::BFUnitLine);
      N_Points1D = qf1D->GetN_QuadPoints();

      BoundaryCondition = BoundaryConditions[n];
      BoundaryValue = BoundaryValues[n];
      
      if(out>2)
      {
        for(unsigned int j=0;j<N_Points1D; j++)
        {
          Output::print("weights1D[", j, "]:", qf1D->get_weight(j));
        }
      }

                                                  // update data base
      N_Edges=cell->GetN_Edges();                 // # edges

      if(out>2)
      {
        for(r=0;r<N_Edges;r++)
        {
          cell->GetVertex(r)->GetCoords(x0, y0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
          if(out>2){Output::print("Local edge r: ", r, " Ecke A ", x0, " ", y0, " Ecke B ", x1, " ", y1);}
        }
      }

      for(r=0;r<N_Edges;r++)                      // loop over all edges of cell
      {                                           // get original coordinates of edge quad. points
        FEDatabase::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][r],
                                   eta1D[BaseFunctCell][r], X1D[r], Y1D[r]);

        for(unsigned int j=0;j<N_Points1D;j++)                 // get values and derivatives in original cell
        {
          FEDatabase::GetOrigValues(RefTrans, xi1D[BaseFunctCell][r][j],
            eta1D[BaseFunctCell][r][j],
            element.GetBaseFunct(),
            Coll, (const TGridCell *)cell,
            xietaval_ref1D[BaseFunctCell][r][j],
            xideriv_ref1D[BaseFunctCell][r][j],
            etaderiv_ref1D[BaseFunctCell][r][j],
            xyval_ref1D[r][j],
            xderiv_ref1D[r][j],
            yderiv_ref1D[r][j]);
        }

      }                                           // endfor r

      FEDatabase::GetOrigFromRef(RefTrans,ref_n,x_pos_ref,y_pos_ref,x_pos,
                                 y_pos);
      for(l=0;l<ref_n;l++)
      {
        FEDatabase::GetOrigValues(RefTrans, x_pos_ref[l],
          y_pos_ref[l],
          element.GetBaseFunct(),
          Coll, (const TGridCell *)cell,
          value_basefunct_ref1D[BaseFunctCell][l],
          xderiv_basefunct_ref1D[BaseFunctCell][l],
          yderiv_basefunct_ref1D[BaseFunctCell][l],
          value_basefunct_ori[l],
          xderiv_basefunct_ori[l],
          yderiv_basefunct_ori[l]);
        // Output::print("Hallo: x_pos_ref[l]: ", x_pos_ref[l], "value_basefunct_ref1D[BaseFunctCell][l]: ", value_basefunct_ref1D[BaseFunctCell][l]);
      }

      for(r=0;r<N_Edges;r++)
      {                                           // For each edge, get the corresponding neighbour cell.
        auto neigh=cell->GetJoint(r)->GetNeighbour(cell);
        //#######################################################################//
        // get coefficients on edges
        //only implemented for coeffs that do not depend on the params
        //#######################################################################//

        //  if(N_Parameters>0)                // get parameters of equ.
        // Parameters->GetParameters(N_Points1D, Coll, cell, i, xi1D[BaseFunctCell][r], eta1D[BaseFunctCell][r], X1D[r], Y1D[r], Param1D);

        if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
        //#######################################################################//
        // If there is a neighbour to the edge, do...
        if(neigh)
        {                                         // Get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
          if(i<q)
          {

            // calculate all needed derivatives of this FE function
                                                  // finite element on neighbour
            auto eleNeigh = fespaces[n]->get_fe(q);
            BaseFunctNeigh = eleNeigh.GetBaseFunct_ID();
            N_Neigh = eleNeigh.GetN_DOF();       // number of basis functions on neighbour
                                                  // dof of current mesh cell on neighbour cell
            DOF_neigh = fespaces[n]->GetGlobalDOF(q);

            std::vector<const FiniteElement*> used_elements(1, &eleNeigh);
            RefTransNeigh = FEDatabase::GetOrig(used_elements, Coll, neigh,
                                                SecondDer, qf_ref_neigh,
                                                qf_orig_neigh);
                                                   
            /* To do: Include this for FE-Spaces of which the polynomial degree
            of the basis fuctions depend on the cell
            RefTrans_neigh = FEDatabase::GetOrig(N_LocalUsedElements,
                                                 LocalUsedElements_neigh,
                                                 Coll, neigh, SecondDer,
                                                 N_Points_neigh, xi_neigh,
                                                 eta_neigh, weights_neigh,
                                                 X_neigh, Y_neigh,
                                                 AbsDetjk_neigh);
            if(N_Parameters>0)                // get parameters of equ.
            Parameters->GetParameters(N_Points_neigh, Coll, neigh, q, xi_neigh, eta_neigh, X_neigh, Y_neigh, Param_neigh);

            if(Coeff)                               // get coefficients of pde in the neighbour cell
            Coeff(N_Points_neigh, X_neigh, Y_neigh, Param_neigh, );

            // prepare 1D quadrature formula in the neighbour cell
            l = eleNeigh.GetBaseFunct2D()->GetPolynomialDegree();
            qf1D_neigh = QuadratureFormulaDatabase::qf_from_degree(
                2*l, BFRefElements::BFUnitLine);
            qf1D_neigh->GetFormulaData(N_Points1D_neigh, weights1D_neigh, zeta_neigh)
            */
            // Get edge of the neighbour cell which is the edge r
            neigh_edge=0;
            while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;

            //RefTransNeigh= eleNeigh.GetRefTransID();          // reftrafo of neighbour
            //FEDatabase::SetCellForRefTrans(neigh,RefTransNeigh);
            //FEDatabase::GetBaseFunct2DFromFE2D(CurrEleNeigh)
            //->MakeRefElementData(LineQuadFormula);

            for(m=0;m<4;m++)                      // arrays for coordinates on neighbour cell
            {
              X1D_neigh[m] = new double[N_Points1D];
              Y1D_neigh[m] = new double[N_Points1D];
            }

            // get original coordinates of edge quad. points of neighbour cell
            FEDatabase::GetOrigFromRef(RefTransNeigh, N_Points1D, 
                                       xi1D[BaseFunctNeigh][neigh_edge],
                                       eta1D[BaseFunctNeigh][neigh_edge],
                                       X1D_neigh[neigh_edge],
                                       Y1D_neigh[neigh_edge]);

            // get values and derivatives on original neighbour cell on edge neigh_edge
            for (unsigned int j=0;j<N_Points1D;j++)
            {                                     //for (k=0; k< MaxN_BaseFunctions2D_Ersatz; k++)
              // Output::print(" xi1D ", xi1D[BaseFunctCell][r][j], " eta1D ", eta1D[BaseFunctCell][r][j], " BaseFunct: ", BaseFunctCell);
              // Output::print(" xi1D ", xi1D[BaseFunctNeigh][neigh_edge][j], " eta1D ", eta1D[BaseFunctCell][r][j], " BaseFunctNeigh: ", BaseFunctNeigh);
              if(out>2){ Output::print("X1D[r][j]: ", X1D[r][j], " Y1D[r][j]: ", Y1D[r][j], " X1D[neigh_edge][j] ", X1D_neigh[neigh_edge][j], " Y1D[neigh_edge][j]: ", Y1D_neigh[neigh_edge][j]);}
            }

            if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0] )
            {
              if(out>2)
              {
                Output::print("Quadrature points on neighbour edge in the correct order.");
              }
              for (unsigned int j=0;j<N_Points1D;j++)
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
                  eta1D[BaseFunctNeigh][neigh_edge][j],
                  eleNeigh.GetBaseFunct(),
                  Coll, (TGridCell *)neigh,
                  xietaval_ref1D[BaseFunctNeigh][neigh_edge][j],
                  xideriv_ref1D[BaseFunctNeigh][neigh_edge][j],
                  etaderiv_ref1D[BaseFunctNeigh][neigh_edge][j],
                  xyval_refNeigh1D[j],
                  xderiv_refNeigh1D[j],
                  yderiv_refNeigh1D[j]);
              }                                   //endfor j

            }                                     //endif
            else
            {
              if(out>2)
              {
                Output::print("Inverse the order of the quadrature oints on neighbour edge !");
              }
              for (unsigned int j=0;j<N_Points1D;j++)
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
                  eta1D[BaseFunctNeigh][neigh_edge][j],
                  eleNeigh.GetBaseFunct(),
                  Coll, (TGridCell *)neigh,
                  xietaval_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  xideriv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  etaderiv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  xyval_refNeigh1D[j],
                  xderiv_refNeigh1D[j],
                  yderiv_refNeigh1D[j]);

              }                                   //endfor j
            }                                     //endelse

            FEDatabase::GetOrigFromRef(RefTransNeigh, ref_n, x_pos_ref,
                                       y_pos_ref, x_pos_neigh, y_pos_neigh);
            for(l=0;l<ref_n;l++)
            {
              FEDatabase::GetOrigValues(
                RefTrans, x_pos_ref[l], y_pos_ref[l],
                eleNeigh.GetBaseFunct(),
                Coll, (TGridCell *)neigh,
                value_basefunct_ref1D[BaseFunctNeigh][l],
                xderiv_basefunct_ref1D[BaseFunctNeigh][l],
                yderiv_basefunct_ref1D[BaseFunctNeigh][l],
                value_basefunct_ori_neigh[l],
                xderiv_basefunct_ori_neigh[l],
                yderiv_basefunct_ori_neigh[l]);
            }

            for(k=0;k<N_; k++)
            {
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    Output::print("Basisfkt: ", DOF[k], "  (x,y)-coordinate: ",
                                  x_pos[l], " ", y_pos[l], " value: ",
                                  value_basefunct_ori[l][k]);
                  }
                }
              }
            }                                     //endfor k

            for(k=0;k<N_Neigh; k++)
            {
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    Output::print("Basisfkt neigh: ", DOF_neigh[k],
                                  "  (x,y)-coordinate: ", x_pos_neigh[l], " ",
                                  y_pos_neigh[l], " value: ",
                                  value_basefunct_ori_neigh[l][k]);
                  }
                }
              }
            }                                     //endfor k

            //Compute the jumps of the basis functions
            //and of their derivatives in the quadrature points on edge r
            //First for the basis functions of cell i

            for(k=0;k<N_; k++)
            {
              dummy = 0;
              l=0;
              // Check if basis function k of cell i is in the FE-Space of neighbour cell q
              while(l<N_Neigh && dummy == 0)
              {
                if(DOF[k] == DOF_neigh[l])dummy=1;
                l++;
              }
              l = l-1;
              // if basis function k of cell i is in the local FE-Space of neighbour cell q do
              if(dummy ==1 )
              {                                   // Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
                for(unsigned int j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][k] = xyval_ref1D[r][j][k]  -  xyval_refNeigh1D[j][l];
                  jump_xderiv[j][k] = xderiv_ref1D[r][j][k] - xderiv_refNeigh1D[j][l];
                  jump_yderiv[j][k] = yderiv_ref1D[r][j][k] - yderiv_refNeigh1D[j][l];
                  if(out>2)
                  {
                    Output::print("x: ", xyval_ref1D[r][j][k],
                                  " xn: ", xyval_refNeigh1D[j][l],
                                  " j= ", jump_xyval[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r,
                                  " qp: ", j);
                    Output::print("xd: ", xderiv_ref1D[r][j][k],
                                  " xdn: ", xderiv_refNeigh1D[j][l],
                                  " jump= ", jump_xderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                    Output::print("yd: ", yderiv_ref1D[r][j][k],
                                  " ydn: ", yderiv_refNeigh1D[j][l],
                                  " j= ", jump_yderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                  }
                }
              }                                   //endif
              // if basis function k of cell i is NOT in the local FE-Space of neighbour cell q do
              if (dummy == 0)
              {
                for(unsigned int j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][k]  = xyval_ref1D[r][j][k] ;
                  jump_xderiv[j][k] = xderiv_ref1D[r][j][k];
                  jump_yderiv[j][k] = yderiv_ref1D[r][j][k];

                  if(out>2)
                  {
                    Output::print("x: ", xyval_ref1D[r][j][k],                                  
                                  " j= ", jump_xyval[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r,
                                  " qp: ", j);
                    Output::print("xd: ", xderiv_ref1D[r][j][k],                                  
                                  " jump= ", jump_xderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                    Output::print("yd: ", yderiv_ref1D[r][j][k],
                                  " j= ", jump_yderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                  }
                }                                 //endfor j
              }                                   //endif
            }                                     //endfor k

            // Then for the basis functions of neighbour cell q
            //which are not in the local FE-Space of cell i
            for(l=0;l<N_Neigh; l++)
            {
              dummy = 0;
              k=0;
              while(k<N_ && dummy == 0 )
              {
                if(DOF_neigh[l] == DOF[k]) dummy=1 ;
                k++;
              }
              k=k-1;
              // If basis function l of neighbour cell q is NOT  in the local FE-Space of cell i do

              if( dummy == 0)
              {

                for(unsigned int j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][l+N_] = -xyval_refNeigh1D[j][l] ;
                  jump_xderiv[j][l+N_]= -xderiv_refNeigh1D[j][l];
                  jump_yderiv[j][l+N_]= -yderiv_refNeigh1D[j][l];
                  if(out>2)
                  {
                    Output::print(
                                  " xn: ", xyval_refNeigh1D[j][l],
                                  " j= ", jump_xyval[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r,
                                  " qp: ", j);
                    Output::print(
                                  " xdn: ", xderiv_refNeigh1D[j][l],
                                  " jump= ", jump_xderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                    Output::print(
                                  " ydn: ", yderiv_refNeigh1D[j][l],
                                  " j= ", jump_yderiv[j][k],
                                  " bf: ", DOF[k], " c: ", i,
                                  " e: ", r, " qp: ", j);
                  }
                }                                 //endfor j
              }                                   //endif
            }                                     //endfor l

            // #################################################################################
            // Compute the edge integrals with the jumps of the basis functions and their derivatives
            // #################################################################################

            //[0][4]: " << Coeffs[0][3]);
            // get vertices of boundary edge
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            if(out>2){  Output::print("Ecke A ", x0, " ", y0, " Ecke B ", x1, " ", y1);}
            // compute length of the boundary edge
            hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            // tangential normal vector to this boundary (normalized)
            //double tx = (x1-x0)/hE;
            //double ty = (y1-y0)/hE;
            tau_par = TDatabase::ParamDB->FACE_SIGMA;
            tau_par = tau_par*hE*hE;
            if (out>2){ Output::print(" tau edge stabilization: ", tau_par);};

            Entries1 = sqmatrices[n]->GetEntries();
            RowPtr1 = sqmatrices[n]->get_row_ptr();
            ColInd1 = sqmatrices[n]->get_vector_columns();
            fespace = sqmatrices[n]->GetFESpace2D();
            ActiveBound = fespace->get_n_active();
            if(out>2)
            {
              Output::print("ActiveBound von sqmatrix ", n, "  :", ActiveBound);
              for(unsigned int j=0;j<N_Points1D; j++)
              {
                Output::print(" b1: ", Coefficients1D[j][5], " b2: ",
                              Coefficients1D[j][6]);
              }
            }
            //edge integrals: test function from cell <-> ansatz function from cell
            if(out>2){ Output::print("testfct. cell<-> ansatzfct. cell Integral:");}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              // Dirichlet node

              if (dof_ii>=ActiveBound)continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (unsigned int j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  // distinguish different types of the cip method
                  // 0 : Burman, Hansbo 2004, eq. (3)
                  // 1:  Burman, Hansbo 2004, eq. (5), without orthogonal term, ES in Section 3
                  // 2: do nothing
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*std::abs(nx*Coefficients1D[j][5] + ny*Coefficients1D[j][6])*
                        (/*Coefficients1D[j][5]**/jump_xderiv[j][jj]*nx 
                       + /*Coefficients1D[j][6]**/jump_yderiv[j][jj]*ny) *
                        (/*Coefficients1D[j][5]**/jump_xderiv[j][ii]*nx 
                       + /*Coefficients1D[j][6]**/jump_yderiv[j][ii]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj]*Coefficients1D[j][1] + jump_yderiv[j][jj]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii]*Coefficients1D[j][1] + jump_yderiv[j][ii]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }

                  w = qf1D->get_weight(j) * hE/2;
                  integral += w*integrant;        // integral on the edge
                }

                // update matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      Output::print("integral: ", integral, " DOF Testfkt: ",
                                    dof_ii, " DOF Ansatzfkt: ", dof_jj);
                    }

                    Entries1[ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  ErrThrow("ERROR in Assemble_edge_integrals (test function cell - ansatz function cell) ");
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)

            //edge integrals: test function from cell <-> ansatz function from neighbour
            if(out>2){ Output::print("testfct. cell<-> ansatzfct. neighbour Integral:");}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound)
                continue;
              for (jj=0;jj<N_Neigh;jj++)
              {
                dof_jj = DOF_neigh[jj];

                //################################################################################################
                //#### Check if ansatz  function jj  of the neighbour cell is in the FE-Space of the cell  ####
                dummy = 0;
                l=0;
                while(l<N_ && dummy == 0)
                {
                  if(dof_jj == DOF[l])dummy=1;
                  l++;
                }
                if (dummy == 1) continue;
                //################################################################################################
                integral=0;
                for (unsigned int j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  // distinguish different types of the cip method
                  // 0 : Burman, Hansbo 2004, eq. (3)
                  // 1:  Burman, Hansbo 2004, eq. (5), without orthogonal term, ES in Section 3
                  // 2: do nothing
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*std::abs(nx*Coefficients1D[j][5] + ny*Coefficients1D[j][6]) *
                        (/*Coefficients1D[j][5]**/jump_xderiv[j][jj+N_]*nx 
                       + /*Coefficients1D[j][6]**/jump_yderiv[j][jj+N_]*ny) *
                        (/*Coefficients1D[j][5]**/jump_xderiv[j][ii]*nx 
                       + /*Coefficients1D[j][6]**/jump_yderiv[j][ii]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj+N_]*Coefficients1D[j][1] + jump_yderiv[j][jj+N_]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii]*Coefficients1D[j][1] + jump_yderiv[j][ii]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }
                  w = qf1D->get_weight(j) * hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      Output::print("integral: ", integral, " DOF Testfkt: ",
                                    dof_ii, " DOF Ansatzfkt: ", dof_jj);
                    }

                    Entries1[ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  Output::print("ERROR in Assemble_edge_integrals (test "
                                "function cell - ansatz function neighbour)");
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)

            //edge integrals: test function from neighbour <-> ansatz function from cell
            if(out>2){ Output::print("testfct. neighbour<-> ansatzfct. cell Integral:");}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound)
                continue;

              //################################################################################################
              //#### Check if test function ii  of the neighbour cell is in the FE-Space of the cell    ####
              dummy = 0;
              l=0;

              while(l<N_ && dummy == 0)
              {
                if(dof_ii == DOF[l])dummy=1;
                l++;
              }
              if (dummy == 1) continue;
              //################################################################################################

              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];

                integral=0;
                for (unsigned int j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*std::abs(nx*Coefficients1D[j][5] + ny*Coefficients1D[j][6])*
                        ( /*Coefficients1D[j][5]**/jump_xderiv[j][jj]*nx 
                        + /*Coefficients1D[j][6]**/jump_yderiv[j][jj]*ny) *
                        ( /*Coefficients1D[j][5]**/jump_xderiv[j][ii+N_]*nx 
                        + /*Coefficients1D[j][6]**/jump_yderiv[j][ii+N_]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj]*Coefficients1D[j][1] + jump_yderiv[j][jj]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii+N_]*Coefficients1D[j][1] + jump_yderiv[j][ii+N_]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }
                  w = qf1D->get_weight(j) * hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      Output::print("integral: ", integral, " DOF testfct: ",
                                    dof_ii, " DOF ansatzfct: ", dof_jj);
                    }

                    Entries1[ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  Output::print("ERROR in Assemble_edge_integrals (test function neighbour - ansatz function cell)");
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)

            //edge integrals: test function from neighbour <-> ansatz function from neighbour
            if(out>2){ Output::print("testfct. neighbour <-> ansatzfct. neighbour Integral:");}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound) continue;

              //################################################################################################
              //#### Check if test function ii  of the neighbour cell is in the FE-Space of the cell     ####
              dummy = 0;
              l=0;
              while(l<N_ && dummy == 0)
              {
                if(dof_ii == DOF[l])dummy=1;
                l++;
              }
              if (dummy == 1) continue;
              //################################################################################################

              for (jj=0;jj<N_Neigh;jj++)
              {
                dof_jj = DOF_neigh[jj];

                //################################################################################################
                //####  Check if ansatz  function jj  of the neighbour cell is in the FE-Space of the cell  ####
                dummy = 0;
                l=0;
                while(l<N_ && dummy == 0)
                {
                  if(dof_jj == DOF[l])dummy=1;
                  l++;
                }
                if (dummy == 1) continue;
                //################################################################################################

                if(dummy == 0)                    //
                  // initialize the boundary integrals
                  integral=0;
                for (unsigned int j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*std::abs(nx*Coefficients1D[j][5] + ny*Coefficients1D[j][6])*
                        ( /*Coefficients1D[j][5]**/jump_xderiv[j][jj+N_]*nx 
                        + /*Coefficients1D[j][6]**/jump_yderiv[j][jj+N_]*ny) *
                        ( /*Coefficients1D[j][5]**/jump_xderiv[j][ii+N_]*nx 
                        + /*Coefficients1D[j][6]**/jump_yderiv[j][ii+N_]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj+N_]*Coefficients1D[j][5] + jump_yderiv[j][jj+N_]*Coefficients1D[j][6]) *
                        (jump_xderiv[j][ii+N_]*Coefficients1D[j][1] + jump_yderiv[j][ii+N_]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }

                  w = qf1D->get_weight(j) * hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      Output::print("integral: ", integral, " DOF testfct: ",
                                    dof_ii, " DOF ansatzfct: ", dof_jj);
                    }

                    Entries1[ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  Output::print("ERROR in Assemble_edge_integrals (test function neighbour - ansatz function neighbour)");
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)

            for (m=0;m<4;m++)
            {
              delete X1D_neigh[m];
              delete Y1D_neigh[m];
            }                                     //endfor m
          }                                       //endif i<q
        }                                         //endif neigh

        else
        {
          weak = TDatabase::ParamDB->WEAK_BC;
          if(weak==1)
          {
            Output::print("weak ", r);
            auto joint = cell->GetJoint(r);
            if(joint->GetType() == BoundaryEdge ||
              joint->GetType() == IsoBoundEdge)
            {
              if(joint->GetType() == BoundaryEdge)
              {
                auto boundedge = (const TBoundEdge *)joint;
                BoundComp = boundedge->GetBoundComp();
                boundedge->GetParameters(t0, t1);
              }
              else
              {
                auto isoboundedge = (const TIsoBoundEdge *)joint;
                BoundComp = isoboundedge->GetBoundComp();
                isoboundedge->GetParameters(t0, t1);
              }
              // get id of the boundary component
              Output::print("weak ", r);
              comp=BoundComp->GetID();
              Output::print("weak ", comp);
              // get type of the boundary condition at the beginning
              // and at the end of the current edge
              if (t0 < t1)
              {
                BoundaryCondition(comp, t0+eps, Cond0);
                BoundaryCondition(comp, t1-eps, Cond1);
              }
              else
              {
                BoundaryCondition(comp, t0-eps, Cond0);
                BoundaryCondition(comp, t1+eps, Cond1);
              }
              Output::print("weak ", r);
              // only one boundary condition per edge allowed
              if(Cond0 == Cond1)
              {
                Output::print("weak ", r);
                if(Cond0 != DIRICHLET)
                  continue;
              }
            }
            Output::print("weak ", r);
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            if(out>2){ Output::print("Ecke A ", x0, " ", y0, " Ecke B ", x1, " ", y1);}
            // compute length of the boundary edge
            hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            // tangential normal vector to this boundary (normalized)
            //double tx = (x1-x0)/hE;
            //double ty = (y1-y0)/hE;
            sigma_par = TDatabase::ParamDB->WEAK_BC_SIGMA;
            if(out>2){Output::print("weak_bc: tau edge stabilization: ", tau_par);};

            Entries1 = sqmatrices[n]->GetEntries();
            RowPtr1 = sqmatrices[n]->get_row_ptr();
            ColInd1 = sqmatrices[n]->get_vector_columns();
            fespace = sqmatrices[n]->GetFESpace2D();
            ActiveBound = fespace->get_n_active();
            if(out>2)
            {
              Output::print("ActiveBound von sqmatrix ", n, "  :", ActiveBound);
              for(unsigned int j=0;j<N_Points1D; j++)
              {
                Output::print(" b1: ", Coeffs[j][1], " b2: ", Coeffs[j][2]);
              }
            }

            //edge integrals: test function from cell <-> ansatz function from cell
            if(out>2){ Output::print("testfct. cell<-> ansatzfct. cell Integral:");}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              // Dirichlet node

              if (dof_ii>=ActiveBound)continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (unsigned int j=0;j<N_Points1D;j++)        // compute edge integral; Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
                {
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += Coeffs[j][0] * (-xyval_ref1D[r][j][ii]*(xderiv_ref1D[r][j][jj]*nx + yderiv_ref1D[r][j][jj]*ny) - xyval_ref1D[r][j][jj]*(xderiv_ref1D[r][j][ii]*nx + yderiv_ref1D[r][j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += sigma_par*Coeffs[j][0]*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj]/hE;
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) < 0) integrant -= (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                  w = qf1D->get_weight(j) * hE/2;
                  integral += w*integrant;        // integral on the edge
                }

                // update matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      Output::print("integral: ", integral, " DOF Testfkt: ",
                                    dof_ii, " DOF Ansatzfkt: ", dof_jj);
                    }

                    Entries1[ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  ErrThrow("ERROR in Assemble_edge_integrals (test function cell - ansatz function cell) ");
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)
          }                                       //endif weak
        }                                         //endfor r (loop over edges)
      }                                           //endfor r (loop over edges)
    }                                             // endfor n (loop over sqmatrices)

    if(weak>=1)
    {
      for(n=0;n<n_rhs;n++)
      {
        auto fespace = ferhs[n];
        ActiveBound = fespace->get_n_active();
        auto element = fespace->get_fe(i);
        CurrentElement = fespace->get_fe_type(i);

        N_ = element.GetN_DOF();
        if(out>2){Output::print("N_BaseFunct ", N_);}

        RHS = rhs[n];
        // find space for this linear form

        // dof of the rhs nodes connected to this cell
        DOF = fespace->GetGlobalDOF(i);
        if(out>2)
        {
          for(k=0; k<N_; k++)
          {
            Output::print("DOF[k]: k: ", k, " DOF: ", DOF[k]);
          }
        }

        BoundaryCondition = BoundaryConditions[n];
        BoundaryValue = BoundaryValues[n];

        N_Joints = cell->GetN_Edges();

        if(out>2)Output::print("N_Joints: ", N_Joints);

        std::vector<const FiniteElement*> used_elements(1, &element);
        RefTrans = FEDatabase::GetOrig(used_elements, Coll, cell, SecondDer,
                                       qf_ref, qf_orig);
        N_Points = qf_orig.GetN_QuadPoints();
        X = qf_orig.get_xi();
        Y = qf_orig.get_eta();
        if(N_Parameters>0)                        // get parameters of equ.
          Parameters->GetParameters(qf_orig, i, Param);
                                                  // get coefficients of pde
        if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);

        for(m=0;m<N_Joints;m++)
        {
          auto joint = cell->GetJoint(m);
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
            joint->GetType() == IsoBoundEdge)
          {
            if(joint->GetType() == BoundaryEdge||
            joint->GetType() == InterfaceJoint)
            {
              auto boundedge = (const TBoundEdge *)joint;
              BoundComp = boundedge->GetBoundComp();
              boundedge->GetParameters(t0, t1);
            }
            else
            {
              auto isoboundedge = (const TIsoBoundEdge *)joint;
              BoundComp = isoboundedge->GetBoundComp();
              isoboundedge->GetParameters(t0, t1);
            }
            // get id of the boundary component
            comp=BoundComp->GetID();
            // get type of the boundary condition at the beginning
            // and at the end of the current edge
            if (t0 < t1)
            {
              BoundaryCondition(comp, t0+eps, Cond0);
              BoundaryCondition(comp, t1-eps, Cond1);
            }
            else
            {
              BoundaryCondition(comp, t0-eps, Cond0);
              BoundaryCondition(comp, t1+eps, Cond1);
            }
            // only one boundary condition per edge allowed
            if(Cond0 == Cond1)
            {
              if(Cond0 == DIRICHLET)
              {
                // get polynomial degree of fe
                if(out>2){Output::print("Edge ", m);}
                l = element.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                auto qf1D = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                if(out>2){Output::print("QuadFormula ", qf1D);}
                N_LinePoints = qf1D->GetN_QuadPoints();

                FEDatabase::GetOrigFromRef(RefTrans, N_Points1D,
                                           xi1D[BaseFunctCell][m],
                                           eta1D[BaseFunctCell][m],
                                           X1D[m], Y1D[m]);

                for(unsigned int j=0;j<N_Points1D;j++)         // get values and derivatives in original cell
                {
                  FEDatabase::GetOrigValues(
                    RefTrans, xi1D[BaseFunctCell][m][j],
                    eta1D[BaseFunctCell][m][j],
                    element.GetBaseFunct(),
                    Coll, (TGridCell *)cell,
                    xietaval_ref1D[BaseFunctCell][m][j],
                    xideriv_ref1D[BaseFunctCell][m][j],
                    etaderiv_ref1D[BaseFunctCell][m][j],
                    xyval_ref1D[m][j],
                    xderiv_ref1D[m][j],
                    yderiv_ref1D[m][j]);
                }

                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute the length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                nx = (y1-y0)/hE;
                ny = (x0-x1)/hE;
                // tangential normal vector to this boundary (normalized)
                //double tx = (x1-x0)/hE;
                //double ty = (y1-y0)/hE;
                sigma_par = TDatabase::ParamDB->WEAK_BC_SIGMA;

                // compute boundary integral: ansatz function from cell with Dirichlet boundary value
                for(k=0;k<N_;k++)
                {
                  l3 = DOF[k];
                  integral=0;
                  for(unsigned int l=0;l<N_LinePoints;l++)
                  {
                    if(out>2){Output::print("l:", l, " normal vector: ", nx, " ", ny);};
                    double zeta = qf1D->get_point(l).x;
                    // get quadrature point on the boundary
                    t = t0 + 0.5*(t1-t0)*(zeta+1);
                    if(out>2){Output::print("Quad formula data: Quadrature point[l]:", zeta, " weight: ", qf1D->get_weight(l));};
                    // get value in this quadrature point (in s)
                    BoundaryValue(comp, t, s);
                    if(out>2)
                    {
                      Output::print("Parameter Quadrature Point: ", t,
                                    " Boundary Value: ", s);
                      Output::print("Function values Quadrature Point: ", k,
                                    " xyval ", xyval_ref1D[m][l][k], " xderiv ",
                                    xderiv_ref1D[m][l][k], " yderiv ",
                                    yderiv_ref1D[m][l][k]);

                    };
                    // multiply value with weights from quadrature formula
                    // and determinant from integral transformation to the
                    // unit edge (-1,1)
                    integrant = 0;
                                                  // contribution of diffusive term
                    integrant += Coeffs[j][0]*(-s*(xderiv_ref1D[m][l][k]*nx + yderiv_ref1D[m][l][k]*ny));
                                                  // penalty term to achieve coercitivity
                    integrant += sigma_par*Coeffs[j][0]*s*xyval_ref1D[m][l][k]/hE;
                    if((Coeffs[l][1]*nx + Coeffs[l][2]*ny) < 0) integrant -= (Coeffs[l][1]*nx + Coeffs[l][2]*ny)*s*xyval_ref1D[m][l][k];
                    integral += hE/2 * qf1D->get_weight(l)*integrant;
                  }
                  if(l3 < ActiveBound) RHS[l3] += integral;
                  else
                  {
                    ErrThrow("Index l3 bigger than ActiveBound!");
                  }
                  if(out>2) Output::print("DOF: ", l3, " integral: ", integral, " RHS: ", RHS[l3]);
                }
              }                                   // endif
            }                                     // endif (Cond0==Cond1)
            else
            {
              ErrThrow("different boundary condition on one edge are not allowed!");
            }
          }                                       // endif (boundary joint)
        }                                         // endfor m (N_Joints)
      }                                           // endfor n (rhs)
    }                                             //endif WEAK
  }                                               // endfor i (loop over cells)

  if (out==2)
      Output::print("free memory ");
  if(N_AllMatrices)
  {
    delete LocMatrices;
    delete Matrices[0];
    delete Matrices;
  }

  delete SecondDer;

  delete UsedElements;
  for (i=0;i<4;i++)
  {
    delete X1D[i];
    delete Y1D[i];
    for (unsigned int j=0;j<N_Points1D;j++)
    {
      delete xyval_ref1D[i][j];
      delete xderiv_ref1D[i][j];
      delete yderiv_ref1D[i][j];
    }
  }
  for(l=0;l<ref_n;l++)
  {
    delete value_basefunct_ori[l];
    delete xderiv_basefunct_ori[l];
    delete yderiv_basefunct_ori[l];
    delete value_basefunct_ori_neigh[l];
    delete xderiv_basefunct_ori_neigh[l];
    delete yderiv_basefunct_ori_neigh[l];
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
      for (j=0;j<ref_n;j++)
      {
	  delete [] value_basefunct_ref1D[i][j];
	  delete [] xderiv_basefunct_ref1D[i][j];
	  delete [] yderiv_basefunct_ref1D[i][j];
      }
      delete [] value_basefunct_ref1D[i];
      delete [] xderiv_basefunct_ref1D[i];
      delete [] yderiv_basefunct_ref1D[i];
  }

  for (unsigned int i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete xderiv_refNeigh1D[i];
    delete yderiv_refNeigh1D[i];
  }

  delete aux;
  delete aux2;
  delete aux4;

  if(n_rhs)
  {
    delete righthand;
    delete LocRhs;
  }

  for(i=0; i < N_BaseFuncts2D; i++)
  {
    for(j=0; j < 4; j++)
      {
	for(m=0; m < MaxN_QuadPoints_1D; m++)
          {
	    delete [] xietaval_ref1D[i][j][m];
	    delete [] xideriv_ref1D[i][j][m];
	    delete [] etaderiv_ref1D[i][j][m];
	  }
	delete [] xietaval_ref1D[i][j];
	delete [] xideriv_ref1D[i][j];
	delete [] etaderiv_ref1D[i][j];
      }
    delete [] xietaval_ref1D[i];
    delete [] xideriv_ref1D[i];
    delete [] etaderiv_ref1D[i];
  }

  delete [] value_basefunct_ref1D;
  delete [] xderiv_basefunct_ref1D;
  delete [] yderiv_basefunct_ref1D;

  delete [] xietaval_ref1D;
  delete [] xideriv_ref1D;
  delete [] etaderiv_ref1D;

  if (weak>=1) // ???
  {
    TDatabase::ParamDB->WEAK_BC = 2;
  }
  
  int N_Rows;
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
//   out = 3;
  if(out>2)
  {
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
          cout << "Matrix: " << setw(5) << i << setw(5) << ColInd[j] << "   ";
          cout << setw(10) << Entries[j] << endl;
        }
      }
      cout << endl;
    }                                             // endfor k
  }                                               //endif

  /*  for(k=0;k<n_matrices;k++)
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

}                                                 // end of Assemble

 
#endif


// HIER /////////////////////////////////////////////////////////////////////////////////////////////////////

void Assemble2D(int n_fespaces, const TFESpace2D** fespaces, int n_sqmatrices,
                TSquareMatrix2D** sqmatrices, int n_matrices,
                TMatrix2D** matrices, int n_rhs, double** rhs,
                const TFESpace2D** ferhs, BoundCondFunct2D** BoundaryConditions,
                BoundValueFunct2D** BoundaryValues, LocalAssembling2D& la,
                bool assemble_dirichlet_rows)
{
  if(n_rhs != la.get_n_rhs())
  {
    ErrThrow("the number of right-hand-sides in Assemble2D does not match that "
             "in the LocalAssembling2D object, ", n_rhs, " != ",
             la.get_n_rhs());
  }
  
  std::vector<const TFESpace2D*> further_fe_spaces;
  for(size_t i = 0, n = la.n_fe_functions(); i < n; ++i)
  {
    // space used in local assembling object
    const TFESpace2D* space = la.get_fe_function(i)->GetFESpace2D().get();
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
    double **Matrices = nullptr;
    double ***LocMatrices = nullptr;
    double **LocRhs = nullptr;
    std::vector<const BaseFunctions*> LocBF(n_fespaces + further_fe_spaces.size());

#ifdef __3D__
    double z0, z1;
    ErrThrow("Assemble2D not working in 3D");
#endif

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    double *righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(int i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  int N_AllMatrices = n_sqmatrices+n_matrices;
  if(N_AllMatrices)
  {
    double * aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(int j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(int i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  //SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
  bool *SecondDer = la.GetNeeds2ndDerivatives();
  
  // ########################################################################
  // loop over all cells
  // ########################################################################
  auto Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  int N_Cells = Coll->GetN_Cells();
  for(int i=0;i<N_Cells;i++) // set the cell indices
  {
    TBaseCell *cell = Coll->GetCell(i);
    cell->SetCellIndex(i);
  }
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector<double> AbsDetjk(MaxN_QuadPoints_2D, 1. );

  for(int i=0;i<N_Cells;i++)
  {
    const TBaseCell *cell = Coll->GetCell(i);

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

    bool is_sdfem = (la.get_disctype() == 2); // SDFEM
    if( is_sdfem || (TDatabase::ParamDB->CELL_MEASURE == 4)) // SDFEM
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
      int N_Edges = cell->GetN_Edges();
      for (int ij=0; ij<N_Edges;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if (N_Edges==3)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    
    la.GetLocalForms(qf_orig, LocBF, cell, i, N_AllMatrices, n_rhs, LocMatrices,
                     LocRhs);

    int N_Joints = cell->GetN_Joints();
    // ####################################################################
    // add local/cellwise matrices to global matrices (ansatz == test)
    // ####################################################################
    for(int j=0;j<n_sqmatrices;j++)
    {
      if(sqmatrices[j] == nullptr)
        continue;
      // find space for this bilinear form
      auto fespace = sqmatrices[j]->GetFESpace2D();
      auto element = fespace->get_fe(i);
      int N_ = element.GetN_DOF();

      double **Matrix = LocMatrices[j];
      int DirichletBound = fespace->get_n_active();
      const int * DOF = fespace->GetGlobalDOF(i);

      /*
      BoundaryCondition = BoundaryConditions[j];
      for(m=0;m<N_Joints;m++)
      {
      joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge ||
           joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge)
          {
            boundedge = (TBoundEdge *)joint;
      BoundComp = boundedge->GetBoundComp();
      boundedge->GetParameters(t0, t1);
      }
      else
      {
      isoboundedge = (TIsoBoundEdge *)joint;
      BoundComp = isoboundedge->GetBoundComp();
      isoboundedge->GetParameters(t0, t1);
      }
      // get id of the boundary component
      comp = BoundComp->GetID();
      // get type of the boundary condition at the beginning
      // and at the end of the current edge
      BoundaryCondition(comp, t0, Cond0);
      //      cout << "bound1" << endl;
      if(Cond0 == ROBIN)
      {
      #ifdef __2D__
      //cout << "robin" << endl;
      // Robin
      lr = element.GetBaseFunct2D()->GetPolynomialDegree();

      // get a suitable line quadrature formula
      qf1D = QuadratureFormulaDatabase::qf_from_degree(
          2*lr, BFRefElements::BFUnitLine);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

      element.GetBaseFunct2D()->MakeRefElementData(qf1D.get_type());

      JointValues=FEDatabase::GetJointDerivatives2D(
      element.GetBaseFunct2D_ID(), qf1D.get_type(), m, D00);
      // get vertices of boundary edge
      cell->GetVertex(m)->GetCoords(x0, y0);
      cell->GetVertex((m+1) % 4)->GetCoords(x1, y1);
      // compute (half of the) length of the boundary edge
      hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
      // cout << "x0: " << x0 << " y0: " << y0 << endl;
      // cout << "x1: " << x1 << " y1: " << y1 << endl;
      // compute boundary integral
      for(n=0;n<N_LinePoints;n++)
      {
      // values of test functions in this quadrature point
      JointValue = JointValues[n];
      // cout << "Zeta  :" << zeta[n] << endl;

      // get quadrature point on the boundary
      for(l=0;l<N_;l++)
      {
      MatrixRow = Matrix[l];
      s = JointValue[l];
      // multiply value with weights from quadrature formula
      // and determinant from integral transformation to the
      // unit edge (-1,1)
      s *= hE * LineWeights[n];
      // !! hold the robin boundary values of
      // !! the function alpha from parameters
      // s *= alpha;
      // update rhs for all test functions
      for(k=0;k<N_;k++)
      MatrixRow[k] += s*JointValue[k]*RobinScale;
      } // endfor l
      } // endfor n
      #endif
      } // end Robin
      } // endif BoundEdge
      } // endfor m
      */

      // add local matrix to global
      for(int m=0;m<N_;m++)
      {
        int l=DOF[m];
        double *MatrixRow = Matrix[m];
        if(l < DirichletBound || assemble_dirichlet_rows)
        {
          for(int k=0;k<N_;k++)
          {
            // DOF[k] is the global index of the k-th local degree of freedom
            // MatrixRow[k] is the assembled value corresponding to the m-th
            // local test function and k-th local ansatz function. That means it
            // corresponds to the l=DOF[m]-th global test function and the 
            // DOF[k]-th global ansatz function
             sqmatrices[j]->add(l, DOF[k], MatrixRow[k]);
          }
        }                                         // endif l
        else
        {
          // Dirichlet node
          sqmatrices[j]->set(l,l,1.0);
        }
      }                                           // endfor m
    }                                             // endfor j

    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(int j=0;j<n_matrices;j++)
    {
      if(matrices[j] == nullptr)
        continue;
      auto test_space = matrices[j]->GetTestSpace2D();
      auto ansatz_space = matrices[j]->GetAnsatzSpace2D();
      auto TestElement = test_space->get_fe(i);
      auto AnsatzElement = ansatz_space->get_fe(i);
      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement.GetID() << endl;
      // cout << "AnsatzElement: " << AnsatzElement.GetID() << endl;

      int N_Test = TestElement.GetN_DOF();
      int N_Ansatz = AnsatzElement.GetN_DOF();

      double **Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;

      double *Entries = matrices[j]->GetEntries();
      const int *RowPtr = matrices[j]->get_row_ptr();
      const int *ColInd = matrices[j]->get_vector_columns();

      const int *TestDOF = test_space->GetGlobalDOF(i);
      const int *AnsatzDOF = ansatz_space->GetGlobalDOF(i);

      // add local matrix to global
      for(int m=0;m<N_Test;m++)
      {
        int l=TestDOF[m];
        double *MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        int end=RowPtr[l+1];
        for(int n=RowPtr[l];n<end;n++)
        {
          for(int k=0;k<N_Ansatz;k++)
          {
            if(AnsatzDOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += MatrixRow[k];
              break;
            }                                   // endif
          }                                     // endfor k
        }                                       // endfor n
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)

    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(int j=0;j<n_rhs;j++)
    {
      const TFESpace2D *fespace = ferhs[j];
      int ActiveBound = fespace->get_n_active();
      auto ele = fespace->get_fe(i);
      BoundCond Cond0, Cond1;

      int N_ = ele.GetN_DOF();

      double *local_rhs = LocRhs[j];
      double *RHS = rhs[j];
      if(RHS == nullptr)
      {
        continue;
      }
      // find space for this linear form

      ActiveBound = fespace->get_n_active();

      // dof of the rhs nodes connected to this cell
      const int * DOF = fespace->GetGlobalDOF(i);

      // add local right-hand side to the global one
      for(int m=0;m<N_;m++)
      {
        int l=DOF[m];
        //cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
          // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
        }
      }                                           // endfor m

      BoundCondFunct2D *BoundaryCondition = BoundaryConditions[j];
      BoundValueFunct2D *BoundaryValue = BoundaryValues[j];
      double t0,t1;
      const TBoundComp *BoundComp;
      int N_EdgePoints;
      const double *EdgePoints;
      double eps=1e-4;
        
      //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
      //  continue;

      auto nf = ele.GetNodalFunctional();

      if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
      {
        ErrThrow("What is SUPERCONVERGENCE_ORDER?? most likely not working");
        /* Superconvergence boundary interpolation */
//         if(nf->GetID() == NF_C_Q_Q2_2D)
//           nf = FEDatabase::GetNodalFunctional2D(NF_S_Q_Q2_2D);
      }

      nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      const FEDescriptor *FEDesc_Obj = ele.GetFEDesc();
       int  N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Edges();

      for(int m=0;m<N_Joints;m++)
      {
        auto joint = cell->GetJoint(m);

        if(joint->GetType() == BoundaryEdge ||
           joint->GetType() == IsoBoundEdge ||
           joint->GetType() == InterfaceJoint)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            auto boundedge = (const TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            auto isoboundedge = (const TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          int comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

            // only one boundary condition per edge allowed
            int l,k,l3;
            double s,t;
            int *EdgeDOF;
            double x0, x1, y0,y1;
            double hE;
            double **JointValues;
            double *JointValue;
            double FunctionalValues[MaxN_BaseFunctions2D];
            double PointValues[MaxN_PointsForNodal2D];
            unsigned int N_LinePoints;
            
            if(Cond0 == Cond1)
            {
                switch(Cond0)
                {
              case DIRICHLET:
                // if DG
                if (N_EdgeDOF==0)
                  break;
                // read boundary values for each quadrature point
                for(l=0;l<N_EdgePoints;l++)
                {
                   s = EdgePoints[l];
                   t = 0.5*(t0*(1-s) + t1*(1+s));
                  BoundaryValue(comp, t, PointValues[l]);
                }                                 // endfor l
                // compute boundary values for each dof on the
                // boundary edge with the nodal functionals

                nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                  FunctionalValues);
                EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                // save boundary values of each dof on the boundary
                // edge in the rhs
                for( l=0;l<N_EdgeDOF;l++)
                {
                  RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                }
                break;

              case NEUMANN:
              {
                // get polynomial degree of fe
                 l = ele.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                N_LinePoints = qf1->GetN_QuadPoints();
                JointValues=FEDatabase::GetJointDerivatives2D(
                  *ele.GetBaseFunct(), *qf1, m, MultiIndex2D::D00);
                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute (half of the) length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(unsigned int l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * qf1->get_weight(l);
                  // update rhs for all test functions
                  for( k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                break;
              }
              case ROBIN:
              {
#ifdef __2D__
                // get polynomial degree of fe
                l = ele.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                N_LinePoints = qf1->GetN_QuadPoints();
                JointValues=FEDatabase::GetJointDerivatives2D(
                  *ele.GetBaseFunct(), *qf1, m, MultiIndex2D::D00);
                // get vertices of boundary edge
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                // compute (half of the) length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(unsigned int l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * qf1->get_weight(l);
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
#endif
                break;
              }
              case SLIP:
                ErrThrow("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !");
                break;
              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble2DSlipBC, see below
                break;

            case DIRICHLET_WEAK:
              // do nothing here
              // everything is done in BoundaryAssemble2D
              break;
            case PERIODIC:
              break;
              default :
                ErrThrow("Unknown boundary condition !");

            }                                     // endswitch Cond0
          }                                       // endif (Cond0==Cond1)
          else
          {
            ErrThrow("different boundary condition on one edge are not allowed!");
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)
 
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(int j=0;j<n_sqmatrices;j++)
  {
    if(sqmatrices[j] == nullptr)
      continue;
    //Modifications added to stiffness matrix because of Hanging rows
    if(assemble_dirichlet_rows)
    {
      //Add the hanging rows
      sqmatrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
      //Set the hanging rows to (-0.5,1,-0.5)
      sqmatrices[j]->correct_hanging_rows();
      //Add the hanging columns
      sqmatrices[j]->ModifyMatrixAccordingToCouplingAFC();
    }
    else 
      sqmatrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
  }                                               // endfor j

  for(int j=0;j<n_matrices;j++)
  {
    if(matrices[j] == nullptr)
      continue;
    // note that there used to be code here which did what 
    // ModifyMatrixAccordingToCoupling does but with the ansatz space.
    matrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
  }                                               // endfor j

  for(int j=0;j<n_rhs;j++)
  {
    auto fespace = ferhs[j];

    double *RHS = rhs[j];
    if(RHS == nullptr)
      continue;

    auto hanging_nodes = fespace->get_sorted_hanging_nodes();
    int hanging_bound = fespace->get_n_active();

    // modify the test functions according to coupling, this corresponds to the
    // matrix method ModifyMatrixAccordingToCoupling
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
    // set right-hand-side to zero in the hanging rows
    for(int i = fespace->get_n_active_non_hanging(); i < hanging_bound; ++i)
    {
      RHS[i] = 0;
    }
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  //For assemble_dirichlet_rows the entries are set after AFC is perfomed
  if(!assemble_dirichlet_rows)
  {
    for(int j=0;j<n_sqmatrices;j++)
    {
      if(sqmatrices[j] == nullptr)
        continue;
      sqmatrices[j]->correct_hanging_rows();
    }
  }

  if(n_rhs)
  {
    delete [] LocRhs[0];
    delete [] LocRhs;
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }





  /* std::cout << "\nOUTPUT:" << std::endl; */
  /*   for(int k=0;k<n_sqmatrices;k++) */
  /*   { */
  /*     cout << endl; */
  /*     cout << "sqmatrix: " << k << endl; */
  /*     auto RowPtr = sqmatrices[k]->get_row_ptr(); */
  /*     auto Entries = sqmatrices[k]->GetEntries(); */
  /*     auto ColInd = sqmatrices[k]->get_vector_columns(); */
  /*     auto N_Rows = sqmatrices[k]->get_n_rows(); */
  /*     auto len = std::to_string(N_Rows).length(); */
  /*     for(int i=0;i<N_Rows;i++) */
  /*     { */
  /*       auto end=RowPtr[i+1]; */
  /*       for(int j=RowPtr[i];j<end;j++) */
  /*       { */
  /*         // cout << j << endl; */
  /*         cout << "Matrix[" << setw(len) << i <<','<< setw(len) << ColInd[j] << "] = "<< setw(15) << Entries[j] << endl; */
  /*       } */
  /*       std::cout <<  std::endl; */
  /*     } */
  /*     cout << endl; */
  /*   }                                             // endfor k */

  /*   for(int k=0;k<n_rhs;k++) */
  /*   { */
  /*     cout << "rhs: " << k << endl; */
  /*     auto N_Rows = ferhs[k]->GetN_DegreesOfFreedom(); */
  /*     auto len = std::to_string(N_Rows).length(); */
  /*     auto RHS=rhs[k]; */
  /*     for(int i=0;i<N_Rows;i++) */
  /*       cout << "Rhs["<< setw(len) << i << "] = " << setw(15) << RHS[i] << endl; */
  /*   } */



}                                                 // end of Assemble




// std::vector<parmoon::Point> Transform_Quad_Points( const TBaseCell* cell, const
//     BFRefElements& ref_element, const ReferenceTransformation_type&
//     ref_trans_id, const int& joint_index, TQuadFormula* quad_form )
// {
//   auto n_quad_pts = quad_form->GetN_QuadPoints();
// 
//   // xi, eta and zeta are the d dimensional points on the reference cell.
//   // Those points are computed below by the transformation of the (d-1)
//   // dimensional quadrature points to the d dimensional reference cell.
//   std::vector<double> xi(n_quad_pts);
//   std::vector<double> eta(n_quad_pts);
// 
//   for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
//   {
//     // s is the 1 d coordinate of the quadrature point
//     auto s = quad_form->get_point(quad_pt_i).x;  // 1D coord of point
//     // xi_eta_2D is s transformed to the reference cell, i.e. 2 dimensional
//     auto xi_eta_2D = transform(ref_element, joint_index, s);
//     xi[quad_pt_i] = xi_eta_2D.x;
//     eta[quad_pt_i] = xi_eta_2D.y;
//   } // endfor quad_pt_i
// 
//   // x, y and z are the points on the original cell, i.e. the points after
//   // transformation of xi, eta and zeta
//   std::vector<double> x(n_quad_pts);
//   std::vector<double> y(n_quad_pts);
//   std::vector<double> z(n_quad_pts);
// 
//   FEDatabase::SetCellForRefTrans(cell, ref_trans_id);
//   FEDatabase::GetOrigFromRef(ref_trans_id, n_quad_pts,
//       xi.data(), eta.data(), x.data(), y.data());
// 
//   int d = 2;
//   std::vector<parmoon::Point> transformed_quad_points(n_quad_pts,
//       parmoon::Point((unsigned int) d-1));
//   for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
//   {
//       transformed_quad_points[quad_pt_i] = parmoon::Point( x[quad_pt_i],
//           y[quad_pt_i] );
//   }
// 
//   return transformed_quad_points;
// }


void Assemble2D_JumpStab(int n_fespaces, const TFESpace2D ** fespaces, 
                         int n_sqmatrices, TSquareMatrix2D ** sqmatrices, 
                         int n_matrices, TMatrix2D ** matrices, 
                         int n_rhs, double ** rhs, const TFESpace2D ** ferhs,
                         BoundCondFunct2D ** BoundaryConditions, 
                         BoundValueFunct2D ** BoundaryValues, 
                         LocalAssembling2D& la, 
                         bool assemble_dirichlet_rows)
{
  if(n_rhs != la.get_n_rhs())
  {
    ErrThrow("the number of right-hand-sides in Assemble2D does not match that "
             "in the LocalAssembling2D object, ", n_rhs, " != ",
             la.get_n_rhs());
  }
  
  std::vector<const TFESpace2D*> further_fe_spaces;
  for(size_t i = 0, n = la.n_fe_functions(); i < n; ++i)
  {
    // space used in local assembling object
    const TFESpace2D* space = la.get_fe_function(i)->GetFESpace2D().get();
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
    double **Matrices = nullptr;
    double ***LocMatrices = nullptr;
    double **LocRhs = nullptr;
    std::vector<const BaseFunctions*> LocBF(n_fespaces + further_fe_spaces.size());

#ifdef __3D__
    double z0, z1;
    ErrThrow("Assemble2D not working in 3D");
#endif

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  if(n_rhs)
  {
    LocRhs = new double* [n_rhs];
    double *righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(int i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  int N_AllMatrices = n_sqmatrices+n_matrices;
  if(N_AllMatrices)
  {
    double * aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(int j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(int i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  //SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
  bool *SecondDer = la.GetNeeds2ndDerivatives();
  
  // ########################################################################
  // loop over all cells
  // ########################################################################
  auto Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  int N_Cells = Coll->GetN_Cells();
  int degree, max_degree;
  for(int i=0;i<N_Cells;i++) // set the cell indices
  {
    TBaseCell *cell = Coll->GetCell(i);
    cell->SetCellIndex(i);
    auto element = fespaces[0]->get_fe(i);
    degree = element.GetBaseFunct()->GetPolynomialDegree();
    max_degree = (degree > max_degree) ? degree : max_degree;
  }
  auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                  2*max_degree, BFRefElements::BFUnitLine);
  auto n_quad_pts = qf1->GetN_QuadPoints();
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  std::vector<double> AbsDetjk(MaxN_QuadPoints_2D, 1. );

  for(int i=0;i<N_Cells;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);

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

    bool is_sdfem = (la.get_disctype() == 2); // SDFEM
    if( is_sdfem || (TDatabase::ParamDB->CELL_MEASURE == 4)) // SDFEM
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
      int N_Edges = cell->GetN_Edges();
      for (int ij=0; ij<N_Edges;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if (N_Edges==3)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }
    FEDatabase::GetOrig(LocalUsedElements, Coll, cell, SecondDer, qf_ref,
                        qf_orig);
    
    la.GetLocalForms(qf_orig, LocBF, cell, i, N_AllMatrices, n_rhs, LocMatrices,
                     LocRhs);

    int N_Joints = cell->GetN_Joints();
    cell->SetNormalOrientation();
    
    std::vector<std::vector<double>> val_bas_fct;
    std::vector<std::vector<double>> val_bas_fct_dx;
    std::vector<std::vector<double>> val_bas_fct_dy;  
    
    std::vector<std::vector<double>> val_bas_fct_neigh;
    std::vector<std::vector<double>> val_bas_fct_neigh_dx;
    std::vector<std::vector<double>> val_bas_fct_neigh_dy;
    std::vector<int> index_map;
    for(int joint_i = 0; joint_i < N_Joints; joint_i++)
    {
      val_bas_fct.resize(n_quad_pts);
      val_bas_fct_dx.resize(n_quad_pts);
      val_bas_fct_dy.resize(n_quad_pts);
      // calculate the values of the basis function and 
      // their derivatives 
      auto element = sqmatrices[0]->GetFESpace2D()->get_fe(i);
      auto n_base_func = element.GetN_DOF();      
      
      auto ref_trans_id = element.GetRefTransID();
      auto ref_trans = FEDatabase::GetRefTrans2D(ref_trans_id);
      ref_trans->SetCell(cell);
      
      auto ref_values = FEDatabase::GetJointDerivatives2D(
        *element.GetBaseFunct(), *qf1, joint_i,MultiIndex2D::D00);
      
      auto ref_values_dxi = FEDatabase::GetJointDerivatives2D(
        *element.GetBaseFunct(), *qf1, joint_i, MultiIndex2D::D10);
      
      auto ref_values_deta = FEDatabase::GetJointDerivatives2D(
        *element.GetBaseFunct(), *qf1, joint_i, MultiIndex2D::D01);
      
      // Transform values of basis functions and their derivatives for each
      // quadrature point to ORIGINAL cell
      auto base_vec_dim = element.GetBaseFunct()->GetBaseVectDim();
      
      std::vector<double> u_orig(n_base_func * base_vec_dim);
      std::vector<double> u_orig_dx(n_base_func * base_vec_dim);
      std::vector<double> u_orig_dy(n_base_func * base_vec_dim);
      
      for (int q_pt_i = 0; q_pt_i < n_quad_pts; q_pt_i++)
      {
        auto quad_pt_x = qf1->get_point(q_pt_i).x;
        
        ref_trans->GetOrigValues(joint_i, quad_pt_x, n_base_func,
                                 ref_values[q_pt_i], ref_values_dxi[q_pt_i],
                                 ref_values_deta[q_pt_i], u_orig.data(),
                                 u_orig_dx.data(),
                                 u_orig_dy.data(), base_vec_dim);
      
        element.GetBaseFunct()->ChangeBF(Coll, cell, u_orig.data());
        element.GetBaseFunct()->ChangeBF(Coll, cell, u_orig_dx.data());
        element.GetBaseFunct()->ChangeBF(Coll, cell, u_orig_dy.data());
        
        val_bas_fct[q_pt_i].resize(n_base_func * base_vec_dim);
        val_bas_fct_dx[q_pt_i].resize(n_base_func * base_vec_dim);
        val_bas_fct_dy[q_pt_i].resize(n_base_func * base_vec_dim);
        
        val_bas_fct[q_pt_i].resize(n_base_func * base_vec_dim);
        val_bas_fct_dx[q_pt_i].resize(n_base_func * base_vec_dim);
        val_bas_fct_dy[q_pt_i].resize(n_base_func * base_vec_dim);
        
        for(int j = 0; j < n_base_func; j++)
        {
          val_bas_fct[q_pt_i][j] = u_orig[j];
          val_bas_fct_dx[q_pt_i][j] = u_orig_dx[j];
          val_bas_fct_dy[q_pt_i][j] = u_orig_dy[j];
        }
      }//endfor q_pt_i

      auto ref_element = element.GetBaseFunct()->GetRefElement();
      auto quad_formula_joint = QuadratureFormulaDatabase::qf_from_degree(
      2 * max_degree, BFRefElements::BFUnitLine);
      
      // auto trans_quad_pts = Transform_Quad_Points(
       //  cell, ref_element, ref_trans_id, joint_i, quad_formula_joint);
      /** @------------------------
       *  @----------------- trans_quad_pts
       */
      // xi, eta and zeta are the d dimensional points on the reference cell.
      // Those points are computed below by the transformation of the (d-1)
      // dimensional quadrature points to the d dimensional reference cell.
      std::vector<double> xi(n_quad_pts);
      std::vector<double> eta(n_quad_pts);

      for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
      {
        // s is the 1 d coordinate of the quadrature point
        auto s = qf1->get_point(quad_pt_i).x;  // 1D coord of point
        // xi_eta_2D is s transformed to the reference cell, i.e. 2 dimensional
        auto xi_eta_2D = transform(ref_element, joint_i, s);
        xi[quad_pt_i] = xi_eta_2D.x;
        eta[quad_pt_i] = xi_eta_2D.y;
      } // endfor quad_pt_i

      // x, y and z are the points on the original cell, i.e. the points after
      // transformation of xi, eta and zeta
      std::vector<double> x(n_quad_pts);
      std::vector<double> y(n_quad_pts);
      std::vector<double> z(n_quad_pts);

      FEDatabase::SetCellForRefTrans(cell, ref_trans_id);
      FEDatabase::GetOrigFromRef(ref_trans_id, n_quad_pts,
          xi.data(), eta.data(), x.data(), y.data());

      int d = 2;
      std::vector<parmoon::Point> transformed_quad_points(n_quad_pts,
          parmoon::Point((unsigned int) d-1));
      for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
      {
          transformed_quad_points[quad_pt_i] = parmoon::Point( x[quad_pt_i],
              y[quad_pt_i] );
      }
      /**
       * @---------------endquadpoints-------------
       */
      
      // check for neighbour
      auto joint = cell->GetJoint(joint_i);
      auto neigh = joint->GetNeighbour(cell);
      bool consider_joint = false;
      int cell_nr_neigh = -1; // dummy 
      
      if(neigh)
      {
        if (neigh->GetClipBoard() >= 0)
        {
          cell_nr_neigh = neigh->GetCellIndex();
          
          if(i < cell_nr_neigh)
          {
            consider_joint = true;
            
            // find out the local joint number corresponding to joint_i
            auto facet_nr_neigh = joint->get_joint_nr_in_cell(neigh);
            auto fe_neigh = fespaces[0]->get_fe(cell_nr_neigh);
            
            // compute the basis function and derivatives for the 
            // neighbour element
            //---------------------------------------------------
            val_bas_fct_neigh.resize(n_quad_pts);
            val_bas_fct_neigh_dx.resize(n_quad_pts);
            val_bas_fct_neigh_dy.resize(n_quad_pts);
            // calculate the values of the basis function and 
            // their derivatives
            auto n_base_func = fe_neigh.GetN_DOF();      
            
            auto ref_trans_id = fe_neigh.GetRefTransID();
            auto ref_trans = FEDatabase::GetRefTrans2D(ref_trans_id);
            ref_trans->SetCell(cell);
            
            auto ref_values = FEDatabase::GetJointDerivatives2D(
              *fe_neigh.GetBaseFunct(), *qf1, facet_nr_neigh,MultiIndex2D::D00);
            
            auto ref_values_dxi = FEDatabase::GetJointDerivatives2D(
              *fe_neigh.GetBaseFunct(), *qf1, facet_nr_neigh, MultiIndex2D::D10);
            
            auto ref_values_deta = FEDatabase::GetJointDerivatives2D(
              *fe_neigh.GetBaseFunct(), *qf1, facet_nr_neigh, MultiIndex2D::D01);
            
            // Transform values of basis functions and their derivatives for each
            // quadrature point to ORIGINAL cell
            auto base_vec_dim = fe_neigh.GetBaseFunct()->GetBaseVectDim();
            
            std::vector<double> u_orig(n_base_func * base_vec_dim);
            std::vector<double> u_orig_dx(n_base_func * base_vec_dim);
            std::vector<double> u_orig_dy(n_base_func * base_vec_dim);
            
            for (int q_pt_i = 0; q_pt_i < n_quad_pts; q_pt_i++)
            {
              auto quad_pt_x = qf1->get_point(q_pt_i).x;
              
              ref_trans->GetOrigValues(facet_nr_neigh, quad_pt_x, n_base_func,
                                      ref_values[q_pt_i], ref_values_dxi[q_pt_i],
                                      ref_values_deta[q_pt_i], u_orig.data(),
                                      u_orig_dx.data(),
                                      u_orig_dy.data(), base_vec_dim);
            
              fe_neigh.GetBaseFunct()->ChangeBF(Coll, cell, u_orig.data());
              fe_neigh.GetBaseFunct()->ChangeBF(Coll, cell, u_orig_dx.data());
              fe_neigh.GetBaseFunct()->ChangeBF(Coll, cell, u_orig_dy.data());
              
              val_bas_fct_neigh[q_pt_i].resize(n_base_func * base_vec_dim);
              val_bas_fct_neigh_dx[q_pt_i].resize(n_base_func * base_vec_dim);
              val_bas_fct_neigh_dy[q_pt_i].resize(n_base_func * base_vec_dim);
              
              val_bas_fct_neigh[q_pt_i].resize(n_base_func * base_vec_dim);
              val_bas_fct_neigh_dx[q_pt_i].resize(n_base_func * base_vec_dim);
              val_bas_fct_neigh_dy[q_pt_i].resize(n_base_func * base_vec_dim);
              
              for(int j = 0; j < n_base_func; j++)
              {
                val_bas_fct_neigh[q_pt_i][j] = u_orig[j];
                val_bas_fct_neigh_dx[q_pt_i][j] = u_orig_dx[j];
                val_bas_fct_neigh_dy[q_pt_i][j] = u_orig_dy[j];
              }
            }//endfor q_pt_i
            //---------------------------------------------------
          }//endif i < cell_nr_neigh
        } //endif neigh->GetClipBoard() >= 0
      } //endif neigh
      else // boundary edge
      {
        
      }
      
    }//endfor joint_i<N_Joints
    // ####################################################################
    // add local/cellwise matrices to global matrices (ansatz == test)
    // ####################################################################
    for(int j=0;j<n_sqmatrices;j++)
    {
      if(sqmatrices[j] == nullptr)
        continue;
      // find space for this bilinear form
      auto fespace = sqmatrices[j]->GetFESpace2D();
      auto element = fespace->get_fe(i);
      int N_ = element.GetN_DOF();

      double **Matrix = LocMatrices[j];
      int DirichletBound = fespace->get_n_active();
      const int * DOF = fespace->GetGlobalDOF(i);

      // add local matrix to global
      for(int m=0;m<N_;m++)
      {
        int l=DOF[m];
        double *MatrixRow = Matrix[m];
        if(l < DirichletBound || assemble_dirichlet_rows)
        {
          for(int k=0;k<N_;k++)
          {
            // DOF[k] is the global index of the k-th local degree of freedom
            // MatrixRow[k] is the assembled value corresponding to the m-th
            // local test function and k-th local ansatz function. That means it
            // corresponds to the l=DOF[m]-th global test function and the 
            // DOF[k]-th global ansatz function
             sqmatrices[j]->add(l, DOF[k], MatrixRow[k]);
          }
        }                                         // endif l
        else
        {
          // Dirichlet node
          sqmatrices[j]->set(l,l,1.0);
        }
      }// endfor m
    }// endfor j
    

    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(int j=0;j<n_matrices;j++)
    {
      if(matrices[j] == nullptr)
        continue;
      auto test_space = matrices[j]->GetTestSpace2D();
      auto ansatz_space = matrices[j]->GetAnsatzSpace2D();
      auto TestElement = test_space->get_fe(i);
      auto AnsatzElement = ansatz_space->get_fe(i);
      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement.GetID() << endl;
      // cout << "AnsatzElement: " << AnsatzElement.GetID() << endl;

      int N_Test = TestElement.GetN_DOF();
      int N_Ansatz = AnsatzElement.GetN_DOF();

      double **Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;

      double *Entries = matrices[j]->GetEntries();
      const int *RowPtr = matrices[j]->get_row_ptr();
      const int *ColInd = matrices[j]->get_vector_columns();

      const int *TestDOF = test_space->GetGlobalDOF(i);
      const int *AnsatzDOF = ansatz_space->GetGlobalDOF(i);

      // add local matrix to global
      for(int m=0;m<N_Test;m++)
      {
        int l=TestDOF[m];
        double *MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        int end=RowPtr[l+1];
        for(int n=RowPtr[l];n<end;n++)
        {
          for(int k=0;k<N_Ansatz;k++)
          {
            if(AnsatzDOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += MatrixRow[k];
              break;
            }                                   // endif
          }                                     // endfor k
        }                                       // endfor n
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)

    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(int j=0;j<n_rhs;j++)
    {
      const TFESpace2D *fespace = ferhs[j];
      int ActiveBound = fespace->get_n_active();
      auto ele = fespace->get_fe(i);
      BoundCond Cond0, Cond1;

      int N_ = ele.GetN_DOF();

      double *local_rhs = LocRhs[j];
      double *RHS = rhs[j];
      if(RHS == nullptr)
      {
        continue;
      }
      // find space for this linear form

      ActiveBound = fespace->get_n_active();

      // dof of the rhs nodes connected to this cell
      const int * DOF = fespace->GetGlobalDOF(i);

      // add local right-hand side to the global one
      for(int m=0;m<N_;m++)
      {
        int l=DOF[m];
        //cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
          // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
        }
      }                                           // endfor m

      BoundCondFunct2D *BoundaryCondition = BoundaryConditions[j];
      BoundValueFunct2D *BoundaryValue = BoundaryValues[j];
      double t0,t1;
      const TBoundComp *BoundComp;
      int N_EdgePoints;
      const double *EdgePoints;
      double eps=1e-4;
        
      //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
      //  continue;

      auto nf = ele.GetNodalFunctional();

      if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
      {
        ErrThrow("What is SUPERCONVERGENCE_ORDER?? most likely not working");
        /* Superconvergence boundary interpolation */
//         if(nf->GetID() == NF_C_Q_Q2_2D)
//           nf = FEDatabase::GetNodalFunctional2D(NF_S_Q_Q2_2D);
      }

      nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      const FEDescriptor *FEDesc_Obj = ele.GetFEDesc();
       int  N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Edges();

      for(int m=0;m<N_Joints;m++)
      {
        auto joint = cell->GetJoint(m);

        if(joint->GetType() == BoundaryEdge ||
           joint->GetType() == IsoBoundEdge ||
           joint->GetType() == InterfaceJoint)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            auto boundedge = (const TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            auto isoboundedge = (const TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          int comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

            // only one boundary condition per edge allowed
            int l,k,l3;
            double s,t;
            int *EdgeDOF;
            double x0, x1, y0,y1;
            double hE;
            double **JointValues;
            double *JointValue;
            double FunctionalValues[MaxN_BaseFunctions2D];
            double PointValues[MaxN_PointsForNodal2D];
            unsigned int N_LinePoints;
            
            if(Cond0 == Cond1)
            {
                switch(Cond0)
                {
              case DIRICHLET:
                // if DG
                if (N_EdgeDOF==0)
                  break;
                // read boundary values for each quadrature point
                for(l=0;l<N_EdgePoints;l++)
                {
                   s = EdgePoints[l];
                   t = 0.5*(t0*(1-s) + t1*(1+s));
                  BoundaryValue(comp, t, PointValues[l]);
                }                                 // endfor l
                // compute boundary values for each dof on the
                // boundary edge with the nodal functionals

                nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                  FunctionalValues);
                EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                // save boundary values of each dof on the boundary
                // edge in the rhs
                for( l=0;l<N_EdgeDOF;l++)
                {
                  RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                }
                break;

              case NEUMANN:
              {
                // get polynomial degree of fe
                 l = ele.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                N_LinePoints = qf1->GetN_QuadPoints();
                JointValues=FEDatabase::GetJointDerivatives2D(
                  *ele.GetBaseFunct(), *qf1, m, MultiIndex2D::D00);
                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute (half of the) length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(unsigned int l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * qf1->get_weight(l);
                  // update rhs for all test functions
                  for( k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
                ele.GetBaseFunct()->ChangeBF(Coll, cell, N_LinePoints,
                                               JointValues);
                break;
              }
              case ROBIN:
              {
#ifdef __2D__
                // get polynomial degree of fe
                l = ele.GetBaseFunct()->GetPolynomialDegree();
                // get a suitable line quadrature formula
                auto qf1 = QuadratureFormulaDatabase::qf_from_degree(
                    2*l, BFRefElements::BFUnitLine);
                N_LinePoints = qf1->GetN_QuadPoints();
                JointValues=FEDatabase::GetJointDerivatives2D(
                  *ele.GetBaseFunct(), *qf1, m, MultiIndex2D::D00);
                // get vertices of boundary edge
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                // compute (half of the) length of the boundary edge
                hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(unsigned int l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(qf1->get_point(l).x+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * qf1->get_weight(l);
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
#endif
                break;
              }
              case SLIP:
                ErrThrow("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !");
                break;
              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble2DSlipBC, see below
                break;

            case DIRICHLET_WEAK:
              // do nothing here
              // everything is done in BoundaryAssemble2D
              break;
            case PERIODIC:
              break;
              default :
                ErrThrow("Unknown boundary condition !");

            }                                     // endswitch Cond0
          }                                       // endif (Cond0==Cond1)
          else
          {
            ErrThrow("different boundary condition on one edge are not allowed!");
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)
 
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(int j=0;j<n_sqmatrices;j++)
  {
    if(sqmatrices[j] == nullptr)
      continue;
    //Modifications added to stiffness matrix because of Hanging rows
    if(assemble_dirichlet_rows)
    {
      //Add the hanging rows
      sqmatrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
      //Set the hanging rows to (-0.5,1,-0.5)
      sqmatrices[j]->correct_hanging_rows();
      //Add the hanging columns
      sqmatrices[j]->ModifyMatrixAccordingToCouplingAFC();
    }
    else 
      sqmatrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
  }                                               // endfor j

  for(int j=0;j<n_matrices;j++)
  {
    if(matrices[j] == nullptr)
      continue;
    // note that there used to be code here which did what 
    // ModifyMatrixAccordingToCoupling does but with the ansatz space.
    matrices[j]->ModifyMatrixAccordingToCoupling(assemble_dirichlet_rows);
  }                                               // endfor j

  for(int j=0;j<n_rhs;j++)
  {
    auto fespace = ferhs[j];

    double *RHS = rhs[j];
    if(RHS == nullptr)
      continue;

    auto hanging_nodes = fespace->get_sorted_hanging_nodes();
    int hanging_bound = fespace->get_n_active();

    // modify the test functions according to coupling, this corresponds to the
    // matrix method ModifyMatrixAccordingToCoupling
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
    // set right-hand-side to zero in the hanging rows
    for(int i = fespace->get_n_active_non_hanging(); i < hanging_bound; ++i)
    {
      RHS[i] = 0;
    }
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  //For assemble_dirichlet_rows the entries are set after AFC is perfomed
  if(!assemble_dirichlet_rows)
  {
    for(int j=0;j<n_sqmatrices;j++)
    {
      if(sqmatrices[j] == nullptr)
        continue;
      sqmatrices[j]->correct_hanging_rows();
    }
  }

  if(n_rhs)
  {
    delete [] LocRhs[0];
    delete [] LocRhs;
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }
}                                                 // end of Assemble




