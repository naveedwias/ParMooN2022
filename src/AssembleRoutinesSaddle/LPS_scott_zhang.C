#include "LPS_scott_zhang.h"
#include "Enumerations_fe.h"
#include "FEDatabase.h"
#include "FEMatrix.h"
#include "FEFunction2D.h"
#include "BaseCell.h"

double ComputePSPGParameterLPS(double hK, bool velocity, double nu,
                               size_t lps_coeff_type, double delta0,
                               double delta1)
{
  double val;
  switch(lps_coeff_type)
  {
    // [DGJN17], Section 3
    case 1:
      val = delta0 * hK * hK / nu;
      break;
    // [DGJN17], Section 4, implicit
    case 2:
    // [DGJN17], Section 4, explicit
    case 3:
      if(velocity)
        val = delta1;
      else
        val = delta0 * hK * hK / nu;
      break;
    // [DGJN17], Section 6, implicit
    case 4:
    // [DGJN17], Section 6, explicit
    case 5:
      if(velocity)
        val = delta1 * hK / nu;
      else
        val = delta0 * hK / nu;
      break;
    default:
      val = delta0 * hK * hK / nu;
      break;
  }
  return val;
}

// ======================================================================
// implementation like Badia, CMAME 2012
// ======================================================================
std::shared_ptr<FEMatrix> LPS_for_pressure_Scott_Zhang(
  const std::shared_ptr<const FEMatrix>& C, bool velocity, double,
  const LPS_parameter_set& lps_ps)
{
  double x[4], y[4];
  double val[6];
  MultiIndex2D der[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                          MultiIndex2D::D01};

  constexpr int MaxN_EntriesPerRow = 400;
  constexpr int N_MaxBaseFuncts = 25;
  constexpr int N_MaxQuadPoints = 100;
  double values_in_quad_points[N_MaxBaseFuncts][3][N_MaxQuadPoints];
  double integral_fct_fct[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_grad_grad[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_fct_gradx[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_fct_grady[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double X_dof[N_MaxBaseFuncts], Y_dof[N_MaxBaseFuncts];

  if(!velocity)
  {
    Output::print("LPS for pressure (Scott/Zhang)");
  }
  else
  {
    Output::print("LPS for velocity (Scott/Zhang)");
  }
  auto pressure_space = C->GetFESpace2D();

  // allocate array for node-to-cell map
  int N_P = pressure_space->get_n_dof();
  int N_Active = pressure_space->get_n_active();
  std::vector<int> node_to_cell(N_P, -1);

  // get collection and number of cells
  auto coll = pressure_space->GetCollection();
  int N_Cells = coll->GetN_Cells();

  // Output::print("set node-to-cell map");
  // loop over all cells
  for(int i = 0; i < N_Cells; i++)
  {
    const int *DOF = pressure_space->GetGlobalDOF(i);
    // # basis functions
    int N_ = pressure_space->get_n_local_dof(i);
    if(N_ > N_MaxBaseFuncts)
    {
      ErrThrow("LPS_for_pressure_Scott_Zhang: N_MaxBaseFuncts too small !! ",
               N_);
    }
    for(int j = 0; j < N_; j++)
    {
      int index = DOF[j];
      if(node_to_cell[index] < 0)
        // assign the cell number
        node_to_cell[index] = i;
    }
  }
  // for (int i = 0; i < N_P; i++)
  //   Output::print("ntc ", i, " ", node_to_cell[i]);

  std::vector<int> RowPtr = C->get_row_array();
  // Output::print("original matrix ", C->get_n_entries());

  // compute new sparsity pattern
  // allocate array for storage
  std::vector<int> EntriesPerRow(N_P * MaxN_EntriesPerRow, -1);
  int N_Entries_new = 0;

  for(int i = 0; i < N_Cells; i++)
  {
    const int * DOF = pressure_space->GetGlobalDOF(i);
    // # basis functions
    int N_ = pressure_space->get_n_local_dof(i);
    // loop over the dofs
    for(int j = 0; j < N_; j++)
    {
      int index_j = DOF[j];
      // the corresponding row
      int start = index_j * MaxN_EntriesPerRow;
      // loop over all dofs of same cell
      for(int k = 0; k < N_; k++)
      {
        int index_k = DOF[k];
        // check whether the index index_j*MaxN_EntriesPerRow;is already in the list
        int l;
        for(l = start; l < start + MaxN_EntriesPerRow; l++)
        {
          if(EntriesPerRow[l] == -1)
          {
            EntriesPerRow[l] = index_k;
            N_Entries_new++;
            break;
          }
          if(EntriesPerRow[l] == index_k)
            break;
        }
        if(l == start + MaxN_EntriesPerRow)
        {
          ErrThrow("LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too "
                   "small !!!");
        }
      }

      // take the cell that is assigned to index_j
      int neigh_no_j = node_to_cell[index_j];
      const int *DOF_neigh_j = pressure_space->GetGlobalDOF(neigh_no_j);
      int N_neigh_j = pressure_space->get_n_local_dof(neigh_no_j);
      // loop over the neighbor dofs
      // this loop is just for safety, it should not give new entries
      for(int jj = 0; jj < N_neigh_j; jj++)
      {
        int index_jj = DOF_neigh_j[jj];
        // the indices of this cell might be coupled to index_k
        for(int k = 0; k < N_; k++)
        {
          int index_k = DOF[k];
          int start = index_k * MaxN_EntriesPerRow;
          // check whether the index is already in the list
          int l;
          for(l = start; l < start + MaxN_EntriesPerRow; l++)
          {
            if(EntriesPerRow[l] == -1)
            {
              EntriesPerRow[l] = index_jj;
              N_Entries_new++;
              break;
            }
            if(EntriesPerRow[l] == index_jj)
              break;
          }
          if(l == start + MaxN_EntriesPerRow)
          {
            ErrThrow("LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too "
                     "small !!!");
          }

          // now check the transposed entry
          start = index_jj * MaxN_EntriesPerRow;
          // check whether the index is already in the list
          for(l = start; l < start + MaxN_EntriesPerRow; l++)
          {
            if(EntriesPerRow[l] == -1)
            {
              EntriesPerRow[l] = index_k;
              N_Entries_new++;
              break;
            }
            if(EntriesPerRow[l] == index_k)
              break;
          }
          if(l == start + MaxN_EntriesPerRow)
          {
            ErrThrow("LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too "
                     "small !!!");
          }
        }
      } // end jj

      // start =8*MaxN_EntriesPerRow;
      // for (int l=start;l<start+MaxN_EntriesPerRow;l++)
      //   Output::print(EntriesPerRow[l]);

      // check the far connections
      for(int k = 0; k < N_; k++)
      {
        // take the cell that is assigned to index_k
        int index_k = DOF[k];
        int neigh_no_k = node_to_cell[index_k];
        const int *DOF_neigh_k = pressure_space->GetGlobalDOF(neigh_no_k);
        int N_neigh_k = pressure_space->get_n_local_dof(neigh_no_k);
        // loop over the dofs assigned to index_j and index_k
        for(int jj = 0; jj < N_neigh_j; jj++)
        {
          int index_jj = DOF_neigh_j[jj];
          int start = index_jj * MaxN_EntriesPerRow;
          for(int kk = 0; kk < N_neigh_k; kk++)
          {
            int index_kk = DOF_neigh_k[kk];
            // check whether the index is already in the list
            int l;
            for(l = start; l < start + MaxN_EntriesPerRow; l++)
            {
              if(EntriesPerRow[l] == -1)
              {
                EntriesPerRow[l] = index_kk;
                N_Entries_new++;
                break;
              }
              if(EntriesPerRow[l] == index_kk)
                break;
            }
            if(l == start + MaxN_EntriesPerRow)
            {
              ErrThrow("LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too "
                       "small !!!");
            }
          } // end kk
        } // end jj
      } // end k
    } // end j
  } // end i
  // Output::print("modified matrix ", N_Entries_new);

  ////////////////////////////////////////////
  // new column pointer and entry pointer
  std::vector<int> ColInd_new(N_Entries_new);
  //EntriesC_new = new double[N_Entries_new];
  //memset(EntriesC_new, 0.0, N_Entries_new*SizeOfDouble);
  // reset and fill the pointers
  // loop over all dofs
  int l = 0;
  RowPtr[0] = 0;
  for(int i = 0; i < N_P; i++)
  {
    int start = i * MaxN_EntriesPerRow;
    int k = 0;
    for(int j = start; j < start + MaxN_EntriesPerRow; j++)
    {
      if(EntriesPerRow[j] > -1)
      {
        ColInd_new[l] = EntriesPerRow[j];
        l++;
        k++;
      }
      else
        break;
    }
    RowPtr[i + 1] = RowPtr[i] + k;
  }

  // for (int i = 0; i <= N_P; i++)
  //   Output::print(RowPtr[i] << " ");
  // Output::print(N_Entries_new);
  // for (int i=0;i<=N_Entries_new;i++)
  //   Output::print(ColInd_new[i] << " ");
  // Output::print(l);
  // define new matrix

  auto sqstructureC = std::make_shared<TStructure>(N_P, N_Entries_new,
                                                   ColInd_new.data(),
                                                   RowPtr.data());
  ColInd_new = sqstructureC->get_columns();
  //sqstructureC = new TSquareStructure2D(N_P, N_Entries_new, ColInd_new, RowPtr, NULL);
  auto C_new = std::make_shared<FEMatrix>(pressure_space, sqstructureC);
  double *EntriesC_new = C_new->GetEntries();

  char PString[] = "pressure";
  std::vector<double> tmp_pres(N_P, 0.);
  TFEFunction2D tmp_pressure(pressure_space, PString, tmp_pres.data());
  /////////////////////////////////////////////////
  // compute matrix entries
  std::vector<double> FEFunctValues(MaxN_EntriesPerRow);
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    double h = cell->GetDiameter();
    double delta = ComputePSPGParameterLPS(h, velocity, 1,
                                           lps_ps.lps_coeff_type, lps_ps.delta0,
                                           lps_ps.delta1);
    //delta = compute_PSPG_delta(0.1, h, 0.001);
    const int *DOF = pressure_space->GetGlobalDOF(i);
    auto element = pressure_space->get_fe(i);
    FE_type CurrentElement = pressure_space->get_fe_type(i);
    // # basis functions
    int N_ = pressure_space->get_n_local_dof(i);
    // basis functions

    bool SecondDer[1];
    SecondDer[0] = false;
    std::vector< const FiniteElement*> used_elements(1, &element);
    // get reference transformation, quadrature points
    FEDatabase::GetOrig(used_elements, coll, cell, SecondDer, qf_ref, qf_orig);
    int N_Points = qf_orig.GetN_QuadPoints();
    
    if(N_Points > N_MaxQuadPoints)
    {
      ErrThrow("LPS_for_pressure_Scott_Zhang: N_MaxQuadPoints too small !! ",
               N_Points);
    }

    // assign geometric positions to the dofs
    // compute vertices
    int N_Edges = cell->GetN_Edges();
    for(int j = 0; j < N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
    }

    switch(CurrentElement)
    {
      case C_P1_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = x[1];
        Y_dof[1] = y[1];
        X_dof[2] = x[2];
        Y_dof[2] = y[2];
        break;
      case C_P2_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = (x[0] + x[1]) / 2.0;
        Y_dof[1] = (y[0] + y[1]) / 2.0;
        X_dof[2] = x[1];
        Y_dof[2] = y[1];
        X_dof[3] = (x[0] + x[2]) / 2.0;
        Y_dof[3] = (y[0] + y[2]) / 2.0;
        X_dof[4] = (x[1] + x[2]) / 2.0;
        Y_dof[4] = (y[1] + y[2]) / 2.0;
        X_dof[5] = x[2];
        Y_dof[5] = y[2];
        break;
      case C_P3_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = 2.0 * x[0] / 3.0 + x[1] / 3.0;
        Y_dof[1] = 2.0 * y[0] / 3.0 + y[1] / 3.0;
        X_dof[2] = x[0] / 3.0 + 2.0 * x[1] / 3.0;
        Y_dof[2] = y[0] / 3.0 + 2.0 * y[1] / 3.0;
        X_dof[3] = x[1];
        Y_dof[3] = y[1];
        X_dof[4] = 2.0 * x[0] / 3.0 + x[2] / 3.0;
        Y_dof[4] = 2.0 * y[0] / 3.0 + y[2] / 3.0;
        X_dof[5] = (x[0] + x[1] + x[2]) / 3.0;
        Y_dof[5] = (y[0] + y[1] + y[2]) / 3.0;
        X_dof[6] = 2.0 * x[1] / 3.0 + x[2] / 3.0;
        Y_dof[6] = 2.0 * y[1] / 3.0 + y[2] / 3.0;
        X_dof[7] = x[0] / 3.0 + 2.0 * x[2] / 3.0;
        Y_dof[7] = y[0] / 3.0 + 2.0 * y[2] / 3.0;
        X_dof[8] = x[1] / 3.0 + 2.0 * x[2] / 3.0;
        Y_dof[8] = y[1] / 3.0 + 2.0 * y[2] / 3.0;
        X_dof[9] = x[2];
        Y_dof[9] = y[2];
        break;
      default:
        ErrThrow("LPS_for_pressure_Scott_Zhang: element not implemented!");
    }

    for(int l = 0; l < N_; l++)
      for(int j = 0; j < N_; j++)
        integral_grad_grad[l][j] = integral_fct_fct[l][j] = integral_fct_gradx[l][j] = integral_fct_grady[l][j] = 0.0;

    for(int l = 0; l < N_; l++)
      for(int j = 0; j < 3; j++)
        for(int k = 0; k < N_Points; k++)
          values_in_quad_points[l][j][k] = 0.0;

    // loop over each basis function
    for(int j = 0; j < N_; j++)
    {
      // reset the function values
      for(int k = 0; k < N_; k++)
        FEFunctValues[k] = 0;
      FEFunctValues[j] = 1.0;

      // compute values for all derivatives
      // in all quadrature points
      // in original mesh cell
      for(int k = 0; k < 3; k++)         // for all derivatives
      {
        // get values in original cell
        auto OrigFEValues = FEDatabase::GetOrigElementValues(
          *element.GetBaseFunct(), der[k]);
        for(int jj = 0; jj < N_Points; jj++) // for all quadrature points
        {
          double * Orig = OrigFEValues[jj]; // value in original cell
          double value = 0;
          for(int l = 0; l < N_; l++) // for all basis functions
            // accumulate value of derivative in point j
            value += FEFunctValues[l] * Orig[l];
          values_in_quad_points[j][k][jj] = value;
          //Output::print(i, " ", j, " ", k, " ", value, " :: ");
        }  // endfor jj
      }                                             // endfor k
    }
    // values in the quad points are stored
    // compute integrals

    // loop over each basis function
    for(int j = 0; j < N_; j++) // test fct
    {
      for(int k = 0; k < N_; k++)
      {
        double value = 0.;
        double value_grad_grad = 0.;
        double value_fct_gradx = 0.;
        double value_fct_grady = 0.;
        // loop over quad points
        for(int jj = 0; jj < N_Points; jj++)
        {
          double weight = qf_orig.get_weight(jj);
          value += values_in_quad_points[j][0][jj] * values_in_quad_points[k][0][jj] * weight;
          value_grad_grad += (values_in_quad_points[j][1][jj]
                                *values_in_quad_points[k][1][jj]
                              +values_in_quad_points[j][2][jj]
                                *values_in_quad_points[k][2][jj])
                             * weight;
          value_fct_gradx += values_in_quad_points[j][0][jj]
                            *values_in_quad_points[k][1][jj]
                            *weight;
          value_fct_grady += values_in_quad_points[j][0][jj]
                            *values_in_quad_points[k][2][jj]
                            *weight;
        }
        integral_grad_grad[j][k] += value_grad_grad;
        integral_fct_fct[j][k] += value;
        integral_fct_gradx[j][k] += value_fct_gradx;
        integral_fct_grady[j][k] += value_fct_grady;
      }
    }

    // for (int j = 0; j < N_; j++)
    // {
    //   for (int k = 0; k < N_; k++)
    //     Output::print(j, " ", k, " ", integral_grad_grad[j][k], " :: ");
    // }

    // fill the matrix with standard Brezzi-Pitkaeranta term
    for(int j = 0; j < N_; j++) // test fct.
    {
      int index_test = DOF[j];
//      if (index_test >= N_Active)
//  continue;
      int start = RowPtr[index_test];
      int end = RowPtr[index_test + 1];
      for(int k = 0; k < N_; k++)
      {
        int index_ansatz = DOF[k];
//  if (index_ansatz >=  N_Active)
//    continue;
        for(int l = start; l < end; l++)
        {
          if(index_ansatz == ColInd_new[l])
          {
            EntriesC_new[l] -= delta * integral_grad_grad[j][k];
            break;
          }
        }
      }
    } //  standard Brezzi-Pitkaeranta term done

    // medium far entries
    for(int j = 0; j < N_; j++) // this is c in Badia (2012)
    {
      int index_j = DOF[j];
      // if(index_j >= N_Active)
      //   continue;
      // take the cell that is assigned to index_j
      int neigh_no_j = node_to_cell[index_j];
      auto neigh_j = coll->GetCell(neigh_no_j);
      const int *DOF_neigh_j = pressure_space->GetGlobalDOF(neigh_no_j);
      int N_neigh_j = pressure_space->get_n_local_dof(neigh_no_j);
      // loop over the neighbor dofs
      for(int jj = 0; jj < N_neigh_j; jj++) // this is a in Badia (2012)
      {
        int index_jj = DOF_neigh_j[jj];
        //  if(index_jj >= N_Active)
        //    continue;
        // memset(tmp_pres,0.0,N_P*SizeOfDouble);
        // set the basis fct. in jj
        tmp_pres[index_jj] = 1.0;
        tmp_pressure.FindGradientLocal(neigh_j, neigh_no_j, X_dof[j], Y_dof[j],
                                       val);
        tmp_pres[index_jj] = 0.0;

        // the indices of this cell might be coupled to index_k
        for(int k = 0; k < N_; k++) // this is d in Badia (2012)
        {
          int index_k = DOF[k];
          // if(index_k >= N_Active)
          //   continue;
          double value = delta * (val[1] * integral_fct_gradx[j][k]
                                 +val[2] * integral_fct_grady[j][k]);

          int start = RowPtr[index_k];
          int end = RowPtr[index_k + 1];
          // check whether the index is already in the list
          for(int l = start; l < end; l++)
          {
            if(index_jj == ColInd_new[l])
            {
              EntriesC_new[l] += value;
              break;
            }
          }
          // now do the transposed entry
          start = RowPtr[index_jj];
          end = RowPtr[index_jj + 1];
          // check whether the index is already in the list
          for(int l = start; l < end; l++)
          {
            if(index_k == ColInd_new[l])
            {
              EntriesC_new[l] += value;
              break;
            }
          }
        }
      } // end jj
    } // end j

    // far entries
    for(int j = 0; j < N_; j++) // this is c in Badia (2012)
    {
      int index_j = DOF[j];
      // if(index_j >= N_Active)
      //   continue;
      // take the cell that is assigned to index_j
      int neigh_no_j = node_to_cell[index_j];
      auto neigh_j = coll->GetCell(neigh_no_j);
      const int *DOF_neigh_j = pressure_space->GetGlobalDOF(neigh_no_j);
      int N_neigh_j = pressure_space->get_n_local_dof(neigh_no_j);
      // loop over the neighbor dofs
      for(int jj = 0; jj < N_neigh_j; jj++) // this is a in Badia (2012)
      {
        int index_jj = DOF_neigh_j[jj];
        // if(index_jj >= N_Active)
        //   continue;
        //memset(tmp_pres,0.0,N_P*SizeOfDouble);
        // set the basis fct. in jj
        tmp_pres[index_jj] = 1.0;
        tmp_pressure.FindGradientLocal(neigh_j, neigh_no_j, X_dof[j], Y_dof[j],
                                       val);
        tmp_pres[index_jj] = 0.0;

        // the indices of this cell might be coupled to index_k
        for(int k = 0; k < N_; k++) // this is d in Badia (2012)
        {
          // take the cell that is assigned to index_k
          int index_k = DOF[k];
          // if (index_k >= N_Active)
          // continue;
          int neigh_no_k = node_to_cell[index_k];
          auto neigh_k = coll->GetCell(neigh_no_k);
          const int *DOF_neigh_k = pressure_space->GetGlobalDOF(neigh_no_k);
          int N_neigh_k = pressure_space->get_n_local_dof(neigh_no_k);
          // loop over the neighbor dofs
          for(int kk = 0; kk < N_neigh_k; kk++) // the b in Badia (2012)
          {
            int index_kk = DOF_neigh_k[kk];
            // if(index_kk >= N_Active)
            //   continue;
            //memset(tmp_pres,0.0,N_P*SizeOfDouble);
            //tmp_pres[index_set_in_tmp_pres] = 0.0;
            // set the basis fct. in kk
            tmp_pres[index_kk] = 1.0;
            tmp_pressure.FindGradientLocal(neigh_k, neigh_no_k, X_dof[k],
                                           Y_dof[k], val + 3);
            tmp_pres[index_kk] = 0.0;
            // the value to add
            double value = delta * (val[1] * val[4] + val[2] * val[5])
                          * integral_fct_fct[j][k];
            //if(index_j == 0 && index_k == 1)
            //{
            //  Output::print(neigh_no_j, " ", X_dof[j], " ", Y_dof[j], " ",
            //                neigh_no_k, " ", X_dof[k], " ", Y_dof[k]);
            //  Output::print(val[1], " ", val[2], " ", val[4], " ", val[5],
            //                " ", integral_fct_fct[j][k], " ", value);
            //}
            int start = RowPtr[index_jj];
            int end = RowPtr[index_jj + 1];
            // check whether the index is already in the list
            for(int l = start; l < end; l++)
            {
              if(index_kk == ColInd_new[l])
              {
                EntriesC_new[l] -= value;
                break;
              }
            }
          }
        }
      } // end jj
    } // end j
  } // end i

  // correct Diriclet boundary conditions
  if(velocity && (0))
  {
    Output::print("dof ", N_P, " active ", N_Active);
    int *row_ptr = C_new->get_row_ptr();
    /*for(int i=0;i<N_Active;i++)
    {
       int start = row_ptr[i];
       int end = row_ptr[i+1];
       for(int j=start;j<end;j++)
       {
         if (ColInd_new[j] >= N_Active)
         {
          EntriesC_new[j] = 0.0;
         }
       }
    }*/
    for(int i = N_Active; i < N_P; i++)
    {
      int start = row_ptr[i];
      int end = row_ptr[i + 1];
      for(int j = start; j < end; j++)
        EntriesC_new[j] = 0.0;
    }
  }
  // Output::print("modified matrix ", N_Entries_new)
  // int *row_ptr = C_new->get_row_ptr();
  // for(int i=0;i<C->get_n_rows();i++)
  // //for(int i=0;i<10;i++)
  // {
  //   start = row_ptr[i];
  //   end = row_ptr[i+1];
  //   for(int j = start; j < end; j++)
  //     Output::print(i, " ", ColInd_new[j], " ", EntriesC_new[j], " :: ");
  // }
  return C_new;
}
