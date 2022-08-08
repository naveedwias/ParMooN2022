#include <LinAlg.h>
#ifdef __2D__
#include <FEFunction2D.h>
#elif __3D__
#include <FEFunction3D.h>
#endif
#include <AlgebraicFluxCorrection.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <BaseCell.h>
#include "FEMatrix.h"

#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

#include "HangingNode.h"

#ifdef _MPI
#include <mpi.h>
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

//! Anonymous namespace for helper methods which do not have to be assessed
//! from the outside.
namespace
{
  /** ********************************************************************** */
  /*
   * Modify matrix A according to Eq.(2.6) from BJK17.M3AS, i.e.,
   * a_ji = 0.0 if a_ij < 0.0, i = 1,...,M and j = M+1,...,N
   */
  void modify_system_matrix(FEMatrix& A)
  {
    Output::print<4>("AFC: Modification for BJK limiter");
// #ifdef _MPI
//     ErrThrow("Algebraic flux correction with MPI still needs some work.");
// #endif
    //Column and Row pointer
    const int* ColInd = A.get_vector_columns();
    const int* RowPtr = A.get_row_ptr();
    double* Entries = A.GetEntries();
    int active_nDof = A.get_n_active_rows();
    int index = 0;
    for(int i = 0; i<active_nDof ; i++)
    {
      for(int j = RowPtr[i]; j < RowPtr[i+1]; j++)
      {
        index = ColInd[j];
        //Check whether j = M+1,...,N
        if(index >= active_nDof && Entries[j] < 0.0)
        {
          for(int jj = RowPtr[index]; jj < RowPtr[index+1]; jj++)
          {
            //Transposed Entry a_ji
            if(ColInd[jj] == i)
            {
              //Set a_ji = 0.0
              Entries[jj] = 0.0;
              break;
            }
          }
        }
      }
    }
  }

/** ************************************************************************ */
  /*!
   * Compute the artificial diffusion matrix to a given stiffness matrix.
   *
   * @param[in] A The stiffness matrix which needs additional
   * diffusion. Must be square.
   * @param[out] matrix_D A vector to write the entries
   * of D, the artificial diffusion matrix to. Must contain only 0.
   * MPI: D will leave this method with matrix-consistency level 0, meaning that
   * only its master rows will be correct. Everything else is left at 0.
   *
   */
  void compute_artificial_diffusion_matrix(
    const FEMatrix& A, FEMatrix& D)
  {
    D.reset();
    
    // catch non-square matrix
    if (!A.is_square() )
    {
      ErrThrow("Matrix must be square!");
    }
    if(A.GetStructure()!=D.GetStructure())
    {
      ErrThrow("A and D should have the same structure!");
    }
    // store number of dofs
    int nDof = A.get_n_rows();
    double* D_Entries=D.GetEntries();

#ifdef _MPI
    const TParFECommunicator3D& comm = A.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // get pointers to columns, rows and entries of matrix A
    const int* ColInd = A.get_vector_columns();
    const int* RowPtr = A.get_row_ptr();
    const double* Entries = A.GetEntries();

    // compute off-diagonal entries of the matrix D
    for(int i=0;i<nDof;i++)                       //row loop
    {
#ifdef _MPI
      // do this only for master rows - they couple only to other masters,
      // slaves and halo1 d.o.f., therefore one can be sure, that if A_ij
      // is on this rank, so is the transposed entry, A_ji
      if(masters[i] != rank)
        continue;
#endif
      for(int l=RowPtr[i];l<RowPtr[i+1];l++)
      {
        int j = ColInd[l];                        //column index
        double k_ij = 0;
        double k_ji = 0;
        if (j!=i)
        {                                         // consider only off-diagonals
          k_ij = Entries[l];
          // now get the transposed entry
          for (int ll=RowPtr[j];ll<RowPtr[j+1];ll++)
          {
            if (ColInd[ll]==i)
            {
              k_ji = Entries[ll];
              break;
            }
          }
          //determine entry D_ij         
          D_Entries[l] = std::min({-k_ij , 0.0 , -k_ji});
        }
      }
    }

    // compute diagonal entries of the matrix D
    for(int i=0;i<nDof;i++)                       //row loop
    {
#ifdef _MPI
      // do this only for master rows
      if(masters[i] != rank)
        continue;
#endif
      double val = 0.0;
      // add all entries of i-th row
      int ll = -1;
      for(int l=RowPtr[i];l<RowPtr[i+1];l++)
      {
        val +=  D_Entries[l];
        int j = ColInd[l];
        if (j==i)                                 //diagonal found
          ll = l;                                 //place of the diagonal entry in the matrix_D entries array
      }
      D_Entries[ll] = -val;
    }
  }

  /** Computation of matrix B in modifed Kuzmin
   * @param[in]  system_matrix Matrix A+D
   * @param[in]  alphas Non symmteric matrix alpha
   * @param[in]  stiffness_Entries Matrix A
   * @param[out] B_Entries Matrix B
   */

  void ComputeMatrixB(FEMatrix& system_matrix, 
                      std::vector<double>& alphas,
                      std::vector<double>& stiffness_Entries, 
                      std::vector<double>& B_Entries)
  {
    Output::print<4>("Computation of Matrix B");
    const int* ColInd = system_matrix.get_vector_columns();
    const int* RowPtr = system_matrix.get_row_ptr();
    const int nDofs = system_matrix.get_n_rows();
    double one_alpha_ij = 0.0, one_alpha_ji = 0.0;

    int index = 0, transposed_entry = 0 ;
    for(int i = 0; i<nDofs; i++)
    {
      for(int j  = RowPtr[i]; j<RowPtr[i+1]; j++)
      {
        index = ColInd[j];
        for(int j1 = RowPtr[index]; j1 <RowPtr[index+1]; j1++)
        {
          if(ColInd[j1] == i)
           {
             transposed_entry = j1;
             break;
           }
        }
        one_alpha_ij = (1-alphas[j])*stiffness_Entries[j];
        one_alpha_ji = (1-alphas[transposed_entry])*stiffness_Entries[transposed_entry];
        B_Entries[j] = -std::max({one_alpha_ij, 0.0, one_alpha_ji});
      }
    }

  }


/** ************************************************************************ */  
  /** Computation of the Jacobian times the fluxes for the Zalesak limiter
   *            used in Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] P_plus, Q_plus, P_minus, Q_minus
   *            R_plus, R_minus: The fluxes used for computation of the 
   *            limiters alpha_ij
   * @param[in] index_i: Denotes the ith row of matrix DF for which summation 
   *                     is required
   * @param[in] index_j: Column number for entry a_ij
   * @param[in] entries_pointer_ij: Denotes the pointer to the array Entries
   * @param[in] D: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */
  double Compute_Jacobian_times_flux_Kuzmin(const FEMatrix& A,
    const double * F ,
    const double * P_plus,
    const double * Q_plus,
    const double * Q_minus,
    const double * P_minus,
    const double * R_minus,
    const double * R_plus,
    const int index_i,
    const int index_j,
    const int entries_pointer_ij,
    const FEMatrix& D)
  {
    double sum = 0, epsinv;
    //DQ_minus, DQ_plus, DP_minus, DP_plus: 
    //             Denotes the Derivative of the raw fluxes.
    double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
    const double eps = 1e-14;
    const int* ColInd = A.get_vector_columns();
    const int* RowPtr = A.get_row_ptr();
    const double* Entries = A.GetEntries();
    int row_i_start, row_i_end, temp=0, index_a_ji=0, index_k=0, index_a_ki=0;
    int row_k_start, row_k_end, index_a_kj=0, index_a_jk=0, a_kj=0;
    
    const double* afc_matrix_D_entries=D.GetEntries();

    epsinv = 1.0/eps;
    // computation of the index of entry a_ji
    for(int mm = RowPtr[index_j] ; mm < RowPtr[index_j+1] ; mm++)
    {
      if(ColInd[mm] == index_i)
      {
        index_a_ji = mm;
        break;
      }
    }

    row_i_start = RowPtr[index_i];
    row_i_end = RowPtr[index_i+1];
    // loop over all dofs that are connected with the matrix entry a_{ij}
    for(int k = row_i_start ; k < row_i_end ; k++)
    {
      index_k=ColInd[k];
      row_k_start=RowPtr[index_k];
      row_k_end=RowPtr[index_k+1];
      index_a_kj=-1;
      //to find the index of the entry a_ki
      for(int jj = row_k_start ; jj < row_k_end ; jj++)
      {
        // index of a_kj
        if(ColInd[jj] == index_j)
          index_a_kj = jj;
        // index of transposed entry a_ki
        if(ColInd[jj] == index_i)
          index_a_ki = jj;
      }                                           //end of loop jj
      //case when a_ik>=a_ki
      if (Entries[k] >= Entries[index_a_ki])
      {
        // cases without contribution to the sum
        if(F[k] > 0 && std::abs(R_plus[index_i]-1) < eps)
          continue;
        if(F[k] < 0 && std::abs(R_minus[index_i]-1) < eps)
          continue;
        if(std::abs(F[k]) < eps)
          continue;
        if ((index_i != index_j) && std::abs(afc_matrix_D_entries[entries_pointer_ij])<1e-14)
          continue;
        // *** first case with contribution ***
        if(F[k] > 0)
        {
          //Calculations for finding derivative of Q_plus
          if (index_i != index_j)
          {
            if(F[entries_pointer_ij] >= 0)
              DQ_plus = 0;
            else
              DQ_plus = -afc_matrix_D_entries[entries_pointer_ij];
          }
          else
          {
            sumD = 0;
            for(int m = row_i_start ; m < row_i_end ; m++)
            {
              if(F[m] < 0)
                sumD += afc_matrix_D_entries[m];
            }
            DQ_plus = sumD;
          }

          // calculations for finding derivative of P_plus
          if (index_i != index_j)
          {
            if(F[entries_pointer_ij] <= 0)
            {
              DP_plus = 0;
            }
            else
            {
              if (Entries[entries_pointer_ij] >= Entries[index_a_ji])
                DP_plus = afc_matrix_D_entries[entries_pointer_ij];
              else
                DP_plus = 0;
            }
          }
          else
          {
            sumD = 0;
            for(int m = row_i_start ; m < row_i_end ; m++)
            {
              int transposed_entry = ColInd[m];
              for(int ll = RowPtr[transposed_entry] ; ll < RowPtr[transposed_entry+1];
                  ll++)
              {
                if(ColInd[ll] == index_i)
                {
                  temp = ll;
                  break;
                }
              }
              if(F[m] > 0 && Entries[m] >= Entries[temp])
                sumD += afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_plus = -sumD;
          }
          double val =  DQ_plus/P_plus[index_i]-Q_plus[index_i]*
                                DP_plus/(P_plus[index_i]*P_plus[index_i]);
          if (std::abs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        // *** second case with contribution ***
        if(F[k]<0)
        {
          // calculations for finding derivative of Q_minus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]<=0)
              DQ_minus=0;
            else
              DQ_minus=-afc_matrix_D_entries[entries_pointer_ij];
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_minus=sumD;
          }
          //Calculations for finding derivative of P_minus
          if (index_i!=index_j)
          {
            if(F[entries_pointer_ij]>=0)
              DP_minus=0;
            else
            {
              if(Entries[entries_pointer_ij]>=Entries[index_a_ji])
                DP_minus=afc_matrix_D_entries[entries_pointer_ij];
              else
                DP_minus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_i_start;m<row_i_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                  ll++)
              {
                if(ColInd[ll]==index_i)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]<0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_minus=-sumD;
          }
          double val = DQ_minus/P_minus[index_i]-Q_minus[index_i]*
                          DP_minus/(P_minus[index_i]*P_minus[index_i]);
          if (std::abs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        ErrThrow("Compute_Jacobian_times_flux_Kuzmin: index pair without case");
      }                                        //end of if case when a_ik>=a_ki
      //case when a_ki>a_ik
      else
      {
        // cases without contribution to the sum
        if(F[k]>0 && std::abs(R_minus[index_k]-1)<eps)
          continue;
        if(F[k]<0 && std::abs(R_plus[index_k]-1)<eps)
          continue;
        if(std::abs(F[k])<eps)
          continue;
        double f_kj, d_kj;
        if(index_a_kj != -1)
        {
          f_kj=F[index_a_kj];
          d_kj=afc_matrix_D_entries[index_a_kj];
          a_kj=Entries[index_a_kj];
        }
        else
        {
          // next k
          continue;
        }     
        if ((index_k != index_j)&& (std::abs(d_kj) < 1e-14))
          continue;
        // computation of the index of entry a_jk
        for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
        {
          if(ColInd[mm]==index_k)
          {
            index_a_jk=mm;
            break;
          }
        }
        // *** first case with contribution ***
        if(F[k]>0)
        {
          //Calculations for finding derivative of Q_minus
          if (index_k!=index_j)
          {
            if(f_kj<=0)
              DQ_minus=0;
            else
              DQ_minus=-d_kj;
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_minus=sumD;
          }
          // calculations for finding derivative of P_minus
          if (index_k!=index_j)
          {
            if(f_kj>=0)
            {
              DP_minus=0;
            }
            else
            {
              if (a_kj>=Entries[index_a_jk])
                DP_minus=d_kj;
              else
                DP_minus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                   ll++)
              {
                if(ColInd[ll]==index_k)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]<0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_minus=-sumD;
          }
          double val =  DQ_minus/P_minus[index_k]-Q_minus[index_k]*
                                 DP_minus/(P_minus[index_k]*P_minus[index_k]);
          if (std::abs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        // *** second case with contribution ***
        if(F[k]<0)
        {
          // calculations for finding derivative of Q_plus
          if (index_k!=index_j)
          {
            if(f_kj>=0)
              DQ_plus=0;
            else
              DQ_plus=-d_kj;
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DQ_plus=sumD;
            //Output::print<2>(DQ_plus);
          }
          //Calculations for finding derivative of P_plus
          if (index_k!=index_j)
          {
            if(f_kj<=0)
              DP_plus=0;
            else
            {
              if(a_kj>=Entries[index_a_jk])
                DP_plus=d_kj;
              else
                DP_plus=0;
            }
          }
          else
          {
            sumD=0;
            for(int m=row_k_start;m<row_k_end;m++)
            {
              int transposed_entry=ColInd[m];
              for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];
                   ll++)
              {
                if(ColInd[ll]==index_k)
                {
                  temp=ll;
                  break;
                }
              }
              if(F[m]>0 && Entries[m]>=Entries[temp])
                sumD+=afc_matrix_D_entries[m];
            }                                     //End of loop m
            DP_plus=-sumD;
          }
          double val = DQ_plus/P_plus[index_k]-Q_plus[index_k]*
                                   DP_plus/(P_plus[index_k]*P_plus[index_k]);
          if (std::abs(val) < epsinv)
            sum += val*F[k];
          continue;
        }
        ErrThrow("Compute_Jacobian_times_flux_Kuzmin: index pair without case");
      }                                      //end of else case when a_ki>a_ik
    }
    return sum;
  }
  
/** ************************************************************************ */  

  /* regularized maximum function given as 
   * max{x,y}=(x+y+std::sqrt((x+y)^2+sigma))/2*/
  double smooth_max(const double x, const double y, const double sigma)
  {
    double value;
    value=(x+y+std::sqrt((x-y)*(x-y)+sigma))*0.5;
    return value;
  }

  /* regularixed minimum function given as 
   * min{x,y}=(x+y-std::sqrt((x-y)^2+sigma))/2*/
  double smooth_min(const double x, const double y, const double sigma)
  {
    double value;
    value=(x+y-std::sqrt((x-y)*(x-y)+sigma))*0.5;
    return value;
  }

  /* derivative of regularized minimum function*/
  double der_smooth_min(const double x, const double y, const double sigma)
  {
    double value;
    value=1+(x-y)/std::sqrt((x-y)*(x-y)+sigma);
    return value;
  }

/** ************************************************************************ */  
  /** Computation of the Jacobian times the fluxes for the regularized Kuzmin 
   *            limiter used in Regularized Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] index_i: Denotes the ith row of matrix DF for which summation 
   *                     is required
   * @param[in] index_j: Column number for entry a_ij
   * @param[in] entries_pointer_ij: Denotes the pointer to the array Entries
   * @param[in] afc_matrix_D_entries: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */
  double Compute_Jacobian_times_flux_Kuzmin_Regularized(const FEMatrix& A,
    const double * F ,
    const int index_i,
    const int index_j,
    const int entries_pointer_ij,
    const FEMatrix& D, 
    const double sigma, 
    const double omega_matrix, const double omega_derivative)
  {
    //Regularized P_plus, P_minus, Q_plus, Q_minus
    double P_plus_i_reg=0.0;
    double P_minus_i_reg=0.0;
    double Q_plus_i_reg=0.0;
    double Q_minus_i_reg=0.0;
    double P_plus_k_reg, P_minus_k_reg, Q_plus_k_reg, Q_minus_k_reg;
    double a_kj, a_jk = 0, a_ji = 0;
    double sum=0;
    
    //DQ_minus, DQ_plus, DP_minus, DP_plus: Denotes the Derivative of the raw fluxes.
    double DQ_plus_i, DQ_minus_i, DP_plus_i,DP_minus_i, DF_ik, 
           DR_plus_i, DR_minus_i;
    double DQ_plus_k, DQ_minus_k, DP_plus_k,DP_minus_k, 
           DR_plus_k, DR_minus_k;
    double R_plus_i, R_minus_i, R_plus_k, R_minus_k;

    const int* ColInd = A.get_vector_columns();
    const int* RowPtr = A.get_row_ptr();
    const double* Entries = A.GetEntries();
    int row_i_start, row_i_end, temp = 0, index_a_ji = 0;
    
    // indicies used in the loop
    int k, m, index_m, m2, m3, mm, row_k_start, row_k_end, index_k;
    int index_a_ki = -1, index_a_kj, index_a_jk = 0;
    
    const double* afc_matrix_D_entries=D.GetEntries();

    // i-th row of sqmatrix
    row_i_start = RowPtr[index_i];
    row_i_end = RowPtr[index_i+1];

    // computation of the index of entry a_ji
    for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
    {
      if(ColInd[mm]==index_i)
      {
        index_a_ji=mm;
        a_ji=Entries[index_a_ji];
        break;
      }
    }

    // computation of regularized Q_plus_i, Q_minus_i, P_plus_i, P_minus_i
    for(m=row_i_start;m<row_i_end;m++)
    {
      Q_minus_i_reg -= smooth_max(0,F[m],sigma);
      Q_plus_i_reg -= smooth_min(0,F[m],sigma);
      index_m=ColInd[m];
      m2=RowPtr[index_m];
      m3=RowPtr[index_m+1];
      //finding index of transposed entry a_im used in P_plus_i and P_minus_i
      for(mm=m2;mm<m3;mm++)
      {
        if(ColInd[mm]==index_i)
          break;
      }                                           //end of loop mm
      if(Entries[mm]>Entries[m])
        continue;
      P_plus_i_reg += smooth_max(0,F[m],sigma);
      P_minus_i_reg += smooth_min(0,F[m],sigma);
    }                                             //end of loop m
    R_plus_i = smooth_min(1,Q_plus_i_reg/P_plus_i_reg,sigma)/2.0;
    R_minus_i = smooth_min(1,Q_minus_i_reg/P_minus_i_reg,sigma)/2.0;

    // computation of derivatives in node i
    if (index_i!=index_j)
    {
      double d_ij = afc_matrix_D_entries[entries_pointer_ij];
      double f_ij = F[entries_pointer_ij];
      double f_ij_f_ij2_sigma = f_ij/std::sqrt(f_ij * f_ij + sigma);

      DQ_plus_i = -(1-f_ij_f_ij2_sigma)*d_ij/2.0;
      DQ_minus_i = -(1+f_ij_f_ij2_sigma)*d_ij/2.0;
      if (Entries[entries_pointer_ij]>=a_ji)
      {
        DP_plus_i = (1+f_ij_f_ij2_sigma)*d_ij/2.0;
        DP_minus_i = (1-f_ij_f_ij2_sigma)*d_ij/2.0;
      }
      else
      {
        DP_plus_i = DP_minus_i = 0;
      }
    }
    else                                          // index_i = index_j
    {
      // DQ_plus, DQ_minus
      DQ_plus_i=DQ_minus_i=0;
      for(int m=row_i_start;m<row_i_end;m++)
      {
        //the case k!=i
        if(ColInd[m]==index_i)
          continue;
        DQ_plus_i += (1-F[m]/std::sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
        DQ_minus_i += (1+F[m]/std::sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
      }
      DQ_plus_i /= 2.0;
      DQ_minus_i /= 2.0;
      // DP_plus, DP_minus
      DP_plus_i=DP_minus_i=0;
      for(int m=row_i_start;m<row_i_end;m++)
      {
        //the case m!=i
        if(ColInd[m]==index_i)
          continue;
        int transposed_entry=ColInd[m];
        for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
        {
          if(ColInd[ll]==index_i)
          {
            temp=ll;
            break;
          }
        }
        if(Entries[m]>=Entries[temp])
        {
          DP_plus_i -= (1+F[m]/std::sqrt(F[m]*F[m]+sigma)) * afc_matrix_D_entries[m];
          DP_minus_i -= (1-F[m]/std::sqrt(F[m]*F[m]+sigma)) * afc_matrix_D_entries[m];
        }
      }                                           //End of loop m
      DP_plus_i /= 2.0;
      DP_minus_i /= 2.0;
    }
    DR_plus_i = ((DQ_plus_i*P_plus_i_reg-Q_plus_i_reg*DP_plus_i)/(P_plus_i_reg*P_plus_i_reg))
      *der_smooth_min(1,Q_plus_i_reg/(P_plus_i_reg),sigma)/2.0;
    DR_minus_i = ((DQ_minus_i*P_minus_i_reg-Q_minus_i_reg*DP_minus_i)/(P_minus_i_reg*P_minus_i_reg))
      *der_smooth_min(1,Q_minus_i_reg/(P_minus_i_reg),sigma)/2.0;

    //loop to find sum over beta_ik
    for(k=row_i_start;k<row_i_end;k++)
    {
      index_k=ColInd[k];
      row_k_start=RowPtr[index_k];
      row_k_end=RowPtr[index_k+1];
      index_a_kj=-1;
      //to find the index of the entry a_ki
      for(int jj=row_k_start;jj<row_k_end;jj++)
      {
        // index of a_kj
        if(ColInd[jj]==index_j)
          index_a_kj=jj;
        // index of transposed entry a_ki
        if(ColInd[jj]==index_i)
          index_a_ki=jj;
      }                                           //end of loop jj

      //case when a_ik>=a_ki
      if (Entries[k]>=Entries[index_a_ki])
      {

        // derivative of F_ik
        if (index_i==index_j)
        {
          if(index_k==index_j)
            DF_ik=0;
          else
            DF_ik=-afc_matrix_D_entries[k];
        }
        else
        {
          if(index_k==index_j)
            DF_ik=afc_matrix_D_entries[entries_pointer_ij];
          else
            DF_ik=0;
        }
        sum += omega_derivative * DR_plus_i*smooth_max(0,F[k],sigma) 
              + omega_matrix * R_plus_i*der_smooth_min(F[k],0,sigma)*DF_ik
              + omega_derivative * DR_minus_i*smooth_min(0,F[k],sigma) 
              + omega_matrix * R_minus_i*der_smooth_min(0,F[k],sigma)*DF_ik;
      }
      else                                        //case when a_ki>a_ik
      {
        double f_kj, d_kj;
        if(index_a_kj != -1)
        {
          f_kj=F[index_a_kj];
          d_kj=afc_matrix_D_entries[index_a_kj];
          a_kj=Entries[index_a_kj];
        }
        else
        {
          // next k
          continue;
        }
        // computation of the index of entry a_jk
        for(int mm=RowPtr[index_j];mm<RowPtr[index_j+1];mm++)
        {
          if(ColInd[mm]==index_k)
          {
            index_a_jk = mm;
            a_jk = Entries[index_a_jk];
            break;
          }
        }

        // computation of regularized Q_plus_k, Q_minus_k, P_plus_k, P_minus_k
        P_plus_k_reg = P_minus_k_reg = Q_plus_k_reg =  Q_minus_k_reg = 0.0;
        for(m=row_k_start;m<row_k_end;m++)
        {
          Q_minus_k_reg -= smooth_max(0,F[m],sigma);
          Q_plus_k_reg -= smooth_min(0,F[m],sigma);
          index_m=ColInd[m];
          m2=RowPtr[index_m];
          m3=RowPtr[index_m+1];
          //finding index of transposed entry a_km used in P_plus_k and P_minus_k
          for(mm=m2;mm<m3;mm++)
          {
            if(ColInd[mm]==index_k)
              break;
          }                                       //end of loop mm
          if(Entries[mm]>Entries[m])
            continue;
          P_plus_k_reg += smooth_max(0,F[m],sigma);
          P_minus_k_reg += smooth_min(0,F[m],sigma);
        }
        R_plus_k = smooth_min(1,Q_plus_k_reg/P_plus_k_reg,sigma)/2.0;
        R_minus_k = smooth_min(1,Q_minus_k_reg/P_minus_k_reg,sigma)/2.0;
        
        // computation of derivatives in node k
        if (index_k!=index_j)
        {
          //double d_kj = afc_matrix_D_entries[entries_pointer_ij];
          //double f_kj = F[entries_pointer_ij];
     
          double f_kj_f_kj2_sigma = f_kj/std::sqrt(f_kj * f_kj + sigma);
          DQ_plus_k = -(1-f_kj_f_kj2_sigma)*d_kj/2.0;
          DQ_minus_k = -(1+f_kj_f_kj2_sigma)*d_kj/2.0;
          
          if (a_kj >= a_jk)
          {
            DP_plus_k = (1+f_kj_f_kj2_sigma)*d_kj/2.0;
            DP_minus_k = (1-f_kj_f_kj2_sigma)*d_kj/2.0;
          }
          else
          {
            DP_plus_k = DP_minus_k = 0;
          }
        }
        else                                          // index_k = index_j
        {
          // DQ_plus, DQ_minus
          DQ_plus_k=DQ_minus_k=0;
          for(int m=row_k_start;m<row_k_end;m++)
          {
            //the case k!=m
            if (ColInd[m]==index_k)
              continue;
            DQ_plus_k += (1-F[m]/std::sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
            DQ_minus_k += (1+F[m]/std::sqrt(F[m]*F[m]+sigma))*afc_matrix_D_entries[m];
          }
          DQ_plus_k /= 2.0;
          DQ_minus_k /= 2.0;
          // DP_plus, DP_minus
          DP_plus_k=DP_minus_k=0;
          for(int m=row_k_start;m<row_k_end;m++)
          {
            //the case m!=k
            if(ColInd[m]==index_k)
              continue;
            int transposed_entry=ColInd[m];
            for(int ll=RowPtr[transposed_entry];ll<RowPtr[transposed_entry+1];ll++)
            {
              if(ColInd[ll]==index_k)
              {
                temp=ll;
                break;
              }
            }
            if(Entries[m]>=Entries[temp])
            {
              DP_plus_k -= (1+F[m]/std::sqrt(F[m]*F[m]+sigma))
                          *afc_matrix_D_entries[m];
              DP_minus_k -= (1-F[m]/std::sqrt(F[m]*F[m]+sigma))
                          *afc_matrix_D_entries[m];
            }
          }                                           //End of loop m
          DP_plus_k /= 2.0;
          DP_minus_k /= 2.0;
        }
        DR_plus_k = ((DQ_plus_k*P_plus_k_reg-Q_plus_k_reg*DP_plus_k)/(P_plus_k_reg*P_plus_k_reg))
          *der_smooth_min(1,Q_plus_k_reg/P_plus_k_reg,sigma)/2.0;
        DR_minus_k = ((DQ_minus_k*P_minus_k_reg-Q_minus_k_reg*DP_minus_k)/(P_minus_k_reg*P_minus_k_reg))
          *der_smooth_min(1,Q_minus_k_reg/P_minus_k_reg,sigma)/2.0;

        // derivative of F_ik
        if (index_k==index_j)
        {
          if(index_i==index_j)
            DF_ik=0;
          else
            DF_ik=-afc_matrix_D_entries[k];
        }
        else
        {
          if(index_i==index_j)
            DF_ik= afc_matrix_D_entries[k];
          else
            DF_ik=0;
        }
        sum += -omega_derivative * DR_plus_k*smooth_max(0,-F[k],sigma) 
               - omega_matrix * R_plus_k*der_smooth_min(-F[k],0,sigma)*DF_ik
               -omega_derivative * DR_minus_k*smooth_min(0,-F[k],sigma) 
               - omega_matrix * R_minus_k*der_smooth_min(0,-F[k],sigma)*DF_ik;
      }
    }
    return sum;
  }

/** ************************************************************************ */  
  
  /** Computation of the Jacobian times the fluxes for the BJK17 limiter
   *            used in regularixed Newton's method
   * @param[in] A: The system matrix A
   * @param[in] F: The Matrix where entries f_ij=d_ij-(uj_-u_i)
   * @param[in] P_plus, Q_plus, P_minus, Q_minus, R_plus, R_minus: 
   *            The fluxes used for computation of the limiters alpha_ij
   * @param[in] umax: Vector where umax_i=max{u_ij, 1<=j<=Ndofs}
   * @param[in] umin: Vector where umin_i=min{u_ij, 1<=j<=Ndofs}
   * @param[in] q: Vector where q_i=gamma_i\sum d_ij
   * @param[in] sol_col_j: entry of sol[col_j]
   * @param[in] row_i: Denotes the ith row of matrix DF for which summation 
   *                   is required
   * @param[in] col_j: Column number for entry a_ij
   * @param[in] entries_pointer: Denotes the pointer to the array Entries
   * @param[in] D: The artificial diffusion matrix
   *
   * @param[out] The summation required in the entry of the Jacobian matrix.
   */

  double Compute_Jacobian_times_flux_BJK17(const FEMatrix& A,
    const double * F ,
    const double * P_plus,
    const double * P_minus,
    const double * Q_plus,
    const double * Q_minus,
    const double * R_plus,
    const double * R_minus,
    const std::vector<double>& umax,
    const std::vector<double>& umin,
    const std::vector<double>& q,
    const double sol_col_j,
    const int row_i,
    const int col_j,
    const int entries_pointer,
    const FEMatrix& D)

  {
    double sum=0, epsinv;
    double DQ_plus, DQ_minus,DP_plus,DP_minus, sumD=0;
    const double eps = 1e-14;
    const int* ColInd = A.get_vector_columns();
    const int* RowPtr = A.get_row_ptr();
    int k0, k1, m0, m1;

    epsinv = 1.0/eps;

    k0=RowPtr[row_i];
    k1=RowPtr[row_i+1];
    m0=RowPtr[col_j];
    m1=RowPtr[col_j+1];
    const double* afc_matrix_D_entries=D.GetEntries();

    for(int k=k0;k<k1;k++)
    {
      // situations that do not lead to contributions
      if(F[k]>0 && std::abs(R_plus[row_i]-1.0) < eps)
        continue;
      if(F[k]<0 && std::abs(R_minus[row_i]-1.0) < eps)
        continue;
      if(std::abs(F[k]) < eps)
        continue;

      int col_k = ColInd[k];
      // *** first condition with a contribution ***
      // if(F[k]>0 && R_plus[row_i]<1)
      if(F[k]>0)
      {
        if(R_plus[row_i]<=R_minus[col_k])
        {
          // computation of derivatives of Q_plus and P_plus
          // NOTE: the index corresponding to umax need not to be uniquely 
          //       defined
          DQ_plus = DP_plus = 0.0;
          if(row_i!=col_j)
          {
            if(std::abs(umax[row_i]-sol_col_j) < eps)
              DQ_plus=-q[row_i];
            if(F[entries_pointer]>0)
              DP_plus=afc_matrix_D_entries[entries_pointer];
          }
          else
          {
            if(std::abs(umax[row_i]-sol_col_j) >= eps)
              DQ_plus=q[row_i];
            sumD=0.0;
            for(int m=k0;m<k1;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_plus=-sumD;
          }
          if (P_plus[row_i] > 0.0)
          {
            double val = DQ_plus/P_plus[row_i]
                       - Q_plus[row_i]*DP_plus/(P_plus[row_i]*P_plus[row_i]);
            if (std::abs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
        else
        {
          // computation of derivatives of Q_minus and P_minus
          DQ_minus = DP_minus = 0.0;
          if(col_k!=col_j)
          {
            if(std::abs(umin[col_k]-sol_col_j) < eps)
              DQ_minus=-q[col_k];
            if(F[entries_pointer]<0)
              DP_minus=afc_matrix_D_entries[entries_pointer];
          }
          // col_k = col_j
          else
          {
            if(std::abs(umin[col_j]-sol_col_j) >= eps)
              DQ_minus=q[col_j];
            sumD=0;
            for(int m=m0;m<m1;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_minus=-sumD;
          }
          if (P_minus[col_k] < 0.0)
          {
            double val = DQ_minus/P_minus[col_k] 
                       - Q_minus[col_k]*DP_minus/(P_minus[col_k]*P_minus[col_k]);
            if (std::abs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
      }
      // *** second condition with a contribution ***
      //if(F[k]<0 && R_minus[row_i]<1)
      if (F[k] < 0)
      {
        if(R_minus[row_i]<=R_plus[col_k])
        {
          DQ_minus = DP_minus = 0.0;
          // computation of derivatives of Q_minus and P_minus
          if(row_i!=col_j)
          {
            if(std::abs(umin[row_i]-sol_col_j) < eps)
              DQ_minus=-q[row_i];
            if(F[entries_pointer]<0)
              DP_minus=afc_matrix_D_entries[entries_pointer];
          }
          else
          {
            if(std::abs(umin[row_i]-sol_col_j) >= eps)
              DQ_minus=q[row_i];
            sumD=0;
            for(int m=k0;m<k1;m++)
            {
              if(F[m]<0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_minus=-sumD;
          }
          if (P_minus[row_i] < 0.0)
          {
            double val = DQ_minus/P_minus[row_i]
                       - Q_minus[row_i]*DP_minus/(P_minus[row_i]*P_minus[row_i]);
            if (std::abs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
        // R_minus[row_i] > R_plus[ColInd[k]]
        else
        {
          DQ_plus = DP_plus = 0.0;
          // computation of derivatives of Q_plus and P_plus
          if(col_k!=col_j)
          {
            if(std::abs(umax[col_k]-sol_col_j)<eps)
              DQ_plus=-q[col_k];
            if(F[entries_pointer]>0)
              DP_plus=afc_matrix_D_entries[entries_pointer];
          }
          // col_k = col_j
          else
          {
            if(std::abs(umax[col_j]-sol_col_j) >= eps)
              DQ_plus=q[col_j];
            sumD=0;
            for(int m=m0;m<m1;m++)
            {
              if(F[m]>0)
                sumD+=afc_matrix_D_entries[m];
            }
            DP_plus=-sumD;
          }
          if (P_plus[col_k] > 0.0)
          {
            double val = DQ_plus/P_plus[col_k] 
                       - Q_plus[col_k]*DP_plus/(P_plus[col_k]*P_plus[col_k]);
            if (std::abs(val) < epsinv)
              sum+=val*F[k];
          }
          continue;
        }
      }
      ErrThrow("Compute_Jacobian_times_flux_BJK17: index pair without case",
               k," ",F[k]);
    }
    return sum;
  }

//Check the sign of point (x0,y0) from the line by points (x1,y1) and (x2, y2)
//This is used to compute the convex hull for BJK limiter in the case for grids with 
//hanging nodes
  int check_sign(const double x0, const double y0, 
                 const double x1, const double y1,
                 const double x2, const double y2)
  {
    double sign_check = 0.0, slope = 0.0;
    slope = (y2 - y1)/(x2 - x1);
    sign_check = (y0 - y1) - slope*(x0 - x1);
    if(sign_check > 0)
      return 1;
    else if(sign_check == 0)
      return 0;
    else
      return -1;
  }

/** ************************************************************************ */  
  
  /**
   * Compute the weights for the linearty preserving limiter from 
   * Barrenechea, John, Knobloch; M3AS 2017. This implementation
   * can be used for grids with hanging nodes
   *
   * @param[in]  FEMatrix system matrix
   * @param[in]  i Row number i
   * @param[out] gamma vector with the weights
   */

  void Compute_Parameter_For_Linearity_Preservation_Hanging(FEMatrix& system_matrix,
                                                    double& gamma, const int i)
  {
    #ifdef __2D__
      Output::print<4>("AFC: compute parameter for linearity preservation for hanging nodes");
      //Set row pointers and column indexes
      const int* ColInd = system_matrix.get_vector_columns();
      const int* RowPtr = system_matrix.get_row_ptr();
      // get fe space
      auto fespace = system_matrix.GetFESpace2D();
      //ActiveBound : Number of non-Dirichlet and non-hanging nodes
      int ActiveBound = fespace->get_n_active_non_hanging();
      //N_Hanging : Number of hanging nodes
      int N_Hanging = fespace->get_n_hanging();
      //HangingBound : Number of non-Dirichlet nodes including hanging nodes
      int HangingBound = ActiveBound + N_Hanging;
      // allocate memory for the parameter
      double parameter_1 = 0.0, parameter_2 = 4711.0;
      //Coordinates of the points that needs to be computed
      double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, x3 = 0.0, y3 = 0.0;
      double dist_max = 0.0, dist_max2 = 0.0, dist_min = 0.0;
      double edge = 0.0, perp_dist = 0.0;
      //Loop over all rows
      //Get position of node x_i
      fespace->GetDOFPosition(i, x0, y0);
      //Loop over columns of row i
      for(int j = RowPtr[i]; j < RowPtr[i+1]; j++)
      {
        int col_j = ColInd[j];
        //If the column is a hanging node column skip
        if(col_j >=ActiveBound && col_j < HangingBound)
          continue;
        //Get position of node x_j
        fespace->GetDOFPosition(col_j, x1, y1);
        //Only non-hanging nodes(includes Dirichlet as well)
        for(int k = j + 1; k< RowPtr[i+1]; k++)
        {
          int col_k = ColInd[k];
          //Skip edge if it includes node x_i
          if(col_j != i && col_k != i)
          {
            //Position of node x_k
            fespace->GetDOFPosition(col_k, x2, y2);
            int m = 0;
            int sign = 0, sign_zero = 0;
            bool different_sign = false;
            //Loop over all entries
            for(int l = RowPtr[i]; l<RowPtr[i+1]; l++)
            {
              int col_l = ColInd[l];
              //Skip if node x_j and x_k
              if(col_l == col_j || col_l == col_k)
                continue;
              else
              {
                fespace->GetDOFPosition(col_l, x3, y3);
                sign+= check_sign(x3, y3, x1, y1, x2, y2);
                if(check_sign(x3, y3, x1, y1, x2, y2) == 0)
                  ++sign_zero;
                ++m;
                if(std::abs(sign)+sign_zero != m)
                {
                  different_sign = true;
                  break;
                }
              }
            }//End for loop l
            if(different_sign)
              continue;
            else
            {
              dist_max = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
              dist_max2 = std::sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2));
              if(dist_max2 > dist_max)
                dist_max = dist_max2;
              //Distance of point x_j and x_k
              edge = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
              //Intermediate calc for finding perpendicular distance between x_i and line x_j and x_k
              perp_dist = std::abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1));
              //Colinear point x_i, x_j, and x_k then skip
              if(perp_dist == 0)
                continue;
              else
                dist_min = perp_dist/edge;
              if(dist_max > parameter_1)
                parameter_1 = dist_max;
              if(dist_min < parameter_2)
                parameter_2 = dist_min;
            }//End of else to check different sign
          }//End for if condition for checking hanging node
        }//End for loop k
      }//End for loop j
      gamma = parameter_1/parameter_2;
       #endif
    #ifdef __3D__
      ErrThrow("The new implementation of Gamma_i not modified for 3D problems!");
    #endif
  }

/** ************************************************************************ */  
  
  /**
   * Compute the weights for the linearty preserving limiter from 
   * Barrenechea, John, Knobloch; M3AS 2017.
   *
   * @param[in] FEMatrix system matrix
   * @param[in] current solution vector
   * @param[out] gamma vector with the weights
   *
   * NOTE: The first step of [BJK17], Rem. 6.2 is not performed, i.e., no 
   *       vertex is shifted. For non-convex patches, the denominator of gamma_i
   *       might become too small, i.e. gamma_i too large. This might lead to a 
   *       smearing of layers but does not change the linearity preservation
   */

  void Compute_Parameter_For_Linearity_Preservation(FEMatrix& system_matrix,
                                                    const std::vector<double>& u,
                                                    std::vector<double>& gamma)
  {
#ifdef __2D__
    int i, j, index, N_Cells, N_Unknowns, N_V;
    double area, edge, x[3], y[3], dist, dist_max, dist_max2;
    FE_type CurrentElement;

    Output::print<4>("AFC: compute parameter for linearity preservation");
    // get fe space
    auto fespace = system_matrix.GetFESpace2D();
    // get collection
    auto Coll = fespace->GetCollection();
    // number of mesh cells
    N_Cells = Coll->GetN_Cells();
    // number of unknowns
    N_Unknowns = u.size();
    // allocate memory for the parameter
    std::vector<double> parameter(N_Unknowns);
    parameter.resize(2*N_Unknowns, 4711.0);

    // loop over the mesh cells
    for(i=0;i<N_Cells;i++)
    {
      auto cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      // get pointer to local dofs
      auto dof = fespace->GetGlobalDOF(i);
      // finite element on the mesh cell
      CurrentElement = fespace->get_fe_type(i);
      // only for P_1
      if (CurrentElement != C_P1_2D_T_A)
      {
        ErrThrow("Compute_Parameter_For_Linearity_Preservation only implemented"
                  "for P_1 !!!");
      }
      area = cell->GetMeasure();
      // loop over the vertices
      for (j=0;j<N_V;j++)
      {
        // get coordinates
        cell->GetVertex(j)->GetCoords(x[j], y[j]);
      }
      for (j=0;j<N_V;j++)
      {
        index = dof[j];
        // check for obtuse angles with opposite edge
        // vertex j+1
        if ((x[(j+2)%N_V] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
          + (y[(j+2)%N_V] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]) <=0)
        {
          dist = std::sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
            + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
          dist_max = std::sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
            + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
        }
        else
        {
          // check for obtuse angles with opposite edge
          // vertex j+2
          if ((x[(j+1)%N_V] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
            + (y[(j+1)%N_V] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]) <=0)
          {
            dist = std::sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
              + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
            dist_max = std::sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
              + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
          }
          else
            // no obtuse angle
          {
            // length opposite egde
            edge = (x[(j+1)%N_V] - x[(j+2)%N_V]) *  (x[(j+1)%N_V] - x[(j+2)%N_V]);
            edge +=  (y[(j+1)%N_V] - y[(j+2)%N_V]) *  (y[(j+1)%N_V] - y[(j+2)%N_V]);
            edge = std::sqrt(edge);
            dist = 2.*area/edge;
            dist_max = std::sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
              + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V]));
            dist_max2 = std::sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
              + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V]));
            if (dist_max2 > dist_max)
              dist_max = dist_max2;
          }
        }
        // compute numerator
        if (dist_max > parameter[index])
          parameter[index] = dist_max;
        // compute denominator
        if (dist < parameter[index+N_Unknowns])
          parameter[index+N_Unknowns] = dist;
      }
    }

    // loop over the unknowns
    for (i=0;i<N_Unknowns;i++)
    {
      gamma[i] = parameter[i]/ parameter[i+N_Unknowns];
    }
#endif
#ifdef __3D__
    int i, j, index, N_Cells, N_Unknowns, N_V;
    //A,B,C,D: Coefficients of the plane oppostie the vertex (x[j],y[j],z[j])
    double dist, dist_max, dist_max1, dist_max2, dist_max3, x[4], y[4], z[4];
    double A, B, C, area;
    FE_type CurrentElement;    
    Output::print<4>("AFC: compute parameter for linearity preservation");
    auto fespace=system_matrix.GetFESpace3D();
    // get collection  
    auto Coll = fespace->GetCollection();
    // number of mesh cells
    N_Cells = Coll->GetN_Cells();
    // number of unknowns
    N_Unknowns = u.size();
    // allocate memory for the parameter
    std::vector<double> parameter(N_Unknowns);
    parameter.resize(2*N_Unknowns, 4711.0);
        
    // loop over the mesh cells
    for(i=0;i<N_Cells;i++)
    {
      auto cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      // get pointer to local dofs
      auto dof = fespace->GetGlobalDOF(i);
      // finite element on the mesh cell
      CurrentElement = fespace->get_fe_type(i);
      // only for P_1
      if (CurrentElement != C_P1_3D_T_A)
      {
        ErrThrow("Compute_Parameter_For_Linearity_Preservation only implemented"
                 "for P_1 !!!");
      }
      area=cell->GetMeasure();
      // loop over the vertices
      for (j=0;j<N_V;j++)
      {
        // get coordinates
        cell->GetVertex(j)->GetCoords(x[j], y[j],z[j]);  
      }    
      for (j=0;j<N_V;j++)
      {
        index = dof[j];
        dist_max1= std::sqrt((x[j] - x[(j+1)%N_V])*(x[(j)] - x[(j+1)%N_V])
        + (y[j] - y[(j+1)%N_V])*(y[(j)] - y[(j+1)%N_V])
        + (z[j] - z[(j+1)%N_V])*(z[(j)]-z[(j+1)%N_V]));
        dist_max2= std::sqrt((x[j] - x[(j+2)%N_V])*(x[(j)] - x[(j+2)%N_V])
        + (y[j] - y[(j+2)%N_V])*(y[(j)] - y[(j+2)%N_V])
        + (z[j] - z[(j+2)%N_V])*(z[(j)]-z[(j+2)%N_V]));
        dist_max3= std::sqrt((x[j] - x[(j+3)%N_V])*(x[(j)] - x[(j+3)%N_V])
        + (y[j] - y[(j+3)%N_V])*(y[(j)] - y[(j+3)%N_V])
        + (z[j] - z[(j+3)%N_V])*(z[(j)]-z[(j+3)%N_V]));
        dist_max=std::max(dist_max1,dist_max2); 
        dist_max=std::max(dist_max, dist_max3);
        A=(z[(j+3)%N_V]-z[(j+1)%N_V])*(y[(j+2)%N_V]-y[(j+1)%N_V])
   -(z[(j+2)%N_V]-z[(j+1)%N_V])*(y[(j+3)%N_V]-y[(j+1)%N_V]);
        B=(z[(j+2)%N_V]-z[(j+1)%N_V])*(x[(j+3)%N_V]-x[(j+1)%N_V])
   -(z[(j+3)%N_V]-z[(j+1)%N_V])*(x[(j+2)%N_V]-x[(j+1)%N_V]);
        C=(y[(j+3)%N_V]-y[(j+1)%N_V])*(x[(j+2)%N_V]-x[(j+1)%N_V])
   -(x[(j+3)%N_V]-x[(j+1)%N_V])*(y[(j+2)%N_V]-y[(j+1)%N_V]);
        
        // OTHER IMPLEMENTATION
        //Equation of a plane with three points (x[j+1],y[j+1],z[j+1]),
  // (x[j+2],y[j+2],z[j+2]) and (x[j+3],y[j+3],z[j+3])
        //D=-[A*x[(j+1)%N_V]+B*y[(j+1)%N_V]+C*z[(j+1)%N_V]];
        //Minimum distance between a point (x[j],y[j],z[j]) and a plane Ax+By+Cz+D=0
        //dist=abs(A*x[j]+B*y[j]+C*z[j]+D)/std::sqrt(A*A+B*B+C*C);
        
        //Idea from BJK17
        dist=(6.0*area)/std::sqrt(A*A+B*B+C*C);
        // compute numerator
        if (dist_max > parameter[index])
          parameter[index] = dist_max;
        // compute denominator
        if (dist < parameter[index+N_Unknowns])
          parameter[index+N_Unknowns] = dist;  
      }
    }
    for (i=0;i<N_Unknowns;i++)
    {
      gamma[i] = parameter[i]/ parameter[i+N_Unknowns];
    }
#endif
  }
  
  /**
   * Monolithic limiter proposed in the paper Dmitri Kuzmin; CMAME 2020.
   * Works both for FEM-FCT and AFC. Extension to FCT for CDR equation
   * is little involved as one needs to apply FCT only to the hyperbolic
   * part and then add the diffusion term. 
   * NOTE: Currently working only for steady-state in ParMooN
   *       Also maybe D can be used instead of system_matrix
   *       for MPI
   * 
   * @param[in] system_matrix Matrix A+D. Used for MPI
   * @param[in] D Artifical Diffusion Matrix
   * @param[in] raw_fluxes f_ij described in papaer
   * @param[in] solution Current solution
   * @param[in] stiffness_Entries Entries of original Galerkin matrix (a_ij)
   * 
   * @param[out] alphas Limiters
   */
  void ComputeMonolithicLimiter(FEMatrix& system_matrix,
                                FEMatrix& D,
                                std::vector<double>& raw_fluxes,
                                std::vector<double>& alphas,
                                const std::vector<double>& solution,
                                std::vector<double>& stiffness_Entries)
  {
#ifdef _MPI
    const TParFECommunicator3D& comm = system_matrix.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    Output::print<2>("Using Monolithic limiter from Kuzmin2020");
    const int * ColInd = system_matrix.get_vector_columns();
    const int * RowPtr = system_matrix.get_row_ptr();
    const int nDofs = system_matrix.get_n_rows();
    const double * D_Entries = D.GetEntries();
    // get pointers to plus number of entries.
    int col_index = 0;
    std::vector<double> umin(nDofs, 0.0), umax(nDofs, 0.0);
    for (int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      umin[i] = umax[i] = solution[i];
      // i-th row of sqmatrix
      int j0 = RowPtr[i];
      int j1 = RowPtr[i+1];
      // check values of the neighbor dofs
      for(int j = j0 ; j<j1 ; j++)
      {
       // column
       col_index = ColInd[j];
       if (solution[col_index] < umin[i])
           umin[i] = solution[col_index];
       if (solution[col_index] > umax[i])
           umax[i] = solution[col_index];
       }
    }
#ifdef _MPI
    comm.consistency_update(umin.data(), 2);
    comm.consistency_update(umax.data(), 2);
#endif
    
//     //corresponding to Eq.(64)
//     std::vector<double> bar_u_ij (N_Entries, 0.0);
//     for(int i = 0; i<nDofs; i++)
//     {
// #ifdef _MPI
//       //skip the non-master rows in MPI case
//       if(masters[i] != rank)
//         continue;
// #endif
//       int j1 = RowPtr[i];
//       int j2 = RowPtr[i+1];
//       for(int j = j1; j<j2; j++)
//       {
//          col_index = ColInd[j];
//          bar_u_ij[j] = D_Entries[j]*(solution[col_index] + solution[i]) 
//                      + stiffness_Entries[j]*(solution[col_index] - solution[i]);
//       }
//     }
    
    //Terms appearing in the brackets of Eq.(46)
    double int_flux_1 = 0.0, int_flux_2 = 0.0, int_f_ij = 0.0, original_flux = 0.0;
    double bar_u_ij, bar_u_ji, diff_product, stiff_product;
    for(int i = 0; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      int_flux_1 = 0.0;
      int_flux_2 = 0.0;
      int_f_ij = 0.0;
      int j1 = RowPtr[i];
      int j2 = RowPtr[i+1];
      for(int j = j1; j<j2 ; j++)
      {
        int transposed_entry_ij = 0;
        col_index = ColInd[j];
        
        for(int k = RowPtr[col_index]; k<RowPtr[col_index + 1]; k++)
        {
          if(ColInd[k] == i)
          {
            transposed_entry_ij = k;
            break;
          }
        }
        
        //corresponding to Eq.(64) multiplied by 2*d_ij
        diff_product = D_Entries[j]*(solution[col_index] + solution[i]);
        stiff_product = stiffness_Entries[j]*(solution[col_index] - solution[i]);
        bar_u_ij = diff_product + stiff_product;
        
        stiff_product = stiffness_Entries[transposed_entry_ij] * 
                        (solution[i] - solution[col_index]);
        bar_u_ji = diff_product + stiff_product;
        
        original_flux = raw_fluxes[j];
        if(raw_fluxes[j] > 0)
        {
          int_flux_1 = (bar_u_ij - 2.0*D_Entries[j]*umax[i]);
          int_flux_2 = (2*D_Entries[j]*umin[col_index] - bar_u_ji );
          if(int_flux_1 > int_flux_2)
            int_f_ij = int_flux_2;
          else
            int_f_ij = int_flux_1;
          
          if(raw_fluxes[j] > int_f_ij)
            raw_fluxes[j] = int_f_ij;
          alphas[j] = raw_fluxes[j]/original_flux;
        }
        else if(raw_fluxes[j] < 0)
        {
          int_flux_1 = ( bar_u_ij - 2.0*D_Entries[j]*umin[i]);
          int_flux_2 = (2*D_Entries[j]*umax[col_index] - bar_u_ji);
          if(int_flux_1 > int_flux_2)
            int_f_ij = int_flux_1;
          else
            int_f_ij = int_flux_2;
          
          if(raw_fluxes[j] < int_f_ij)
            raw_fluxes[j] = int_f_ij;
          
          alphas[j] = raw_fluxes[j]/original_flux;
        }
        /**
        * This case should not arise for stead-state. If the values are very small
        * they can be cosidered as numerical approximation errors.
        */
//           if(std::fabs(alphas[j]) < 1e-13 && std::fabs(alphas[j]) != 0)
//           {
//             Output::print<4>("Negative alphas(", i,",",ColInd[j],") = ", alphas[j]);
//             alphas[j] = 0.0;
//           }
        //alphas[j] = 1.0;
     }
    }
  }
  /**
  * Monolithic limiter proposed in the paper Dmitri Kuzmin; CMAME 2020 for balanced equations
  * The original matrix is epsilon*B+A, where B is the diffusion matrix, and A is the convective
  * +reactive part
  * 
  * @param[in] system_matrix Matrix (epsilon*B+A)
  * @param[in] diffusion_matrix Diffusion Matrix (A)
  * @param[out] raw_fluxes f_ij described in paper
  * @param[in] solution Current solution
  * @param[in] poisson_sol Solution to the problem -Delta(q)=f
  * @param[in] coeffs Bilinear coeffecients
  * @param[in] is_not_afc_fixed_point_rhs If FPR or not
  * 
  * @param[out] alphas Limiters
  */
  void ComputeMonolithicLimiterSteady(FEMatrix& system_matrix,
                                      FEMatrix& diffusion_matrix,
                                      std::vector<double>& raw_fluxes,
                                      std::vector<double>& alphas,
                                      const std::vector<double>& solution,
                                      BlockVector& poisson_sol,
                                      double* coeffs,
                                      const int is_not_afc_fixed_point_rhs)
  {
#ifdef _MPI
    const TParFECommunicator3D& comm = system_matrix.GetFESpace3D()->get_communicator();
    const int* masters = comm.GetMaster();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    Output::print<3>("Using Monolithic limiter from Kuzmin2020 for Balance laws");
    //Get pointers on the column index, row pointer, and Entries
    const int * ColInd = system_matrix.get_vector_columns();
    const int * RowPtr = system_matrix.get_row_ptr();
    const int nDofs = system_matrix.get_n_rows();
    const int N_Entries = system_matrix.get_n_entries();
    //Entries of epsilon*B + A
    double * matrix_Entries = system_matrix.GetEntries();
    //Entries of matrix B
    const double * diffusion_matrix_Entries = diffusion_matrix.GetEntries();
    //Epsilon
    double diffusion_coeffecient = coeffs[0];
    //f_ij
    raw_fluxes.resize(N_Entries, 0.0);
    
    //These computations needs to be done exactly once
    if(is_not_afc_fixed_point_rhs)
    {
      //Convective part of the matrix (A)
      //a_ij = (epsilon*b_ij + a_ij) - epsilon*b_ij
      convective_matrix.resize(N_Entries, 0.0);
      for(int i = 0; i< N_Entries; i++)
        convective_matrix[i] = matrix_Entries[i] - diffusion_coeffecient*diffusion_matrix_Entries[i];
      AlgebraicFluxCorrection::correct_dirichlet_hanging_rows(diffusion_matrix);
      //Compute the artificial diffusion matrix(D), the reason it's here because we only apply
      //the artifical diffusion to the convective part and we assemble it only once
      artificial_diffusion_entries.resize(N_Entries, 0.0);
      Output::print<3>("Computing artifical diffusion matrix, D");
      for(int i = 0 ; i<nDofs; i++)
      {
        for(int j = RowPtr[i]; j<RowPtr[i+1]; j++)
        {
          int col_index = ColInd[j];
          double k_ij = 0.0;
          double k_ji = 0.0;
          if(col_index != i)
          {
            k_ij = convective_matrix[j];
            //now get the transposed entry
            for(int ll = RowPtr[col_index] ; ll < RowPtr[col_index + 1]; ll++)
            {
              if(ColInd[ll] == i)
              {
                k_ji = convective_matrix[ll];
                break;
              }
            }
            //Determine entry d_ij, max{|b_ij|, |b_ji|}
            artificial_diffusion_entries[j] = -std::max({k_ij, k_ji, 0.0});
          }
        }//end of loop j
      }//end of loop i
      //compute diagonal entries of matrix D
      for(int i = 0; i<nDofs; i++)
      {
        double val = 0.0;
        //add all entries of i-th row
        int entry_found = -1;
        for(int j = RowPtr[i]; j<RowPtr[i+1]; j++)
        {
          val += artificial_diffusion_entries[j];
          int index = ColInd[j];
          if(index == i)
            entry_found = j;
        }
        artificial_diffusion_entries[entry_found] = -val;
      }//end of loop i
    }//end of if case

    
    //Compute u_min and u_max without bar states
    std::vector<double> umin(nDofs, 0.0), umax(nDofs, 0.0);
    for (int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      umin[i] = umax[i] = solution[i];
      // i-th row of sqmatrix
      int j0 = RowPtr[i];
      int j1 = RowPtr[i+1];
      // check values of the neighbor dofs
      for(int j = j0 ; j<j1 ; j++)
      {
         // column
         int col_index = ColInd[j];
         if (solution[col_index] < umin[i])
             umin[i] = solution[col_index];
         if (solution[col_index] > umax[i])
             umax[i] = solution[col_index];
      }//end of loop j
    }//end of loop i
#ifdef _MPI
    comm.consistency_update(umin.data(), 2);
    comm.consistency_update(umax.data(), 2);
#endif
    
    //compute and store the bar states
    std::vector<double> scaled_bar_u_ij(N_Entries, 0.0);
    //Compute bar states
    //scaled_bar_u_ij = 2*d_ij(u_i + u_j) - a_ij(u_j-u_i) + b_ij(q_j-q_i)
    for(int i = 0; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      int j1 = RowPtr[i];
      int j2 = RowPtr[i+1];
      for(int j = j1; j<j2 ; j++)
      {
        int col_index = ColInd[j];
        double diff_product = 0.0, stiff_product = 0.0, poisson_product = 0.0;
        if(i == col_index)
          continue;
        //(b_ij)(q_j-q_i)
        poisson_product = diffusion_matrix_Entries[j]*(poisson_sol[col_index] - poisson_sol[i]);
        //(d_ij)(u_i+u_j)
        diff_product = artificial_diffusion_entries[j]*(solution[col_index] + solution[i]);
        //(a_ij)(u_j-u_i)
        stiff_product = convective_matrix[j]*(solution[col_index] - solution[i]);
        
        scaled_bar_u_ij[j] = diff_product + stiff_product - poisson_product;
      }//end of loop j
    }//end of loop i

    //New implementation started
    std::vector<double> sum_bar_states(nDofs, 0.0);
    for(int i = 0; i<nDofs; i++)
    {
      int j1 = RowPtr[i], j2 = RowPtr[i+1];
      double sum = 0.0;
      int diagonal_entry = 0;
      for(int j = j1; j<j2; j++)
      {
        int col_index = ColInd[j];
        if (i == col_index)
        {
          diagonal_entry = j;
          continue;
        }//end of if
        sum += scaled_bar_u_ij[j];
      }//end of loop j
      sum_bar_states[i] = sum/(-2*artificial_diffusion_entries[diagonal_entry]);
    }//end of loop i

    //Find minimum and maximum of sum_bar_states
    for(int i = 0; i<nDofs; i++)
    {
      int j1 = RowPtr[i], j2 = RowPtr[i+1];
      double min = 10000, max = -10000;
      for(int j = j1; j<j2; j++)
      {
        int col_index = ColInd[j];
        if(sum_bar_states[col_index] < min)
          min = sum_bar_states[col_index];
        if(sum_bar_states[col_index] > max)
          max = sum_bar_states[col_index];
      }//end of loop j
      if(min < umin[i])
        umin[i] = min;
      if(max < umax[i])
        umax[i] = max;
    }//end of loop i
    
    //Find correct fluxes
    for(int i = 0; i < nDofs; i++)
    {
      for(int j = RowPtr[i]; j < RowPtr[i+1]; j++)
      {      
        int col_index = ColInd[j];
        //Terms appearing in the brackets of Eq.(46)
        double int_flux_1 = 0.0, int_flux_2 = 0.0, f_ij_max = 0.0, f_ij_min = 0.0;
        //Transposed Entry
        int transposed_entry_ij = 0;
        for(int k = RowPtr[col_index]; k<RowPtr[col_index + 1]; k++)
        {
          if(ColInd[k] == i)
          {
            transposed_entry_ij = k;
            break;
          }
        }//end of loop k
        double right_hand_side = 0.0;
        if(i == col_index)
          continue;
        //(b_ij)(q_j-q_i)
        right_hand_side = diffusion_matrix_Entries[j]*(poisson_sol[col_index] - poisson_sol[i]);
        //f_ij = d_ij(u_j-u_i)
        raw_fluxes[j] = artificial_diffusion_entries[j]*(solution[col_index]- solution[i]);
        if(raw_fluxes[j] > 0)
        {
          int_flux_1 = scaled_bar_u_ij[j] - 2.0*artificial_diffusion_entries[j]*umax[i];
          int_flux_2 = 2.0*artificial_diffusion_entries[j]*umin[col_index] - scaled_bar_u_ij[transposed_entry_ij];
          if(int_flux_1 > int_flux_2)
            f_ij_max = int_flux_2;
          else
            f_ij_max = int_flux_1;
          
          raw_fluxes[j] = std::min(f_ij_max, raw_fluxes[j]);
          raw_fluxes[j] = std::max(0.0, raw_fluxes[j]);
        }
        else if(raw_fluxes[j] < 0)
        {
          int_flux_1 = scaled_bar_u_ij[j] - 2.0*artificial_diffusion_entries[j]*umin[i];
          int_flux_2 = 2.0*artificial_diffusion_entries[j]*umax[col_index] - scaled_bar_u_ij[transposed_entry_ij];
          if(int_flux_1 > int_flux_2)
            f_ij_min = int_flux_1;
          else
            f_ij_min = int_flux_2;
          raw_fluxes[j] = std::max(f_ij_min, raw_fluxes[j]);
          raw_fluxes[j] = std::min(0.0, raw_fluxes[j]);
        }
        raw_fluxes[j] += right_hand_side;
      }//end of loop j
    }//end of loop i
    //Add the artificial diffusion matrix
    if(is_not_afc_fixed_point_rhs)
    {
      //epsilon*A + B +D
      for(int i = 0; i<nDofs; i++)
      {
        for(int j= RowPtr[i]; j<RowPtr[i+1]; j++)
        {
          //int col_index = ColInd[j];
          matrix_Entries[j] += artificial_diffusion_entries[j];
        }
      }
    }

    
  }//end of Monolithic limiter
}



/** ************************************************************************ */

ParameterDatabase AlgebraicFluxCorrection::default_afc_database()
{
  ParameterDatabase db("default algebraic flux correction database");

  // Type of AFC to be applied.
  db.add("algebraic_flux_correction", "none", " Chose which type of afc to use.",
    {"none", "afc", "fem-fct-cn"}
  );
    
  db.add("compute_cut_lines", "no", "Compute the width of layer (yes, no)", 
          {"no", "yes"});

    // type of the limiter to be applied for steady-state case
  db.add("afc_limiter", "kuzmin", "Choose an afc limiter. Options are"
    "kuzmin, BJK17, constant, monolithic, monolithic_steady, MUAS, MUAS_Kno21, MUAS_MAX, MUAS_MAX_ABS", 
    {"kuzmin", "BJK17", "constant", "monolithic", "monolithic_steady", "MUAS", "MUAS_Kno21", 
    "MUAS_MAX", "MUAS_MAX_ABS"}
  );

  // type of the limiter to be applied for steady-state case
  db.add("afc_initial_iterate", "afc_zero", "Choose the initial iterate"
    "afc_zero, galerkin, supg, upwind", {"afc_zero", "galerkin", 
       "supg", "upwind"}
  );

  // iteration scheme for the afc methods
  db.add("afc_iteration_scheme", "fixed_point_matrix", 
         "Choose an iteration scheme for the afc methods. Options are"
    "fixed_point_rhs, fixed_point_matrix, newton, newton_regu",
    {"fixed_point_rhs", "fixed_point_matrix", "newton", 
      "newton_regu", "newton_no_damp"}
  );
  
  //maximum number of iterations in the non linear loop
  db.add("afc_nonlinloop_maxit", (size_t) 1 ,
         "Maximal number of iterations for the nonlinear loop in AFC."
         "Must be a value between 0 and 100000.", (size_t)  0, (size_t)100000);
  
  //tolerance for non linear loop in AFC scheme
  db.add("afc_nonlinloop_epsilon", 1e-10, 
         "Stopping criterion for the nonlinear loop in AFC."
         "Must be a value between 1e-20 an 1e20.", 1e-20, 1e20);
  

  //Constants related to Dynamic Damping from [JK08]
  db.add("afc_nonlinloop_damping_factor", 1.0, "A damping parameter"
          "for the nonlinear loop in AFC. Must be a value between 1" 
          "(no damping) and 0 (no update).", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_max", 1.0, 
         "Maximal number for afc_nonlinloop_damping_factor."
         "Only changed internally", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_min", 0.01, 
         "Minimal number for afc_nonlinloop_damping_factor."
         "Intended to stay constant", 0.0,1.0);
  
  db.add("afc_nonlinloop_damping_factor_max_global", 1.0, 
         "Maximal number for afc_nonlinloop_damping_factor."
         "for complete iteration", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_increase", 1.1, 
         "Increase factor for afc_nonlinloop_damping_factor"
         "Intended to stay constant", 1.0,2.0);        // c_2 in [JK08]

  db.add("afc_nonlinloop_damping_factor_decrease", 0.5, 
         "Decrease factor for afc_nonlinloop_damping_factor"
         "Intended to stay constant", 0.0,1.0);

  db.add("afc_nonlinloop_damping_factor_max_increase", 1.001, 
         "Increase factor for afc_nonlinloop_damping_factor_max."
         "Intended to stay constant", 1.0,2.0);        // c_3 in [JK08]

  db.add("afc_nonlinloop_damping_factor_max_decrease", 0.9, 
         "Decrease factor for afc_nonlinloop_damping_factor_max."
         "Intended to stay constant", 0.0,1.0);        // c_4 in [JK08]

  db.add("afc_nonlinloop_damping_factor_min_tol", 1.001, 
         "Tolerance for afc_nonlinloop_damping_factor_min."
         "Intended to stay constant", 1.0,2.0);        // c_1 in [JK08]

  db.add("afc_nonlinloop_damping_factor_constant", "no", 
         "Whether or not afc_nonlinloop_damping_factor"
         "should be constant",                         //
  {
    "no", "yes"
  });
  
  //Anderson acceleration parameters
  db.add("afc_nonlinloop_anderson_acc", "no", 
         "Whether or not to use Anderson acceleration", {"no", "yes"}
  );
  
  db.add("afc_nonlinloop_anderson_acc_vec", (size_t) 10, 
         "Number of vectors in Anderson acceleration",  
         (size_t)  1, (size_t)  1000);
  
  db.add("afc_nonlinloop_anderson_acc_start", (size_t) 0, 
         "Starting iterate for Anderson acceleration",  
         (size_t)  0, (size_t)  1000);
  
  //parameters related to fixed point matrix and Newton method
  db.add("afc_newton_regu_sigma", 1e-8, "Penalty for regularized Newton method"
    "is scaled with fourth power of mesh width", 0.0,1.0);

  db.add("afc_fixed_point_matrix_weight", 0.75, 
         "weight for implicit fixed point method"
         "fixed_point_matrix", 0.0,1.0);

  db.add("afc_fixed_point_derivative_weight", 0.1, 
         "weight for contribution of the derivative"
         "in Newton's method", 0.0,1.0);
  
  db.add("afc_fixed_point_derivative_weight_factor", 0.0, 
         "factor for weight for contribution of the derivative"
         "in Newton's method", 0.0,1.0);
   
  db.add("afc_damping_bound_newton", 0.02, 
         "bound for switching from fixed_point_rhs to Newton", 0.02, 1.0);
 
  db.add("afc_change_method_threshold", 1e-5, 
         "threshold for changing to Newton's method", 0.0,100.0);

  /*db.add("afc_nonlinloop_switch_to_newton_scheme", 1, 
   * "scheme for switching to formal Newton's method"
   *    "1 - standard, 10 - first fpr (for BAIL proceedings), 
   *    11 - first fpm (for BAIL proceedings)"
   *    "20 - first fpr with reswitch (for BAIL proceedings), 
   *    21 - first fpm with reswitch (for BAIL proceedings)"
   *    ,{1,10,11,20,21});*/
  
  //projection to admissible values. Idea from [BB17.CMAME] 
  db.add("afc_project_to_admissible_values", "no", 
         "Whether or not to project intermediate iterates to admissible values",
         {"no", "yes"});
   
  return db;
}


// ////////////////////////////////////////////////////////////////////////////
// Implementation of the methods in the namespace AlgebraicFluxCorrection    //
// ////////////////////////////////////////////////////////////////////////////

void AlgebraicFluxCorrection::steady_state_algorithm(
FEMatrix& system_matrix,
FEMatrix& diffusion_matrix, 
const std::vector<double>& sol,
std::vector<double>& rhs,
BlockVector& poisson_sol,
const std::vector<int>& neum_to_diri,
double* coeffs,
FEMatrix& D,
FEMatrix& D_B,
std::vector<double>& gamma,
std::vector<double>& alphas,
bool compute_D_and_gamma,
const ParameterDatabase& db,
Limiter limiter,
Iteration_Scheme it_scheme,
const int is_not_afc_fixed_point_rhs)
{
  Output::print<4>("AFC: enter steady_state_algorithm");
  //catch non-square matrix
  if (!system_matrix.is_square())
  {
    ErrThrow("System matrix must be square for AFC!");
  }
  
#ifdef _MPI
  const TParFECommunicator3D& comm = system_matrix.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  // store the total number of dofs
  int nDofs = system_matrix.get_n_rows();
// #ifdef _MPI
//     auto space = system_matrix.GetFESpace3D();
//     nDofs = space->get_communicator().get_n_global_dof();
// #endif
  // heritage style index declaration
  int i,j0,j1,j2,j3,jj,index;
  double alpha_ij;

  // get pointers to columns, rows and entries of matrix A
  const int * ColInd = system_matrix.get_vector_columns();
  const int * RowPtr = system_matrix.get_row_ptr();
  // non-const, matrix entries get modified!
  double * Entries = system_matrix.GetEntries();
  double * D_Entries = D.GetEntries();
  double * D_B_Entries = D_B.GetEntries();


  int N_Entries = system_matrix.get_n_entries();
  // allocate memory for matrix F and flux limiters
  double* F = new double[N_Entries+6*nDofs];
  memset(F, 0, (N_Entries+6*nDofs)*sizeof(double));
  double* P_plus = F + N_Entries;
  double* P_minus = P_plus + nDofs;
  double* Q_plus = P_minus + nDofs;
  double* Q_minus = Q_plus + nDofs;
  double* R_plus = Q_minus + nDofs;
  double* R_minus = R_plus + nDofs;
  
  raw_fluxes.resize(N_Entries, 0.0);
  
  std::vector<double> umin, umax, q;  
 
  //B_Entries : Matrix B for MUAS-type methods
  //modified_flux: Fluxes computed for MUAS-type methods
  std::vector<double> B_Entries, modified_flux;
  B_Entries.resize(N_Entries, 0.0);
  modified_flux.resize(N_Entries, 0.0);
  //hanging node related
  //determine the lower bound for hanging nodes
  size_t hangingLowBound;
  int N_hanging;
#ifdef __3D__
  
  N_hanging = system_matrix.GetFESpace3D()->get_n_hanging();
  hangingLowBound = system_matrix.GetFESpace3D()->get_n_active_non_hanging();
  auto fespace = system_matrix.GetFESpace3D();
#elif __2D__
  
  N_hanging = system_matrix.GetFESpace2D()->get_n_hanging();
  hangingLowBound = system_matrix.GetFESpace2D()->get_n_active_non_hanging();
  auto fespace = system_matrix.GetFESpace2D();
#endif

  if ((it_scheme == Iteration_Scheme::NEWTON)||
    (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
    ||(it_scheme == Iteration_Scheme::NEWTON_REGU)
    ||(it_scheme == Iteration_Scheme::FIXEDPOINT_RHS))
    alphas.resize(N_Entries+1,0.0);
  // compute entries of the artificial diffusion matrix D
  // TODO make matrix D an actual TMatrix and not only an entries vector
  if (compute_D_and_gamma)
  {
    if(limiter == Limiter::BJK17)
      modify_system_matrix(system_matrix);	  
    Output::print<4>("AFC: compute matrix D");
    if(limiter != Limiter::MONOLITHIC_STEADY)
      compute_artificial_diffusion_matrix(system_matrix, D);
    if (limiter == Limiter::BJK17)
    {
      Output::print<4>("AFC: compute vector gamma");
      Compute_Parameter_For_Linearity_Preservation(system_matrix, sol, gamma);
      //In case of hanging nodes the computation of gamma_i takes a lot of time, because one
      //needs to loop over the entire row twice. One way around is to compute the gamma_i for
      //all i and then only for those i's which have some entry because of hanging DOFs.
      //This isn't the most efficient way but it gets the work done. Also, use of hanging nodes
      //might not be used in the future (specifically for BJK limiter)
      if(N_hanging != 0)
      {
        for(int i = 0; i<nDofs; i++)
        {
          for(int j = RowPtr[i]; j < RowPtr[i+1] ; j++)
          {
            int col_j = ColInd[j];
            if( col_j >= (int)hangingLowBound && col_j < (int)(hangingLowBound + N_hanging))
            {
                Compute_Parameter_For_Linearity_Preservation_Hanging(system_matrix, gamma[i], i);
                break;
            }//End of if case
          }//End of loop j
        }//End of loop i
      }
    }
  }
  
  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  //for fixed_point_rhs and iteration>1 we don't need to add matrix D again.
  if (is_not_afc_fixed_point_rhs)
  {
    original_stiffness_Entries = system_matrix.get_entries();
    //for balance laws the matrix D changes as it only contains the convective part
    if(limiter != Limiter::MONOLITHIC_STEADY)
      system_matrix += D;
  }
  /*Previous Implementation
    Daxpy(N_Entries, 1.0, &afc_matrix_D_entries[0], Entries);*/
  
  // allocate and fill arrays for linearity preserving limiter
  // from Barrenechea, John, Knobloch M3AS (2017)
  
  if (limiter == Limiter::BJK17)
  {
    umin.resize(nDofs,0.0);
    umax.resize(nDofs,0.0);
    q.resize(nDofs,0.0);

    for (int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      umin[i] = umax[i] = sol[i];
      q[i] = 0;
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // check values of the neighbor dofs
      for(int j = j0 ; j < j1 ; j++)
      {
        // column
        index = ColInd[j];
        if (sol[index] < umin[i])
          umin[i] = sol[index];
        if (sol[index] > umax[i])
          umax[i] = sol[index];
        if (i != index)
          q[i] += D_Entries[j];
      }
      q[i] *= gamma[i];
    }
  }

  // compute matrix F
  // loop over all rows
  for(int i = 0 ; i < nDofs ; i++)
  {
#ifdef _MPI
    //skip the non-master rows in MPI case
    if(masters[i] != rank)
        continue;
#endif
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(int j = j0 ; j < j1 ; j++)
    {
      // column
      index = ColInd[j];
      // d_ij (u_j - u_i)
      F[j] = D_Entries[j] * (sol[index]-sol[i]);
      raw_fluxes[j] = F[j];
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  // linearity preserving limiter from [BJK17]
  if (limiter == Limiter::BJK17)
  {
    for(int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // loop over the neighbors
      for(int j = j0 ; j < j1 ; j++)
      {
        // column
        index = ColInd[j];
        // diagonal
        if (index == i)
          continue;
        if (F[j] > 0)
          P_plus[i] += F[j];
        if (F[j] < 0)
          P_minus[i] += F[j];
      }
      Q_plus[i] = q[i]*(sol[i]-umax[i]);
      Q_minus[i] = q[i]*(sol[i]-umin[i]);
    }
  }

  // MUAS-type limiters
  // NOTE: MPI not working with MUAS-type limiters
  if ((limiter == Limiter::MUAS)||(limiter == Limiter::MUAS_Kno21)||
        (limiter == Limiter::MUAS_MAX)||(limiter == Limiter::MUAS_MAX_ABS))
  {
    for(int i = 0 ; i < nDofs ; i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(int j = j0 ; j < j1 ; j++)
      {
        index = ColInd[j];
        // only the active part of the matrix
        if(original_stiffness_Entries[j] > 0)
        {
          if(sol[i]-sol[index] > 0)
            P_plus[i] += original_stiffness_Entries[j]*(sol[i]-sol[index]);
          else
            P_minus[i] += original_stiffness_Entries[j]*(sol[i]-sol[index]);
        }
        if (limiter == Limiter::MUAS)
	{
            if ((sol[index]-sol[i]) > 0)
               Q_plus[i] += std::abs(original_stiffness_Entries[j])*(sol[index]-sol[i]);
            else
               Q_minus[i] += std::abs(original_stiffness_Entries[j])*(sol[index]-sol[i]);
	}
	if (limiter == Limiter::MUAS_Kno21)
        {
            if (F[j] > 0)
               Q_minus[i] -= F[j];
            if (F[j] < 0)
               Q_plus[i] -= F[j];
	} 
	if (limiter == Limiter::MUAS_MAX)
        {
            double val_q;
	    int jjj0,jjj1;
            val_q = std::abs(original_stiffness_Entries[j]);
	    // find transposed entry
	    jjj0 = RowPtr[index];
            jjj1 = RowPtr[index+1];
            for(int jjj = jjj0 ; jjj < jjj1 ; jjj++)
	    {
	       if ((ColInd[jjj]) == i)
               {
		  if (original_stiffness_Entries[jjj] > val_q)
		  {
		     val_q = original_stiffness_Entries[jjj];
		     break;
		  }	  
	       }
	    }	    
            if ((sol[index]-sol[i]) > 0)
               Q_plus[i] += val_q*(sol[index]-sol[i]);
            else
               Q_minus[i] += val_q*(sol[index]-sol[i]);
        }
        if (limiter == Limiter::MUAS_MAX_ABS)
        {
            double val_q;
            int jjj0,jjj1;
            val_q = std::abs(original_stiffness_Entries[j]);
            // find transposed entry
            jjj0 = RowPtr[index];
            jjj1 = RowPtr[index+1];
            for(int jjj = jjj0 ; jjj < jjj1 ; jjj++)
            {
               if ((ColInd[jjj]) == i)
               {
                  if (std::abs(original_stiffness_Entries[jjj]) > val_q)
                  {  
                     val_q = std::abs(original_stiffness_Entries[jjj]);
                     break;
                  }
               }       
            } 
            if ((sol[index]-sol[i]) > 0)
               Q_plus[i] += val_q*(sol[index]-sol[i]);
            else
               Q_minus[i] += val_q*(sol[index]-sol[i]);
        }
      }
    }
  }

  // Kuzmin limiter
  if (limiter == Limiter::KUZMIN)
  {
    for(int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(int j = j0 ; j < j1 ; j++)
      {
        // VJ: removed these lines 18/05/09 since they affect 
        // just diagonal entries
        // if ((it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) 
        //      && (Entries[j] > 0))
        // continue;
        // column
        index = ColInd[j];
        // check transposed entry -> jj
        // diagonal
        if (index == i)
          continue; 
        j2 = RowPtr[index];
        j3 = RowPtr[index+1];
        for (jj=j2;jj<j3;jj++)
        {
          if (ColInd[jj] == i)
          {
            break;
          }
        }
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        // note that a[jj] > a[j] <==> a[jj]+d[jj] > a[j]+d[j] 
        //  since d is a symmetric matrix
#ifndef _MPI
        if (Entries[jj] > Entries[j])
          continue;
        // only the active part of the matrix
        if (F[j] > 0)
        {
          P_plus[i] += F[j];
          if (index < nDofs)
            Q_plus[index] += F[j];// commented by Partl
          else
            ErrThrow("The index cannot be bigger than nDOFs");
          Q_minus[i] -= F[j];
        }
        if (F[j] < 0)
        {
          P_minus[i] += F[j];
          Q_plus[i] -= F[j];
          if (index < nDofs)
            Q_minus[index] +=  F[j];
          else
            ErrThrow("The index cannot be bigger than nDOFs");// commented by Partl
        }
#else
        if (original_stiffness_Entries[jj] + D_Entries[j] <= Entries[j])
        {
          if (F[j] > 0)
          {
            P_plus[i] += F[j];
            Q_minus[i] -= F[j];
          }
          if (F[j] < 0)
          {
            P_minus[i] += F[j];
            Q_plus[i] -= F[j];
          }
        }
        else
        {
          double Fjj = D_Entries[j] * ( sol[i] - sol[index] );
          
          if(Fjj  > 0)
            Q_plus[i] += Fjj;
          
          if(Fjj < 0)
            Q_minus[i] +=  Fjj;
        }
#endif
        
        
      }                                           // end loop j
    }
  }

  // apply the nodal correction factor evaluated at the upwind node i
  // loop over all nodes

  // original but discontinuous proposal
  // other propsals in MooNMD
  for(int i = 0; i < nDofs ; i++)
  {
#ifdef _MPI
    //skip the non-master rows in MPI case
    if(masters[i] != rank)
       continue;
#endif
    // initialization
    R_plus[i] = R_minus[i] = 1.0;
    if (std::abs(P_plus[i])>0)
    {
      R_plus[i] = Q_plus[i]/P_plus[i];
      if (R_plus[i] >1)
        R_plus[i] = 1;
    }
    if (std::abs(P_minus[i])>0)
    {
      R_minus[i] = Q_minus[i]/P_minus[i];
      if (R_minus[i] >1)
        R_minus[i] = 1;
    }
  }

  // treat Dirichlet nodes
  for (int j = 0; j < (int) neum_to_diri.size() ; j++)
  {
    i = neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }
  
#ifdef _MPI
  // put R_plus and R_minus into level 2 consistency
  comm.consistency_update(R_plus, 2);
  comm.consistency_update(R_minus, 2);
#endif

  // store the flux limiters
  // loop over all rows
  for(int i = 0; i < nDofs ; i++)
  {
#ifdef _MPI
    //skip the non-master rows in MPI case
    if(masters[i] != rank)
      continue;
#endif
    int diag_index = i;
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    // loop over the columns
    for(int j = j0 ; j < j1 ; j++)
    {
      // column
      index = ColInd[j];
      // diagonal entry
      if (index == i)
      {
        diag_index = j;
        continue;
      }
      /*// this should not be happen
      if (Entries[j] > 0)
      {
        ErrThrow("positive non-diagonal entry in AFC ", i, " ", j, " ",
          Entries[j]);
      }
      if ((it_scheme == Iteration_Scheme::FIXEDPOINT_RHS) && (Entries[j] > 0))
      {
        ErrThrow("positive entry in AFC ", i, " ", j, " ", Entries[j]);
      }*/

      // original, symmetric application
      if (limiter == Limiter::BJK17)
      {
        alpha_ij = 1;
        if (F[j] > 0)
        {
          alpha_ij = R_plus[i];
          if (R_minus[index] < alpha_ij)
            alpha_ij = R_minus[index];
        }
        if (F[j] < 0)
        {
          alpha_ij = R_minus[i];
          if (R_plus[index] < alpha_ij)
            alpha_ij = R_plus[index];
        }

        //storage of limiters
        alphas[j]=alpha_ij;
      }

      // MUAS-type methods
      if ((limiter == Limiter::MUAS)||(limiter == Limiter::MUAS_Kno21)||
        (limiter == Limiter::MUAS_MAX)||(limiter == Limiter::MUAS_MAX_ABS))
      {
        alpha_ij = 1;
        if(sol[i]-sol[index] > 0)
          alpha_ij = R_plus[i];
        else if (sol[i]-sol[index] < 0)
          alpha_ij = R_minus[i];
        
        //Storage of limiters
        alphas[j] = alpha_ij;
      }

      if (limiter == Limiter::KUZMIN)
      {
        // check transposed entry
        j2 = RowPtr[index];
        j3 = RowPtr[index+1];
        for (jj=j2;jj<j3;jj++)
        {
          // index of transposed entry
          if (ColInd[jj]==i)
          {
            break;
          }
        }
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        // commented by Partl
#ifndef _MPI
        if (Entries[jj] > Entries[j])
          continue;

        alpha_ij = 1;
        if (F[j] > 0)
        {
          alpha_ij = R_plus[i];
        }
        if (F[j] < 0)
        {
          alpha_ij = R_minus[i];
        }
        // storage of limiters
        alphas[j] = alpha_ij;
        alphas[jj] = alpha_ij;
#else
        if(original_stiffness_Entries[jj] + D_Entries[j] <= Entries[j])
        {
          alpha_ij = 1;
          if (F[j] > 0)
          {
            alpha_ij = R_plus[i];
          }
          if (F[j] < 0)
          {
            alpha_ij = R_minus[i];
          }
          alphas[j] = alpha_ij;
        }
        else
        {
          double Fjj = D_Entries[j] * ( sol[i] - sol[index] );
          
          alpha_ij = 1;
          if (Fjj  > 0)
          {
            alpha_ij = R_plus[index];
          }
          if (Fjj < 0)
          {
            alpha_ij = R_minus[index];
          }
          
          alphas[j] = alpha_ij;
        }
#endif
      }                                           //End of Kuzmin limiter
    }                                             //End of For loop j
    alphas[diag_index]=1.0;
  }                                               //End of loop i
  if(db["afc_limiter"].is("monolithic"))
  {
    ComputeMonolithicLimiter(system_matrix, D, raw_fluxes, alphas ,sol,
                             original_stiffness_Entries);
  }

  //New implementation started for Monolithic limiter for balance laws
  if(db["afc_limiter"].is("monolithic_steady"))
  {
    ComputeMonolithicLimiterSteady(system_matrix, diffusion_matrix, 
                                   raw_fluxes, alphas ,sol,
                                   poisson_sol,
                                   coeffs, 
                                   is_not_afc_fixed_point_rhs);
  }

  if ((limiter == Limiter::MUAS)||(limiter == Limiter::MUAS_Kno21)||
	(limiter == Limiter::MUAS_MAX)||(limiter == Limiter::MUAS_MAX_ABS))
  {
    ComputeMatrixB(system_matrix, alphas, original_stiffness_Entries, 
                  B_Entries);
    for(int i = 0; i<nDofs; i++)
    {
      for(int j = RowPtr[i]; j<RowPtr[i+1]; j++)
      {
        index = ColInd[j];
        D_B_Entries[j] = D_Entries[j] - B_Entries[j];
        modified_flux[j] = D_B_Entries[j]*(sol[index] - sol[i]);
        //Set limiters to one now for a Posteriori estimator. This will cancel out the
        //edge stabilisation term
        alphas[j] = 1.0;
      }
    }
  }
  
  //Update in the RHS for different limiters
  if (it_scheme == Iteration_Scheme::FIXEDPOINT_RHS)
  {
#ifndef _MPI
    int index, jj;
#endif
    for(int i = 0; i < nDofs; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      int j1 = RowPtr[i];
      int j2 = RowPtr[i+1];
      for(int j = j1; j<j2; j++)
      {
        if(limiter == Limiter::KUZMIN)
        {
#ifndef _MPI         
          index = ColInd[j];
          if(index == i)
              continue;
          // check transposed entry
          int j3 = RowPtr[index];
          int j4 = RowPtr[index+1];
          for (jj = j3 ; jj < j4 ; jj++)
          {
            // index of transposed entry
            if (ColInd[jj]==i)
              break;
          }
          // check upwind condition
          // this ensures that the 'link' between i and index is treated only once
          if (Entries[jj] > Entries[j])
            continue;
          rhs[i] += alphas[j]*F[j];
          // update rhs wrt to current column
          // note that F[j] = -F[jj] and 
          //alpha_j = alpha_jj (symmetry of alpha matrix)
          if (index < nDofs)// commented by Partl
            rhs[index] -= alphas[j]*F[j];
          else
            ErrThrow("The index cannot be bigger than nDOFs");
#else
          rhs[i] += alphas[j]*F[j];    
#endif
        }
        else if(limiter == Limiter::BJK17)
          rhs[i] += alphas[j]*F[j];
        else if(limiter == Limiter::MONOLITHIC || limiter == Limiter::MONOLITHIC_STEADY)
          rhs[i] += raw_fluxes[j];
        else if((limiter == Limiter::MUAS)||(limiter == Limiter::MUAS_Kno21)||
        (limiter == Limiter::MUAS_MAX)||(limiter == Limiter::MUAS_MAX_ABS))
          rhs[i] += modified_flux[j];
        else
          ErrThrow("Unknown limiter to the AFC scheme!");
      }
    }
  }

  if (it_scheme == Iteration_Scheme::FIXEDPOINT_MATRIX)
  {
    int j3, j4, index_j;
    double tau0 =  (double)db["afc_fixed_point_matrix_weight"];
    double tau1 = 1.0-tau0;
    // update matrix
    for(int i = 0; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      j3 = RowPtr[i];
      j4 = RowPtr[i+1];
      
      // loop over the columns of the matrix
      for(int j = j3 ; j < j4 ; j++)
      {
        index_j=ColInd[j];
        //Non-Diagonal Entries
        if(i!=index_j)
        {
          Entries[j]-=alphas[j]*D_Entries[j]*tau0;
          // if (Entries[j] > 0)
          //  Output::print<2>("non-diag pos ", Entries[j]);
          rhs[i] += alphas[j]*D_Entries[j]*sol[index_j]*tau1;
        }
        else
        {
          double sum_flux=0.0;
          for(int jj=j3;jj<j4;jj++)
          {
            if(ColInd[jj]!=i)
              sum_flux += alphas[jj]*D_Entries[jj];
          }
          Entries[j]+=sum_flux*tau0;
          //if (Entries[j] < 0)
          //  Output::print<2>("diag neg ", Entries[j]);

//           rhs[i] += alphas[j]*D_Entries[j]*sol[index_j]*tau1;
          rhs[i] -= sum_flux*sol[i]*tau1;
        }
      }
    }                                             //Formation of matrix complete
  }

  if ((it_scheme == Iteration_Scheme::NEWTON)
      || (it_scheme == Iteration_Scheme::NEWTON_REGU))
  {
    int j3, j4, index_j;
    std::vector<double> df(N_Entries,0.0);
    double tau0 = (double)db["afc_fixed_point_matrix_weight"];
    double tau1 = (double)db["afc_fixed_point_derivative_weight"]
                 *(double)db["afc_fixed_point_derivative_weight_factor"];
    if (db["afc_iteration_scheme"].is("newton_no_damp"))
    {
      if ((double)db["afc_fixed_point_derivative_weight_factor"]>1e-3)
      {
        tau0 =  tau1 = 1.0;
      }
    }
    Output::print<4>("tau0/tau1 ", tau0 , " " , tau1);
    
    // compute first part of rhs
    // with old matrix
    for(int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      j3 = RowPtr[i];
      j4 = RowPtr[i+1];
      for(int j = j3 ; j < j4 ; j++)
      {
        index_j = ColInd[j];
        rhs[i] += -(Entries[j])*sol[index_j]+alphas[j]*F[j];
      }
    }

    Output::print<4>("AFC: computing Jacobian");
    // compute Jacobian, store on matrix entries
    if (limiter == Limiter::KUZMIN)
    {
      if (it_scheme == Iteration_Scheme::NEWTON)
      {
        // Computation of Matrix DF
        // loop over the degrees of freedom
        for(int i = 0 ; i < nDofs ; i++)
        {
#ifdef _MPI
         //skip the non-master rows in MPI case
         if(masters[i] != rank)
           continue;
#endif
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          // loop over the columns of the matrix
          for(int j=j3;j<j4;j++)
          {
            index_j=ColInd[j];
            //Non-Diagonal Entries
            if(i!=index_j)
            {
              //Derivative Product is a function which returns the summation 
              //inside the DF matrix
              df[j]=Entries[j]-tau0*alphas[j]*D_Entries[j]
                -tau1*Compute_Jacobian_times_flux_Kuzmin(system_matrix,
                                                         F,P_plus,Q_plus,
                                                         Q_minus,P_minus,
                                                         R_minus,R_plus,i,
                                                         index_j,j,
                                                         D);
            }
            else
            {
              double sum_flux=0.0;
              for(int jj=j3;jj<j4;jj++)
              {
                if(ColInd[jj]!=i)
                  sum_flux += alphas[jj]*D_Entries[jj];
              }
              df[j]=Entries[j]+tau0*sum_flux
                -tau1*Compute_Jacobian_times_flux_Kuzmin(system_matrix,
                                                         F,P_plus,Q_plus,
                                                         Q_minus,P_minus,
                                                         R_minus,R_plus,i,
                                                         index_j,j,
                                                         D);
            }
          }               //end of loop j
        }                 //end of loop i     //Formation of matrix DF complete
      }                   //end of Newton case

      if (it_scheme == Iteration_Scheme::NEWTON_REGU)
      {
        // Computation of Matrix DF
        // loop over the degrees of freedom
        double sigma = (double)db["afc_newton_regu_sigma"];
        
        for(int i = 0 ; i < nDofs ; i++)
        {
#ifdef _MPI
         //skip the non-master rows in MPI case
          if(masters[i] != rank)
            continue;
#endif
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          // loop over the columns of the matrix
          for(int j=j3;j<j4;j++)
          {
            index_j=ColInd[j];
            // Entries contains already a+d
            df[j]=Entries[j]
             -Compute_Jacobian_times_flux_Kuzmin_Regularized(system_matrix,
                                                             F,i,index_j,j,
                                                             D,
                                                             sigma, tau0, tau1);
          }
        }                //end of loop i      //Formation of matrix DF complete
      }                  //end of Newton_regu
    }                    //end of KUZMIN limiter case

    if(limiter == Limiter::BJK17)
    {
      if (it_scheme == Iteration_Scheme::NEWTON)
      {
        //Computation of matrix DF
        //Loop over the degrees of freedom
        for(int i = 0 ; i < nDofs ; i++)
        {
#ifdef _MPI
          //skip the non-master rows in MPI case
          if(masters[i] != rank)
           continue;
#endif
          j3=RowPtr[i];
          j4=RowPtr[i+1];
          //Loop over the columns
          for(int j = j3 ; j < j4 ; j++)
          {
            index_j=ColInd[j];
            //Non-Diagonal Entries
            if(i!=index_j)
            {
              //Derivative Product is a function which returns the summation inside the DF matrix
              df[j]=Entries[j]
                   -tau0*alphas[j]*D_Entries[j]
                   -tau1*Compute_Jacobian_times_flux_BJK17(system_matrix, 
                                                           F, P_plus,P_minus,
                                                           Q_plus,Q_minus,
                                                           R_plus,R_minus,
                                                           umax,umin,q,
                                                           sol[index_j],i,
                                                           index_j,j,
                                                           D);
            }
            else
            {
              double sum_flux=0.0;
              for(int jj = j3 ; jj < j4 ; jj++)
              {
                if(ColInd[jj] != i)
                  sum_flux +=alphas[jj]*D_Entries[jj];
              }
              df[j]=Entries[j]+tau0*sum_flux
                   -tau1*Compute_Jacobian_times_flux_BJK17(system_matrix, 
                                                           F, P_plus,P_minus,
                                                           Q_plus,Q_minus,
                                                           R_plus,R_minus,
                                                           umax,umin,q,
                                                           sol[index_j],i,
                                                           index_j,j,
                                                           D);
            }
          }                //end of loop j
        }                  //end of loop i  //Formation of matrix DF complete
      }                    //end of Newton case
    }                      //end of BJK17 limiter case

    // update rhs for system to be solved
    for(int i = 0 ; i < nDofs ; i++)
    {
#ifdef _MPI
      //skip the non-master rows in MPI case
      if(masters[i] != rank)
        continue;
#endif
      j3=RowPtr[i];
      j4=RowPtr[i+1];
      for(int j = j3 ; j < j4 ; j++)
      {
        index_j = ColInd[j];
        rhs[i] += df[j]*sol[index_j];           //Using DF instead of matrix A
        Entries[j]=df[j];
      }
    }
  }

  delete [] F;                                    
}                                               //End of function definition

void AlgebraicFluxCorrection::correct_dirichlet_hanging_rows(FEMatrix& MatrixA)
{
  //hold pointers to row, kcol, entries array
  const int* RowPtr_A      = MatrixA.get_row_ptr();
  const int* KCol_A        = MatrixA.get_vector_columns();
  double* Entries_A  = MatrixA.GetEntries();

  //determine first and one-after-last dirichlet rows
  size_t diriHighBound = MatrixA.get_n_rows();
  size_t diriLowBound = diriHighBound - MatrixA.GetTestSpace()->get_n_dirichlet();
  
  // loop over rows and set them to unity-vectors
  for (size_t rowIndex = diriLowBound;
    rowIndex < diriHighBound ;++rowIndex)
  {
    int l0 = RowPtr_A[rowIndex];
    int l1 = RowPtr_A[rowIndex+1];
    for (int l=l0;l<l1;l++)
    {
      // diagonal entry
      if (KCol_A[l]== (int) rowIndex)
        Entries_A[l] = 1;
      else
        Entries_A[l] = 0;
    }
  }
  
  MatrixA.correct_hanging_rows();
}

/** *********************************************************************** */

void AlgebraicFluxCorrection::correct_dirichlet_hanging_rhs
(FEMatrix& Matrix, BlockVector& RHS)
{
  size_t hangingLowBound;
  int N_hanging;
#ifdef __3D__
  N_hanging = Matrix.GetFESpace3D()->get_n_hanging();
  hangingLowBound = Matrix.GetFESpace3D()->get_n_active_non_hanging();
#elif __2D__
  N_hanging = Matrix.GetFESpace2D()->get_n_hanging();
  hangingLowBound = Matrix.GetFESpace2D()->get_n_active_non_hanging();
#endif
  for(int i = 0; i<N_hanging ; i++)
    RHS[hangingLowBound+i] = 0.0;
}

/** ************************************************************************ */

void AlgebraicFluxCorrection::AFC_Compute_New_Iterate(
  const BlockVector& old_solution,
  BlockVector& new_solution,
  const ParameterDatabase& db,
  FEMatrix& matrix)
{
  new_solution.add_scaled(old_solution,-1.0);
  new_solution.scale(db["afc_nonlinloop_damping_factor"]);
  new_solution.add_scaled(old_solution,1.0);

#ifdef __2D__
  auto fespace = matrix.GetFESpace2D();
#elif __3D__
  auto fespace = matrix.GetFESpace3D();
#endif

  int N_Dirichlet = fespace->get_n_dirichlet();
  size_t diriLowBound = fespace->get_n_dof() - N_Dirichlet;
  for(int i=0; i<N_Dirichlet; i++)
  {
    new_solution[diriLowBound+i] = old_solution[diriLowBound+i];
  }

  //Note: Projection only possible when the maximum and minimum are already known.
  if (db["afc_project_to_admissible_values"].is("yes"))
  {
    Output::print<2>("  projection to admissible values");

    int len = new_solution.length();
    for (int i=0;i<len;i++)
    {
      if (new_solution[i]<0)
        new_solution[i] = 0;
      if (new_solution[i]>1)
        new_solution[i] = 1;
    }
  }
}

