// =======================================================================
// @(#)Assembler4.C        1.16
//
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Alfonso Caiazzo & Laura Blank
//
// History:     start of implementation 10.08.2016
//
// =======================================================================

#include <Assembler4.h>
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

#include <string.h>
#include <stdlib.h>
#include <vector>



//================================================================================
Assembler4::Assembler4()
 : qf_ref(QuadratureFormula_type::BaryCenterTria), qf_orig(qf_ref) // dummy type
{
  Coll = nullptr;
  square_matrices.resize(0);
  rectangular_matrices.resize(0);
  rhs_blocks.resize(0);
  n_square_matrices = 0;
  n_rectangular_matrices = 0;
  n_rhs_blocks = 0;
  n_all_matrices = 0;
  maximum_number_base_functions = 0;
}

//================================================================================
void Assembler4::init(BlockFEMatrix &M,
    BlockVector &rhs,
    std::vector<const TFESpace2D*>& fespaces,
    std::vector<const TFESpace2D*>& ferhs)
{
  // set the collection
  this->Coll = fespaces[0]->GetCollection();

  this->maximum_number_base_functions = 0;

  /// this loop sets the CellIndex. Probably it should be done somewhere else?
  for (int i = 0; i < this->Coll->GetN_Cells(); i++)
  {
    TBaseCell *cell = this->Coll->GetCell(i);
    cell->SetCellIndex(i);

    // find the maximum number of basis functions used by the current collection
    for (unsigned int j = 0; j < fespaces.size(); j++)
    {
      auto fe = fespaces[j]->get_fe(i);
      if (fe.GetN_DOF() > this->maximum_number_base_functions)
      {
        this->maximum_number_base_functions = fe.GetN_DOF();
      }
    }
  }
  // --------------------------------------------------------------------
  Output::print("Assembler: maximum_number_base_functions = ", maximum_number_base_functions);


  // set vector of blocks
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  // Note: This class is supposed to work only for NSType 14 at the moment
  if(blocks.size() == 9)
  {
    this->n_square_matrices = 5;
    square_matrices.resize(n_square_matrices);
    square_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    square_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
    square_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
    square_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
    square_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

    this->n_rectangular_matrices = 4;
    rectangular_matrices.resize(n_rectangular_matrices);
    rectangular_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get());
    rectangular_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
    rectangular_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
    rectangular_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
  }
  else if(blocks.size() == 1)
  {
    this->n_square_matrices = 1;
    square_matrices.resize(n_square_matrices);
    square_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

    this->n_rectangular_matrices = 0;
    //rectangular_matrices.resize(n_rectangular_matrices);
  }
  else
  {
    Output::print("Assembler should be only used with NSType 14 Matrices");
    ErrThrow(" --> Wrong blocks.size() = ", blocks.size());
  }

  this->n_rhs_blocks = ferhs.size();
  rhs_blocks.resize(n_rhs_blocks);

  for (size_t k = 0; k < n_rhs_blocks; k++)
  {
    rhs_blocks[k] = rhs.block(k);
  }

  this->n_all_matrices = this->n_square_matrices + this->n_rectangular_matrices;
}


//================================================================================
//================================================================================
void Assembler4::Assemble2D(BlockFEMatrix &M,
    BlockVector &b_rhs,
    std::vector<const TFESpace2D*>& fespaces,
    std::vector<const TFESpace2D*>& ferhs,
    const Example2D& example,
    const std::vector<std::shared_ptr<LocalAssembling2D>>& la_list,
    int AssemblePhaseID)
{
#ifdef __3D__
  ErrThrow("Assembler4::Assembler4() not yet working in 3D");
#endif

  // set matrices and rhs blocks
  this->init(M, b_rhs, fespaces, ferhs);

  // LocRhs: an array of pointers (size: number of rhs)
  double **LocRhs;
  // righthand: a big pointer
  double *righthand;
  if (n_rhs_blocks)
  {
    LocRhs = new double* [n_rhs_blocks];
    righthand = new double [n_rhs_blocks* maximum_number_base_functions];
    for (size_t i = 0; i < n_rhs_blocks; i++)
      LocRhs[i] = righthand + i * maximum_number_base_functions;
  }

  // set pointers to matrices
  double *aux;
  ///@todo all these pointers are used but might not be initialized
  double **Matrices;
  double ***LocMatrices;

  if (this->n_all_matrices)
  {
    aux = new double
        [this->n_all_matrices * maximum_number_base_functions * maximum_number_base_functions];
    Matrices = new double* [this->n_all_matrices*maximum_number_base_functions];

    for (int j = 0; j < (this->n_all_matrices * maximum_number_base_functions); j++)
      Matrices[j] = aux + j * maximum_number_base_functions;

    // LocMatrices[i] = a matrix of size NBaseFct x NBaseFct
    LocMatrices = new double** [this->n_all_matrices];

    for (int i = 0; i < this->n_all_matrices; i++)
      LocMatrices[i] = Matrices + i * maximum_number_base_functions;
  }

  //================================================================================
  // loop over all cells
  //================================================================================
  for (int i = 0; i < this->Coll->GetN_Cells(); i++)
  {
    TBaseCell *cell = this->Coll->GetCell(i);

    // only for multiphase flows
    if ((AssemblePhaseID >= 0) &&
        (AssemblePhaseID != cell->GetPhase_ID()) )
      continue;
    for (size_t num_localAss = 0; num_localAss < la_list.size(); num_localAss++)
    {
      //assemble the selected local form on the i-th cell
      this->assemble_local_system(fespaces, i, LocMatrices, LocRhs, la_list[num_localAss]);

      // add local/cellwise matrices to global matrices
      //(ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
      this->add_local_to_global_matrix(i, LocMatrices, Matrices);

      // add local/cellwise right-hand sides to global right-hand sides
      this->add_local_to_global_rhs(i, ferhs, LocRhs, example);

      // boundary condition part
      this->impose_boundary_conditions(i, ferhs, example);
    }
  }


  // modify matrix according to coupling of hanging nodes
  // this part is only relevant to the case with hanging nodes and
  // it could be put in a separate function
  this->handle_hanging_nodes(ferhs);

  // delete pointers
  if (n_rhs_blocks)
  {
    delete [] righthand;
    delete [] LocRhs;
  }

  if (this->n_all_matrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }
}


//================================================================================
// 
//================================================================================
void Assembler4::assemble_local_system(std::vector <const TFESpace2D*>& fespaces,
    int i, double ***LocMatrices, double **LocRhs,
    std::shared_ptr<LocalAssembling2D> la)
{
  TBaseCell *cell = this->Coll->GetCell(i);

  // find local (on this cell) used elements
  std::vector<const FiniteElement*> LocalUsedElements;
  int n_fespaces = fespaces.size();
  std::vector<const BaseFunctions*> LocBF(n_fespaces);

  for (int j = 0; j < n_fespaces; j++)
  {
    auto& element = fespaces[j]->get_fe(i);
    LocalUsedElements[j] = &element;
    LocBF[j] = element.GetBaseFunct();
  }

  // calculate values on original element
  bool *SecondDer = la->GetNeeds2ndDerivatives();
  FEDatabase::GetOrig(LocalUsedElements, this->Coll, cell, SecondDer, qf_ref,
                      qf_orig);

  bool is_sdfem = (la->get_disctype() == 2); // SDFEM

  if (is_sdfem || (TDatabase::ParamDB->CELL_MEASURE == 4))
  {
    TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
    int N_Edges = cell->GetN_Edges();
    for (int ij = 0; ij < N_Edges; ij++)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
      TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    }
    if (N_Edges==3)
      TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
  }

  la->GetLocalForms(qf_orig, LocBF, cell, i, n_all_matrices, n_rhs_blocks,
                    LocMatrices, LocRhs);
}

//================================================================================
// impose boundary conditions
//================================================================================
void Assembler4::impose_boundary_conditions(int i_cell,
    std::vector<const TFESpace2D*>& ferhs,
    const Example2D& example)
{
  TBaseCell *cell = this->Coll->GetCell(i_cell);

  for (unsigned int j = 0; j < rhs_blocks.size(); j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    const int *DOF = ferhs[j]->GetGlobalDOF(i_cell);
    double *RHS = this->rhs_blocks[j];

    BoundCondFunct2D *BoundaryCondition = fespace->get_boundary_condition();
    BoundValueFunct2D *BoundaryValue = example.get_bd()[j];
    auto ele = fespace->get_fe(i_cell);
    double t0, t1;
    const TBoundComp *BoundComp;
    int N_EdgePoints;
    const double *EdgePoints;
    double eps = 1e-4;

    auto nf = ele.GetNodalFunctional();

    if (TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
    {
      ErrThrow("What is SUPERCONVERGENCE_ORDER?? most likely not working");
      /* Superconvergence boundary interpolation */
//       if (nf->GetID() == NF_C_Q_Q2_2D)
//         nf = FEDatabase::GetNodalFunctional2D(NF_S_Q_Q2_2D);
    }

    nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

    const FEDescriptor *FEDesc_Obj = ele.GetFEDesc();
    int  N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
    // setting Dirichlet boundary condition
    BoundCond Cond0, Cond1;

    for (int m = 0; m < cell->GetN_Joints(); m++)
    {
      TJoint *joint = cell->GetJoint(m);

      if (joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge ||
          joint->GetType() == InterfaceJoint)
      {
        if (joint->GetType() == BoundaryEdge||
            joint->GetType() == InterfaceJoint)
        {
          TBoundEdge *boundedge = (TBoundEdge *)joint;
          BoundComp = boundedge->GetBoundComp();
          boundedge->GetParameters(t0, t1);
        }
        else
        {
          TIsoBoundEdge *isoboundedge = (TIsoBoundEdge *)joint;
          BoundComp = isoboundedge->GetBoundComp();
          isoboundedge->GetParameters(t0, t1);
        }
        // get id of the boundary component
        int comp = BoundComp->GetID();
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

        double FunctionalValues[MaxN_BaseFunctions2D];
        double PointValues[MaxN_PointsForNodal2D];

        // only one boundary condition per edge allowed
        if(Cond0 == Cond1)
        {
          switch(Cond0)
          {
          case DIRICHLET:
          {
            // if DG
            if (N_EdgeDOF == 0)
              break;
            // read boundary values for each quadrature point
            for (int l = 0; l < N_EdgePoints; l++)
            {
              double s = EdgePoints[l];
              double t = 0.5 * (t0 * (1-s) + t1 * (1+s));
              BoundaryValue(comp, t, PointValues[l]);
            }

            // compute boundary values for each dof on the
            // boundary edge with the nodal functionals
            nf->GetEdgeFunctionals(this->Coll, cell, m, PointValues,
                FunctionalValues);
            int *EdgeDOF = FEDesc_Obj->GetJointDOF(m);
            // save boundary values of each dof on the boundary
            // edge in the rhs
            for (int l = 0; l < N_EdgeDOF; l++)
            {
              RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
            }
            break;
          }
          case NEUMANN:
          case ROBIN:
            Output::print<4>(" ** Assembler: WARNING: NEUMANN and ROBIN boundary conditions are not supported. ** ");
            Output::print<4>(" ** You can try to use the old Assemble2D(...) instead ** ");
            break;
          case DIRICHLET_WEAK:
            break;
          case PERIODIC:
            break;
          default :
            ErrThrow("Unknown boundary condition !");
          }                                     // endswitch Cond0
        }                                       // end if (Cond0 == Cond1)
        else
        {
          ErrThrow("different boundary condition on one edge are not allowed!");
        }
      }                                         // endif (boundary joint)
    }
  }
}


//===================================================================
// add local/cellwise matrices to global matrices
//(ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
//===================================================================
void Assembler4::add_local_to_global_matrix(int i,
    double ***LocMatrices,
    double **Matrices)
{
  bool assemble_dirichlet_rows = false;
  // square matrices
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    double **Matrix = LocMatrices[j];
    auto fespace = square_matrices[j]->GetFESpace2D();
    int DirichletBound = fespace->get_n_active_non_hanging();

    /**
         DOF[k] is the global index of the k-th local degree of freedom
         MatrixRow[k] is the assembled value corresponding to the m-th
         local test function and k-th local ansatz function. That means it
         corresponds to the l=DOF[m]-th global test function and the
         DOF[k]-th global ansatz function
     */
    //int *DOF = GlobalNumbers[j] + BeginIndex[j][i];
    ///@todo use a vector<int> dof = fespace.local_to_global_DOFs(i,cell)
    const int *DOF = fespace->GetGlobalDOF(i);
    int n_local_dof = fespace->get_n_local_dof(i);
    /**
         The list of degrees of freedom consists of 3 parts:
         0,...,ActiveBounds-1: active nodes (no hanging, no Dirichlet)
         ActiveBounds, DirichletBounds-1: Hanging nodes
         DirichletBound,...end : the DOF with DIRICHLET condition
     */

    // add local matrix to global matrix
    for (int m = 0; m < n_local_dof; m++)
    {
      // active DOF
      if (DOF[m] < DirichletBound || assemble_dirichlet_rows)
      {
        for (int k = 0; k < n_local_dof; k++)
          square_matrices[j]->add(DOF[m], DOF[k], Matrix[m][k]);
      }
      else
      {
        //Dirichlet standard treatment: 1.0 on the diagonal
        square_matrices[j]->set(DOF[m], DOF[m], 1.0);
      }
    }                                           // endfor m
  }                                             // endfor j




  // rectangular matrices
  for (int j = 0; j < this->n_rectangular_matrices; j++)
  {
    auto TestElement = rectangular_matrices[j]->GetTestSpace2D()->get_fe(i);
    auto AnsatzElement = rectangular_matrices[j]->GetAnsatzSpace2D()->get_fe(i);

    int N_Test = TestElement.GetN_DOF();
    int N_Ansatz = AnsatzElement.GetN_DOF();

    double **Matrix = Matrices + (j+this->n_square_matrices) * maximum_number_base_functions;

    double *Entries = rectangular_matrices[j]->GetEntries();
    const int *RowPtr = rectangular_matrices[j]->get_row_ptr();
    const int *ColInd = rectangular_matrices[j]->get_vector_columns();

    //int *TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
    //int *AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];
    const int *TestDOF = rectangular_matrices[j]->GetTestSpace()->GetGlobalDOF(i);
    const int *AnsatzDOF = rectangular_matrices[j]->GetAnsatzSpace()->GetGlobalDOF(i);

    // add local matrix to global
    for (int m = 0; m < N_Test; m++)
    {
      int end = RowPtr[ TestDOF[m] + 1];

      for (int n = RowPtr[ TestDOF[m] ]; n < end; n++)
      {
        for (int k = 0; k < N_Ansatz; k++)
        {
          if (AnsatzDOF[k] == ColInd[n])
          {
            Entries[n] += Matrix[m][k];
            break;
          }
        }
      }
    }
  }
}


//================================================================================
// add local right-hand sides to global right-hand side
//================================================================================
void Assembler4::add_local_to_global_rhs(int i,
    std::vector<const TFESpace2D*>& ferhs,
    double **LocRhs,
    const Example2D&)
{
  for (unsigned int j = 0; j < rhs_blocks.size(); j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    int N_ = fespace->get_n_local_dof(i);

    // get the pointer to j-th RHS from the global righthandside
    // local_rhs = LocRhs[j] ?
    //double *local_rhs = righthand+j*maximum_number_base_functions;
    double *RHS_block_j = rhs_blocks[j];
    // find space for this linear form

    int DirichletBound = fespace->get_n_active_non_hanging();

    // dof of the rhs nodes connected to this cell
    const int *DOF = ferhs[j]->GetGlobalDOF(i);
    //RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

    // add local right-hand side to the global one
    for (int m = 0; m < N_; m++)
    {
      int l = DOF[m];
      if (l < DirichletBound)
      {
        // node l is inner node
        RHS_block_j[l] += LocRhs[j][m];
      }
    }
  }
}



//================================================================================
// hanging nodes
//================================================================================
void Assembler4::handle_hanging_nodes(std::vector<const TFESpace2D*>& ferhs)
{
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    square_matrices[j]->ModifyMatrixAccordingToCoupling(false);
  }

  for (int j = 0; j < this->n_rectangular_matrices; j++)
  {
    // hanging nodes in test space
    // note that there used to be code here which did what 
    // ModifyMatrixAccordingToCoupling does but with the ansatz space.
    rectangular_matrices[j]->ModifyMatrixAccordingToCoupling(false);
  }

  for (size_t j = 0; j < n_rhs_blocks; j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    double *RHS = rhs_blocks[j];

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
  }


  //================================================================================
  // write coupling into matrix
  //================================================================================
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    square_matrices[j]->correct_hanging_rows();
  }
}


//================================================================================
Assembler4::~Assembler4(){}




