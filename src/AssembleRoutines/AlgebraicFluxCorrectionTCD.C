#include "AlgebraicFluxCorrectionTCD.h"
#include "math.h"


ParameterDatabase default_afc_tcd_database()
{
  ParameterDatabase db("default afc database for TCD");
  
  db.add("algebraic_flux_correction", "none", " Chose which type of afc to use.",
         {"none", "afc", "fem-fct-cn"} );

  db.add("afc_fct_scheme", "linear", " Chose which type of scheme to use for FCT.",
         {"linear", "non-linear", "explicit"} );
  
  db.add("theta_scheme", 0.5, "Choose theta for time stepping scheme. Options "
         "are 1 (Backward Euler), 0.5 (C-N)", 0.0,1.0);
  
  db.add("compute_cut_lines", "no", "Compute the width of layer (yes, no)", 
         {"no", "yes"} );
  
  db.add("afc_prelimiter", 0, "Choose an afc flux prelimiting scheme. Options "
         "are 0 (none), 1 (min-mod), 2 (grad-direction), 3 (both)", {0,1,2,3} );
  
  db.add("afc_limiter", "zalesak", "Choose an afc limiter. Options are"
         "zalesak, monolithic", {"zalesak", "monolithic"} );
  
  db.add("afc_nonlinloop_maxit", (size_t) 1 ,
         "Maximal number of iterations for the nonlinear loop in AFC."
         "Must be a value between 0 and 100000.", (size_t)  0, (size_t)100000 );
  
  db.add("grad_direction_tolerance", 1.0e-16, 
         "Tolerance for comparison with zero in grad-direction limiter."
         "Must be a value between 0.0 an 1e-10.", 0.0, 1.0e-10 );
  
  db.add("afc_nonlinloop_epsilon", 1e-10, 
         "Stopping criterion for the nonlinear loop in AFC."
         "Must be a value between 1e-20 an 1e20.", 1e-20, 1e20 );
  
  db.add("projection_onto_admissible_interval", "no",
         "Project solution onto given interval of admissible values (yes, no)", 
         {"no", "yes"} );
  
  db.add("admissible_minimum", 0.0, "Lower endpoint of admissible interval", -1.0e-10, 1.0e10);
  
  db.add("admissible_maximum", 1.0, "Upper endpoint of admissible interval", -1.0e-10, 1.0e10);
  
  db.add("solDotUsedInZalesak", false, "Whether or not to use high-order stabilization"
  " via time derivatives in definition of anti-diffusive fluxes."
  " Relevant only for nonlinear scheme.", {true,false});
  
  return db;
}


template <int d>
RowsColsRecord<d>::RowsColsRecord(const Params& params)
: indicesOfDiagEntries(params.massMatPtr->get_n_rows()),
  rows(params.massMatPtr->get_row_ptr()),
  cols(params.massMatPtr->get_vector_columns()),
  massMat(*params.massMatPtr)
{
  classifyRows();
  storeIndicesOfDiagEntries();
}


template <int d>
void RowsColsRecord<d>::classifyRows()
{
  recordMyMasterRows();
  recordMyActiveMasterRows();
  recordMyNonActiveMasterRows();
}


template <int d>
void RowsColsRecord<d>::recordMyMasterRows()
{
#ifdef _MPI
  const TParFECommunicator3D& comm = massMat.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  myMasterRows.resize(comm.GetN_Master());
  
  int i = -1;
  
  for(int row = 0; row < massMat.get_n_rows(); row++)
  {
    if(masters[row] == rank)
    {
      i++;
      myMasterRows[i] = row;
    }
  }
#else 
  myMasterRows.resize(massMat.get_n_rows());
  
  for(unsigned int row = 0; row < myMasterRows.size(); row++)
    myMasterRows[row] = row;
#endif
}


template <int d>
void RowsColsRecord<d>::recordMyActiveMasterRows()
{
#ifdef _MPI
  const TParFECommunicator3D& comm = massMat.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  myActiveMasterRows.resize(comm.GetN_Master());
  
  int i = -1;
  
  for(int row = 0; row < massMat.get_n_active_rows(); row++)
  {
    if(masters[row] == rank)
    {
      i++;
      myActiveMasterRows[i] = row;
    }
  }
  
  int numberOfMyActiveMasterRows = i + 1;
  myActiveMasterRows.resize(numberOfMyActiveMasterRows);
#else 
  myActiveMasterRows.resize(massMat.get_n_active_rows());
  
  for(unsigned int row = 0; row < myActiveMasterRows.size(); row++)
    myActiveMasterRows[row] = row;
#endif
}


template <int d>
void RowsColsRecord<d>::recordMyNonActiveMasterRows()
{
#ifdef _MPI
  const TParFECommunicator3D& comm = massMat.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  myNonActiveMasterRows.resize(comm.GetN_Master());
  
  int i = -1;
  
  for(int row = massMat.get_n_active_rows(); row < massMat.get_n_rows(); row++)
  {
    if(masters[row] == rank)
    {
      i++;
      myNonActiveMasterRows[i] = row;
    }
  }
  
  int numberOfMyNonActiveMasterRows = i + 1;
  myNonActiveMasterRows.resize(numberOfMyNonActiveMasterRows);
#else
  myNonActiveMasterRows.resize( massMat.get_n_rows() - massMat.get_n_active_rows() );
  
  int i = -1;
  for(int row = massMat.get_n_active_rows(); row < massMat.get_n_rows(); row++)
  {
    i++;
    myNonActiveMasterRows[i] = row;
  }
#endif
}


template <int d>
void RowsColsRecord<d>::storeIndicesOfDiagEntries()
{
  const int* rows = massMat.get_row_ptr();
  const int* cols = massMat.get_vector_columns();
  
  for(int row : myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if(row == col)
        indicesOfDiagEntries[row] = i;
    }
  }
}


template <int d>
TCD_BC_setter<d>::TCD_BC_setter(const Params& params)
: boundaryValuesVector(*params.rhsPtr),
  boundaryValuesFunction(*params.boundaryValuesFunction),
  boundaryConditionFunction(*params.boundaryConditionFunction),
  boundaryValuesFEFunctionEntries(params.massMatPtr->get_n_rows()),
  boundaryValuesFEFunction( params.solFESpacePtr,
                            "boundaryValuesFEFct",
                            boundaryValuesFEFunctionEntries.data() ),
  boundaryTimeDerivativeFct(params.boundaryTimeDerivativeFct),
  boundaryTimeDerivativeVector(*params.rhsPtr),
  boundaryTimeDerivativeFEFctEntries(params.massMatPtr->get_n_rows()),
  boundaryTimeDerivativeFEFct( params.solFESpacePtr,
                               "timeDerivativeFEFct",
                               boundaryTimeDerivativeFEFctEntries.data() ),
  timeStep(params.timeStep)
{}


template <int d>
void TCD_BC_setter<d>::setBoundaryValuesVectorForTime(double time)
{
  double originalTime = TDatabase::TimeDB->CURRENTTIME;
  
  TDatabase::TimeDB->CURRENTTIME = time;
  boundaryValuesFEFunction.SetDirichletBC(boundaryConditionFunction,
                                          boundaryValuesFunction);
  TDatabase::TimeDB->CURRENTTIME = originalTime;
  
  memcpy( boundaryValuesVector.get_entries(),
          boundaryValuesFEFunction.GetValues(),
          sizeof(double)*boundaryValuesFEFunction.GetLength() );
}


template <int d>
void TCD_BC_setter<d>::setBoundaryTimeDerivativeVectorForTime(double time)
{
  if(boundaryTimeDerivativeFct)
  {
    double originalTime = TDatabase::TimeDB->CURRENTTIME;
    
    TDatabase::TimeDB->CURRENTTIME = time;
    boundaryTimeDerivativeFEFct.SetDirichletBC(boundaryConditionFunction,
                                               boundaryTimeDerivativeFct);
    TDatabase::TimeDB->CURRENTTIME = originalTime;
    
    memcpy( boundaryTimeDerivativeVector.get_entries(),
            boundaryTimeDerivativeFEFct.GetValues(),
            sizeof(double)*boundaryTimeDerivativeFEFct.GetLength() );
  }
  else
  {
    setBoundaryValuesVectorForTime(time);
    boundaryTimeDerivativeVector.copy_nonactive(boundaryValuesVector);
  
    if(time - timeStep >= 0.0)
    {
      setBoundaryValuesVectorForTime(time - timeStep);
      boundaryTimeDerivativeVector.addScaledNonActive(boundaryValuesVector, -1.0);
      boundaryTimeDerivativeVector.scaleNonActive(1.0/timeStep);
    }
    else
    {
      setBoundaryValuesVectorForTime(time + timeStep);
      boundaryTimeDerivativeVector.addScaledNonActive(boundaryValuesVector, -1.0);
      boundaryTimeDerivativeVector.scaleNonActive(-1.0/timeStep);
    }
  }
}


template <int d>
TCD_limiter<d>::TCD_limiter(const Params& params)
: rowsColsRecord(*params.RowsColsRecordPtr),
  fluxes(*params.fluxesPtr),
  diffMat(*params.diffMatPtr),
  lumpMatDiagonal(*params.lumpMatDiagonalPtr),
  sol(*params.solPtr),
  bcSetter(*params.tcdBcSetterPtr)
{}


template <int d>
TCD_limiter_zalesak<d>::TCD_limiter_zalesak(const Params& params)
: TCD_limiter<d>(params),
  solutionTime(0.0),
  P_plus(params.massMatPtr->get_n_rows(), 0.0),
  P_minus(params.massMatPtr->get_n_rows(), 0.0),
  Q_plus(params.massMatPtr->get_n_rows(), 0.0),
  Q_minus(params.massMatPtr->get_n_rows(), 0.0),
  R_plus(params.massMatPtr->get_n_rows(), 0.0),
  R_minus(params.massMatPtr->get_n_rows(), 0.0),
  timeStep(params.timeStep),
  solInterm(*params.solOldPtr),
  solOld(*params.solOldPtr),
  rhsOld(*params.rhsOldPtr),
  stiffMat(*params.stiffMatPtr),
  toleranceForZero(params.toleranceForZero),
  timeSteppingScheme(params.timeSteppingScheme),
  prelimiter(params.prelimiter),
  oneMinusTheta(1.0 - params.theta)
{}



template <int d>
void TCD_limiter_zalesak<d>::calculateSolIntermExplicit(const BlockVector& solDot, const BlockVector& sol)
{
  for(int row : TCD_limiter<d>::rowsColsRecord.myActiveMasterRows)
    solInterm[row] = sol[row] + timeStep * solDot[row];
  
  TCD_limiter<d>::bcSetter.setBoundaryValuesVectorForTime(solutionTime);
  TCD_limiter<d>::bcSetter.setBCsIn(solInterm);
  
  #ifdef _MPI
  TCD_limiter<d>::updateConsistency(solInterm.get_entries(), 2);
  #endif
}


template <int d>
void TCD_limiter_zalesak<d>::prelimitFluxes()
{
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  
  const double* diffMatEntries = TCD_limiter<d>::diffMat.GetEntries();
  double* fluxes = TCD_limiter<d>::fluxes.data();
  
  double minModValue;
  double solIntermDiff;
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if( row == col )
        continue;
  
      switch(prelimiter)
      {
        case NONE:
          break;
          
        case BOTH:
        case MIN_MOD:
          minModValue = - diffMatEntries[i]*( solInterm[row] - solInterm[col] );
          fluxes[i] = minMod(2.0 * minModValue, fluxes[i]);
          break;
        
        case GRAD_DIRECTION:
          solIntermDiff = solInterm[col] - solInterm[row];
          
          if(makesSenceToCompareSignOf(solIntermDiff))
            if( fluxes[i]*( solIntermDiff ) > 0.0 )
              fluxes[i] = 0.0;
          break;
      }
    }
  }
}


template <int d>
void TCD_limiter_zalesak<d>::limitFluxes()
{
  if(prelimiter != NONE)
    prelimitFluxes();
  
  calculateRs();
  updateFluxes();
}


template <int d>
void TCD_limiter_zalesak<d>::calculateRs()
{
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  std::vector<double>& fluxes = TCD_limiter<d>::fluxes;
  
  memset(P_plus.data(), 0, P_plus.size()*sizeof(double));
  memset(P_minus.data(), 0, P_minus.size()*sizeof(double));
  memset(Q_plus.data(), 0, Q_plus.size()*sizeof(double));
  memset(Q_minus.data(), 0, Q_minus.size()*sizeof(double));
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myActiveMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if( row == col )
        continue;
      
      if(fluxes[i] > 0.0)
        P_plus[row] += fluxes[i];
      else
        P_minus[row] += fluxes[i];
      
      double solDiff = solInterm[col] - solInterm[row];
      
      if(solDiff > Q_plus[row])
        Q_plus[row] = solDiff;
      
      if(solDiff < Q_minus[row])
        Q_minus[row] = solDiff;
    }
    
    if(P_plus[row] == 0.0)
      R_plus[row] = 1.0;
    else
      R_plus[row] = std::min( 1.0, (TCD_limiter<d>::lumpMatDiagonal[row]*Q_plus[row])/( P_plus[row] * timeStep ) );
    
    if(P_minus[row] == 0.0)
      R_minus[row] = 1.0;
    else
      R_minus[row] = std::min( 1.0, (TCD_limiter<d>::lumpMatDiagonal[row]*Q_minus[row])/( P_minus[row] * timeStep ) );
  }
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myNonActiveMasterRows)
  {
    R_plus[row] = 1.0;
    R_minus[row] = 1.0;
  }
  
  #ifdef _MPI
  TCD_limiter<d>::updateConsistency(R_plus.data(), 2);
  TCD_limiter<d>::updateConsistency(R_minus.data(), 2);
  #endif
}


template <int d>
void TCD_limiter_zalesak<d>::updateFluxes()
{
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  
  double alpha;
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if(TCD_limiter<d>::fluxes[i] > 0.0)
        alpha = std::min( R_plus[row], R_minus[col] );
      else
        alpha = std::min( R_minus[row], R_plus[col] );
      
      TCD_limiter<d>::fluxes[i] *= alpha;
    }
  }
}


template <int d>
void TCD_limiter_zalesak<d>::calculateSolInterm()
{
  switch(timeSteppingScheme)
  {
    case CRANK_NICOLSON:
      calculateSolIntermNonOneTheta();
      break;
    case BACKWARD_EULER:
      calculateSolIntermBackwardEuler();
      break;
    default:
      ErrThrow("Unsupported type of timeSteppingScheme.");
  }
  
  #ifdef _MPI
  TCD_limiter<d>::updateConsistency(solInterm.get_entries(), 2);
  #endif
}


template <int d>
void TCD_limiter_zalesak<d>::calculateSolIntermNonOneTheta()
{
  calculateSolIntermNonOneThetaActiveRows();
  calculateSolIntermNonOneThetaNonActiveRows();
}


template <int d>
void TCD_limiter_zalesak<d>::calculateSolIntermNonOneThetaActiveRows()
{
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  const double* stiffMatEntries = stiffMat.GetEntries();
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myActiveMasterRows)
  {
    double matrixRowVectorProduct = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      matrixRowVectorProduct += stiffMatEntries[i]*solOld[col];
    }
    
    double bracket = rhsOld[row] - matrixRowVectorProduct;
    
    solInterm[row] = solOld[row] +
                     oneMinusTheta*timeStep*bracket/TCD_limiter<d>::lumpMatDiagonal[row];
  }
}


template <int d>
void TCD_limiter_zalesak<d>::calculateSolIntermNonOneThetaNonActiveRows()
{
  double timeInterm = TDatabase::TimeDB->CURRENTTIME
                      - oneMinusTheta*timeStep;
  
  TCD_limiter<d>::bcSetter.setBoundaryValuesVectorForTime(timeInterm);
  TCD_limiter<d>::bcSetter.setBCsIn(solInterm);
}


template <int d>
void TCD_limiter_zalesak<d>::limitFluxesExplicit(const BlockVector& solDot, const BlockVector& sol)
{
  calculateSolIntermExplicit(solDot, sol);
  
  calculateRs();
  updateFluxes();
}


template <int d>
TCD_limiter_monolithic<d>::TCD_limiter_monolithic(const Params& params)
: TCD_limiter<d>(params),
  solMax(*params.rhsPtr),
  solMin(*params.rhsPtr),
  origStiffMat(*params.stiffMatPtr),
  diffPartMatOrStiffWithoutDiffMat(*params.stiffMatPtr)
{}



template <int d>
void TCD_limiter_monolithic<d>::determineExtremesOf(const BlockVector& sol)
{
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myMasterRows)
  {
    double maxValue = sol[row];
    double minValue = sol[row];
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
       if(maxValue < sol[col])
         maxValue = sol[col];
       
       if(minValue > sol[col])
         minValue = sol[col];
    }
    
    solMax[row] = maxValue;
    solMin[row] = minValue;
  }
  
  #ifdef _MPI
  TCD_limiter<d>::updateConsistency(solMin.get_entries(), 2);
  TCD_limiter<d>::updateConsistency(solMax.get_entries(), 2);
  #endif
}


template <int d>
void TCD_limiter_monolithic<d>::limitFluxesUsing(const BlockVector& sol)
{
  determineExtremesOf(sol);
  
  const int* rows = TCD_limiter<d>::rowsColsRecord.rows;
  const int* cols = TCD_limiter<d>::rowsColsRecord.cols;
  
  const double* diffMatEntries = TCD_limiter<d>::diffMat.GetEntries();
  const double* origStiffMatEntries = origStiffMat.GetEntries();
  double* fluxes = TCD_limiter<d>::fluxes.data();
  
  double val1, val2;
  double diffSol_ij, diffSol_ji;
  double diffProduct, stiffProduct;
  
  for(int row : TCD_limiter<d>::rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      diffProduct = diffMatEntries[i]*( sol[col] + sol[row] );
      stiffProduct = origStiffMatEntries[i]*( sol[col] - sol[row] );
      diffSol_ij = diffProduct + stiffProduct;
       
      stiffProduct = origStiffMat.get(col, row)*( sol[row] - sol[col] );
      diffSol_ji = diffProduct + stiffProduct;
        
      if(fluxes[i] > 0.0)
      {
        val1 = diffSol_ij - 2.0*diffMatEntries[i]*solMax[row];
        val2 = 2.0*diffMatEntries[i]*solMin[col] - diffSol_ji;
       
        fluxes[i] = minFluxMin( fluxes[i], val1, val2 );
      }
      else
      {
        if(fluxes[i] < 0.0)
        {
          val1 = diffSol_ij - 2.0*diffMatEntries[i]*solMin[row];
          val2 = 2.0*diffMatEntries[i]*solMax[col] - diffSol_ji;
          
          fluxes[i] = maxFluxMax( fluxes[i], val1, val2 );
        }
      }
    }
  }
}


template <int d>
AFC_TCD<d>::AFC_TCD(Params& params)
: projectionOntoAdmissibleValues(params.projectionOntoAdmissibleValues),
  rowsColsRecord(*params.RowsColsRecordPtr),
  bcSetter(*params.tcdBcSetterPtr),
  lumpMatDiagonal(params.massMatPtr->get_n_rows(), 0.0),
  diffMat(*params.stiffMatPtr),
  stiffMat(*params.stiffMatPtr),
  massMat(*params.massMatPtr),
  solOld(*params.solOldPtr),
  rhs(*params.rhsPtr),
  solDot(*params.solOldPtr),
  fluxes(params.stiffMatPtr->get_n_entries(), 0.0),
  limiter(params.limiter),
  timeStep(params.timeStep),
  timeSteppingScheme(params.timeSteppingScheme),
  complexity(params.complexity),
  admissibleMinimum(params.admissibleMinimum),
  admissibleMaximum(params.admissibleMaximum)
{
  params.lumpMatDiagonalPtr = &lumpMatDiagonal;
  params.diffMatPtr = &diffMat;
  params.fluxesPtr = &fluxes;
  
  switch(params.limiter)
  {
    case ZALESAK:
      limiterPtr = new TCD_limiter_zalesak<d>(params);
      break;
      
    case MONOLITHIC:
      limiterPtr = new TCD_limiter_monolithic<d>(params);
      break;
  };
  
  calculateLumpMatDiagonal();
}


template <int d>
void AFC_TCD<d>::projectOntoAdmissibleValues(BlockVector& solution)
{
  for(int row : rowsColsRecord.myActiveMasterRows)
  {
    if(solution[row] > admissibleMaximum)
      solution[row] = admissibleMaximum;
    else if(solution[row] < admissibleMinimum)
      solution[row] = admissibleMinimum;
  }
}


template <int d>
void AFC_TCD<d>::calculateDiffMatFrom(const FEMatrix& someStiffMat)
{
  double* diffMatEntries = diffMat.GetEntries();
  const double* stiffMatEntries = someStiffMat.GetEntries();
  const int* rows = rowsColsRecord.rows;
  const int* cols = rowsColsRecord.cols;
  double stiffMat_ij, stiffMat_ji, minimum;
  
  for(int row : rowsColsRecord.myMasterRows)
  {
    double rowSum = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if(row != col)
      {
        stiffMat_ij = stiffMatEntries[i];
        stiffMat_ji = someStiffMat.get(col, row);
        minimum = std::min({ - stiffMat_ij , 0.0 , - stiffMat_ji });
        
        diffMatEntries[i] = minimum;
        rowSum += minimum;
      }
    }
    
    int diagonal = rowsColsRecord.indicesOfDiagEntries[row];
    
    diffMatEntries[diagonal] = - rowSum;
  }
}


template <int d>
void AFC_TCD<d>::calculateFluxes(const BlockVector& solDot, const BlockVector& sol)
{
  const int* rows = rowsColsRecord.rows;
  const int* cols = rowsColsRecord.cols;
  const double* massMatEntries = massMat.GetEntries();
  const double* diffMatEntries = diffMat.GetEntries();
  
  double massProduct, diffProduct;
  
  for(int row : rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      massProduct = massMatEntries[i]*( solDot[row] - solDot[col] );
      diffProduct = diffMatEntries[i]*( sol[row] - sol[col] );
      
      fluxes[i] = massProduct - diffProduct;
    }
  }
}


template <int d>
void AFC_TCD<d>::calculateLumpMatDiagonal()
{
  const double* massMatEntries = massMat.GetEntries();
  const int* rows = rowsColsRecord.rows;
  
  for(int row : rowsColsRecord.myMasterRows)
    for(int i = rows[row]; i < rows[row+1]; i++)
      lumpMatDiagonal[row] += massMatEntries[i];
}


template <int d>
void AFC_TCD<d>::calculateSolDot(const BlockVector& sol,
                                 const BlockVector& rhs)
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  const int* cols = AFC_TCD<d>::rowsColsRecord.cols;
  const double* stiffMatEntries;
  
  if(AFC_TCD<d>::limiter == MONOLITHIC && AFC_TCD<d>::complexity == NONLINEAR)
  {
    TCD_limiter_monolithic<d>* limiterPtr;
    limiterPtr = static_cast<TCD_limiter_monolithic<d>*>( AFC_TCD<d>::limiterPtr);
    stiffMatEntries = limiterPtr->getDiffPartMatOrStiffWithoutDiffMat().GetEntries();
  }
  else // zalesak or EXPLICIT MONOLITHIC
    stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    double sumStiffMatSol = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      sumStiffMatSol += stiffMatEntries[i]*sol[col];
    }
    
    solDot[row] = ( rhs[row] - sumStiffMatSol ) / AFC_TCD<d>::lumpMatDiagonal[row];
  }
  
  AFC_TCD<d>::bcSetter.setBTimeDerivativesIn(solDot);
  
  #ifdef _MPI
  AFC_TCD<d>::updateConsistency(solDot.get_entries(), 2);
  #endif
}


template <int d>
AFC_TCD_explicit<d>::AFC_TCD_explicit(Params& params)
: AFC_TCD<d>(params),
  doAssembling(false),
  solutionTime(0.0),
  solRK(*params.solOldPtr),
  rungeKuttaCoeff(*params.solOldPtr)
{}


template <int d>
BlockVector& AFC_TCD_explicit<d>::getMonolithicRungeKuttaCoeff()
{
  TCD_limiter_monolithic<d>* limiterPtr;
  limiterPtr = static_cast<TCD_limiter_monolithic<d>*>( AFC_TCD<d>::limiterPtr);
  
  limiterPtr->setOrigStiffMat(AFC_TCD<d>::stiffMat);
  
  AFC_TCD<d>::calculateDiffMatFrom(AFC_TCD<d>::stiffMat);
  AFC_TCD<d>::calculateStiffMatFromDiffMat(AFC_TCD<d>::stiffMat);
  
  AFC_TCD<d>::bcSetter.setBoundaryTimeDerivativeVectorForTime(TDatabase::TimeDB->CURRENTTIME);
  AFC_TCD<d>::calculateSolDot(solRK, AFC_TCD<d>::rhs);
  AFC_TCD<d>::calculateFluxes(AFC_TCD<d>::solDot, solRK);
  limiterPtr->limitFluxesUsing(solRK);
  
  FEMatrix& diffPartMat = limiterPtr->getDiffPartMatOrStiffWithoutDiffMat();
  addDiffContributionIntoDiffAndStiffMat(diffPartMat);
  
  assembleRightHandSideRows();
  checkCFLCondition();
  
  return rungeKuttaCoeff;
}


template <int d>
void AFC_TCD_explicit<d>::addDiffContributionIntoDiffAndStiffMat(const FEMatrix& diffPartMat)
{
  double* stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  const double* diffMatEntries = AFC_TCD<d>::diffMat.GetEntries();
  const double* diffPartMatEntries = diffPartMat.GetEntries();
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  
  for(int row = 0; row < AFC_TCD<d>::stiffMat.get_n_rows(); row++)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      stiffMatEntries[i] += diffPartMatEntries[i] - diffMatEntries[i];
    }
  }
  
  AFC_TCD<d>::calculateDiffMatFrom(AFC_TCD<d>::stiffMat);
  AFC_TCD<d>::calculateStiffMatFromDiffMat(AFC_TCD<d>::stiffMat);
}


template <int d>
void AFC_TCD_explicit<d>::assembleRightHandSideRows()
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  const int* cols = AFC_TCD<d>::rowsColsRecord.cols;
  const double* stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    double sumOverNeighbours = 0.0;
    double stiffMatProduct = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      stiffMatProduct += stiffMatEntries[i]*solRK[col];
      sumOverNeighbours += AFC_TCD<d>::fluxes[i];
    }
    
    rungeKuttaCoeff[row] = AFC_TCD<d>::rhs[row] + sumOverNeighbours
                           - stiffMatProduct;
    rungeKuttaCoeff[row] *= AFC_TCD<d>::timeStep/AFC_TCD<d>::lumpMatDiagonal[row];
  }
}


template <int d>
void AFC_TCD_explicit<d>::checkCFLCondition()
{
  double product;
  double maxProduct = 0.0;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myMasterRows)
  {
    product = AFC_TCD<d>::timeStep*AFC_TCD<d>::diffMat.get(row, row);
    product /= AFC_TCD<d>::lumpMatDiagonal[row];
    
    if(product > 0.5)
    {
      Output::print("WARNING: CFL condition violated for row ",
                    row, ". Product is ", product);
      
      if(product > maxProduct)
        maxProduct = product;
    }
  }
  
  if(maxProduct > 0.0)
  {
    double upperBoundForTimeStep = 0.5 / ( maxProduct / AFC_TCD<d>::timeStep );
    Output::print("HINT: Time step needs to be less than or equal to ", upperBoundForTimeStep, " to meet CFL condition.");
  }
}


template <int d>
BlockVector& AFC_TCD_explicit<d>::getZalesakRungeKuttaCoeff()
{
  AFC_TCD<d>::calculateDiffMatFrom(AFC_TCD<d>::stiffMat);
  AFC_TCD<d>::calculateStiffMatFromDiffMat(AFC_TCD<d>::stiffMat);
  
  AFC_TCD<d>::bcSetter.setBoundaryTimeDerivativeVectorForTime(TDatabase::TimeDB->CURRENTTIME);
  AFC_TCD<d>::calculateSolDot(solRK, AFC_TCD<d>::rhs);
  
  TCD_limiter_zalesak<d>* limiterPtr;
  limiterPtr = static_cast<TCD_limiter_zalesak<d>*>(AFC_TCD<d>::limiterPtr);
  
  AFC_TCD<d>::calculateFluxes(AFC_TCD<d>::solDot, solRK);
  limiterPtr->limitFluxesExplicit(AFC_TCD<d>::solDot, solRK);
  
  assembleRightHandSideRows();
  checkCFLCondition();
  
  return rungeKuttaCoeff;
}


template <int d>
AFC_TCD_implicit<d>::AFC_TCD_implicit(Params& params)
: AFC_TCD<d>(params),
  solPrevNonlin(*params.solOldPtr),
  rhsOld(*params.rhsOldPtr),
  rhsFinal(*params.rhsFinalPtr),
  termWithMatrixForRhsFinal(*params.rhsPtr),
  fluxOldPartOfRhsFinal(*params.rhsPtr),
  iterationNumber(*params.interationNumberPtr),
  theta(params.theta),
  oneMinusTheta(1.0 - params.theta)
{}


template <int d>
void AFC_TCD_implicit<d>::calculateMatAndRhs()
{
  switch( AFC_TCD<d>::limiter)
  {
    case ZALESAK:
      calculateZalesakMatAndRhs();
      break;
      
    case MONOLITHIC:
      calculateMonolithicMatAndRhs();
      break;
  };
}



template <int d>
void AFC_TCD_implicit<d>::calculateZalesakMatAndRhs()
{
  TCD_limiter_zalesak<d>* limiterPtr;
  limiterPtr = static_cast<TCD_limiter_zalesak<d>*>(AFC_TCD<d>::limiterPtr);
  
  if( notInNonlinearLoopAfterFirstIteration() )
  {
    AFC_TCD<d>::calculateDiffMatFrom(AFC_TCD<d>::stiffMat);
    AFC_TCD<d>::calculateStiffMatFromDiffMat(AFC_TCD<d>::stiffMat);
    
    limiterPtr->calculateSolInterm();
    
    if(AFC_TCD<d>::complexity == LINEAR)
      setSolPrevNonlinForLinearScheme(limiterPtr->getSolInterm());
  }
  
  calculateImplicitZalesakFluxes();
  limiterPtr->limitFluxes();
  
  if(AFC_TCD<d>::timeSteppingScheme != BACKWARD_EULER)
    checkCFLCondition();
  
  calculateZalesakRhsFinal();
  
  if( notInNonlinearLoopAfterFirstIteration() )
    calculateStiffMatFinal();
}


template <int d>
bool AFC_TCD_implicit<d>::notInNonlinearLoopAfterFirstIteration()
{
  if(AFC_TCD<d>::complexity == NONLINEAR)
    if(iterationNumber > 0)
      return false;
  
  return true;
}


template <int d>
void AFC_TCD_implicit<d>::calculateImplicitZalesakFluxes()
{
  const double* massMatEntries = AFC_TCD<d>::massMat.GetEntries();
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  const int* cols = AFC_TCD<d>::rowsColsRecord.cols;
  const double* diffMatEntries = AFC_TCD<d>::diffMat.GetEntries();
  const double* solOldEntries = AFC_TCD<d>::solOld.get_entries();
  
  double massProduct, diffProduct;
  double solOldDiff, solPrevNonlinDiff;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      if( row == col )
        continue;
      
      solOldDiff = solOldEntries[row] - solOldEntries[col];
      solPrevNonlinDiff = solPrevNonlin[row] - solPrevNonlin[col];
      
      massProduct = massMatEntries[i]*( solPrevNonlinDiff - solOldDiff ) / AFC_TCD<d>::timeStep;
      
      diffProduct = diffMatEntries[i] * theta * solPrevNonlinDiff
                    + diffMatEntries[i] * oneMinusTheta *solOldDiff;
      
      AFC_TCD<d>::fluxes[i] = massProduct - diffProduct;
    }
  }
}


template <int d>
void AFC_TCD_implicit<d>::checkCFLCondition()
{
  double* stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  double cflMin = AFC_TCD<d>::timeStep*oneMinusTheta;
  double cfl;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    int i = AFC_TCD<d>::rowsColsRecord.indicesOfDiagEntries[row];
    
    cfl = AFC_TCD<d>::lumpMatDiagonal[row]/stiffMatEntries[i];

     if(cfl < cflMin)
        cflMin = cfl;
  }
  
  if(cflMin < AFC_TCD<d>::timeStep*oneMinusTheta)
  {
    double upperBoundForTimeStep = cflMin / oneMinusTheta;
    
    Output::print("WARNING: CFL condition violated");
    Output::print("HINT: Time step needs to be less than or equal to ",
                  upperBoundForTimeStep, " to meet CFL condition.");
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateZalesakRhsFinal()
{
  if(AFC_TCD<d>::limiter == ZALESAK)
  {
    TCD_limiter_zalesak<d>* zalesakLimiterPtr;
    zalesakLimiterPtr = static_cast<TCD_limiter_zalesak<d>*>(AFC_TCD<d>::limiterPtr);
        calculateZalesakRhsFinalActiveRowsWith (zalesakLimiterPtr->getSolInterm());
  }
  
    calculateZalesakRhsFinalNonActiveRows();
}


template <int d>
void AFC_TCD_implicit<d>::calculateZalesakRhsFinalActiveRowsWith (const BlockVector& solInterm)
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  double* rhsFinalEntries = rhsFinal.get_entries();
  const double* rhsEntries = AFC_TCD<d>::rhs.get_entries();
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    double finalFlux = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
      finalFlux += AFC_TCD<d>::fluxes[i];
    
    rhsFinalEntries[row] = solInterm[row] * AFC_TCD<d>::lumpMatDiagonal[row]
                           + theta * AFC_TCD<d>::timeStep * rhsEntries[row]
                           + AFC_TCD<d>::timeStep * finalFlux;
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateZalesakRhsFinalNonActiveRows()
{
  AFC_TCD<d>::bcSetter.setBoundaryValuesVectorForTime(TDatabase::TimeDB->CURRENTTIME);
  AFC_TCD<d>::bcSetter.setBCsIn(rhsFinal);
}


template <int d>
void AFC_TCD_implicit<d>::calculateStiffMatFinal()
{
  calculateStiffMatFinalActiveRows();
  calculateStiffMatFinalNonActiveRows();
  AFC_TCD<d>::stiffMat.correct_hanging_rows();
}


template <int d>
void AFC_TCD_implicit<d>::calculateStiffMatFinalActiveRows()
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  double* stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
      stiffMatEntries[i] = theta*AFC_TCD<d>::timeStep*stiffMatEntries[i];
    
    int diagonal = AFC_TCD<d>::rowsColsRecord.indicesOfDiagEntries[row];
    stiffMatEntries[diagonal] += AFC_TCD<d>::lumpMatDiagonal[row];
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateStiffMatFinalNonActiveRows()
{
  AFC_TCD<d>::stiffMat.resetNonActive();
  AFC_TCD<d>::stiffMat.set_dirichlet_diagonals();
}


template <int d>
void AFC_TCD_implicit<d>::calculateMonolithicMatAndRhs()
{
  TCD_limiter_monolithic<d>* limiterPtr;
  limiterPtr = static_cast<TCD_limiter_monolithic<d>*>(AFC_TCD<d>::limiterPtr);
  
  if( notInNonlinearLoopAfterFirstIteration() )
  {
    limiterPtr->setOrigStiffMat( AFC_TCD<d>::stiffMat);
    FEMatrix& diffPartMatOrStiffWithoutDiffMat = limiterPtr->getDiffPartMatOrStiffWithoutDiffMat();
    
    switchMats(AFC_TCD<d>::stiffMat, diffPartMatOrStiffWithoutDiffMat);
    
    AFC_TCD<d>::stiffMat += &diffPartMatOrStiffWithoutDiffMat;
    AFC_TCD<d>::calculateDiffMatFrom(AFC_TCD<d>::stiffMat);
    AFC_TCD<d>::calculateStiffMatFromDiffMat(AFC_TCD<d>::stiffMat);
    
    checkCFLCondition();
    
    calculateTermWithMatrixForFinalRhsFinal();
    
    calculateStiffMatFinal();
        calculateZalesakRhsFinalNonActiveRows();
    
    AFC_TCD<d>::calculateDiffMatFrom(diffPartMatOrStiffWithoutDiffMat);
    AFC_TCD<d>::calculateStiffMatFromDiffMat(diffPartMatOrStiffWithoutDiffMat);
    
    AFC_TCD<d>::bcSetter.setBoundaryTimeDerivativeVectorForTime(TDatabase::TimeDB->CURRENTTIME - AFC_TCD<d>::timeStep);
    
    AFC_TCD<d>::calculateSolDot(AFC_TCD<d>::solOld, rhsOld);
    AFC_TCD<d>::calculateFluxes(AFC_TCD<d>::solDot, AFC_TCD<d>::solOld);
    limiterPtr->limitFluxesUsing(AFC_TCD<d>::solOld);
    
    calculateFluxesOldContributionToRhsFinal();
    
    AFC_TCD<d>::bcSetter.setBoundaryTimeDerivativeVectorForTime(TDatabase::TimeDB->CURRENTTIME);
  }
  
  AFC_TCD<d>::calculateSolDot(solPrevNonlin, AFC_TCD<d>::rhs);
  AFC_TCD<d>::calculateFluxes(AFC_TCD<d>::solDot, solPrevNonlin);
  limiterPtr->limitFluxesUsing(solPrevNonlin);
  
    calculateMonolithicRhsFinalActiveRows();
}


template <int d>
void AFC_TCD_implicit<d>::switchMats(FEMatrix& mat1, FEMatrix& mat2)
{
  double* mat1Entries = mat1.GetEntries();
  double* mat2Entries = mat2.GetEntries();
  
  for(int i = 0; i < mat1.get_n_entries(); i++)
  {
    double temp = mat1Entries[i];
    
    mat1Entries[i] = mat2Entries[i];
    mat2Entries[i] = temp;
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateTermWithMatrixForFinalRhsFinal()
{
  const double* stiffMatEntries = AFC_TCD<d>::stiffMat.GetEntries();
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  const int* cols = AFC_TCD<d>::rowsColsRecord.cols;
  const double* solOldEntries = AFC_TCD<d>::solOld.get_entries();
  
  double timeStep = AFC_TCD<d>::timeStep;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    double matrixRowVectorProduct = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      matrixRowVectorProduct += stiffMatEntries[i]*solOldEntries[col];
    }
    
    termWithMatrixForRhsFinal[row] = AFC_TCD<d>::lumpMatDiagonal[row] * solOldEntries[row]
                                   - oneMinusTheta * timeStep * matrixRowVectorProduct;
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateFluxesOldContributionToRhsFinal()
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
     fluxOldPartOfRhsFinal[row] = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
       fluxOldPartOfRhsFinal[row] += oneMinusTheta * AFC_TCD<d>::timeStep * AFC_TCD<d>::fluxes[i];
    }
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateImplicitMonolithicFluxes(const BlockVector& sol)
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  const int* cols = AFC_TCD<d>::rowsColsRecord.cols;
  const double* massMatEntries = AFC_TCD<d>::massMat.GetEntries();
  const double* diffMatEntries = AFC_TCD<d>::diffMat.GetEntries();
  
  double massProduct, diffProduct;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myMasterRows)
  {
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      int col = cols[i];
      
      massProduct = massMatEntries[i]*( AFC_TCD<d>::solDot[row] - AFC_TCD<d>::solDot[col] );
      diffProduct = diffMatEntries[i]*( sol[row] - sol[col] );
      
      AFC_TCD<d>::fluxes[i] = massProduct - diffProduct;
    }
  }
}


template <int d>
void AFC_TCD_implicit<d>::calculateMonolithicRhsFinalActiveRows()
{
  const int* rows = AFC_TCD<d>::rowsColsRecord.rows;
  double timeStep = AFC_TCD<d>::timeStep;
  
  for(int row : AFC_TCD<d>::rowsColsRecord.myActiveMasterRows)
  {
    double finalFlux = 0.0;
    
    for(int i = rows[row]; i < rows[row+1]; i++)
    {
      finalFlux += theta * timeStep * AFC_TCD<d>::fluxes[i];
    }
    
    double source = theta * timeStep * AFC_TCD<d>::rhs[row]
                  + oneMinusTheta * timeStep * rhsOld[row];
           
    rhsFinal[row] = termWithMatrixForRhsFinal[row] + source + finalFlux + fluxOldPartOfRhsFinal[row];
  }
}



#ifdef __3D__
template class RowsColsRecord<3>;
template class TCD_BC_setter<3>;
template class TCD_limiter<3>;
template class TCD_limiter_zalesak<3>;
template class TCD_limiter_monolithic<3>;
template class AFC_TCD<3>;
template class AFC_TCD_implicit<3>;
template class AFC_TCD_explicit<3>;
#else
template class RowsColsRecord<2>;
template class TCD_BC_setter<2>;
template class TCD_limiter<2>;
template class TCD_limiter_zalesak<2>;
template class TCD_limiter_monolithic<2>;
template class AFC_TCD<2>;
template class AFC_TCD_implicit<2>;
template class AFC_TCD_explicit<2>;
#endif
