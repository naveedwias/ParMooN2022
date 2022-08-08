#include <Database.h>
#include <FEMatrix.h>
#include <BlockVector.h>


#include <Solver.h>


#include <algorithm>
#include "templateNames.h"
#ifdef __2D__
#include "FEFunction2D.h"
#else
#include "FEFunction3D.h"
#endif
#include <string.h>

#ifdef _MPI
#include <mpi.h>
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#include "MooNMD_Io.h"


enum AFC_TCD_limiter
{
  ZALESAK,
  MONOLITHIC
};

enum AFC_TCD_Complexity
{
  LINEAR = 0,
  NONLINEAR = 1,
  EXPLICIT = 2
};

enum AFC_TCD_TimeSteppingScheme
{
  CRANK_NICOLSON = 0,
  BACKWARD_EULER = 1,
  RUNGE_KUTTA_HEUN = 2
};

enum AFC_TCD_Prelimiter
{
  NONE = 0,
  MIN_MOD = 1,
  GRAD_DIRECTION = 2,
  BOTH = 3
};

template<int d>
struct RowsColsRecord;

template<int d>
class TCD_BC_setter;


/** @brief Object that is used for passing necessary
 * information to constructor of AFC_TCD class. */
template<int d>
struct AFC_TCD_params
{
  using FEFunction = typename Template_names<d>::FEFunction;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction
    = typename Template_names<d>::BoundaryConditionFunction;
  using FESpace = typename Template_names<d>::FESpace;
  
  FEMatrix* stiffMatPtr;
  const FEMatrix* massMatPtr;
  const BlockVector* rhsPtr;
  const BlockVector* rhsOldPtr;
  BlockVector* rhsFinalPtr;
  const BlockVector* solOldPtr;
  double oldTime;
  double timeStep;
  double theta;
  const unsigned int* interationNumberPtr;
  const BoundaryValuesFunction* boundaryValuesFunction;
  const BoundaryConditionFunction* boundaryConditionFunction;
  std::shared_ptr<const FESpace> solFESpacePtr;
  AFC_TCD_limiter limiter;
  AFC_TCD_Complexity complexity;
  AFC_TCD_TimeSteppingScheme timeSteppingScheme;
  AFC_TCD_Prelimiter prelimiter;
  double toleranceForZero;
  bool projectionOntoAdmissibleValues;
  double admissibleMinimum;
  double admissibleMaximum;
  const BoundaryValuesFunction* boundaryTimeDerivativeFct;
  bool solDotUsedInZalesak;
  
  std::vector<double> *fluxesPtr;
  const FEMatrix *diffMatPtr;
  const std::vector<double> *lumpMatDiagonalPtr;
  const BlockVector* solPtr;
  
  RowsColsRecord<d>* RowsColsRecordPtr = nullptr;
  TCD_BC_setter<d>* tcdBcSetterPtr = nullptr;
  
  void checkParamsConsistency(const ParameterDatabase& db)
  {
    if( db["afc_fct_scheme"].is("explicit")
        && db["afc_limiter"].is("monolithic")
        && db["time_discretization"].is("Runge_Kutta_Heun") )
      return;
    
    if( db["afc_fct_scheme"].is("explicit")
        && db["afc_limiter"].is("zalesak")
        && db["time_discretization"].is("Runge_Kutta_Heun") )
      return;
    
    if( db["afc_fct_scheme"].is("non-linear")
        && db["afc_limiter"].is("monolithic")
        && ( db["time_discretization"].is("crank_nicolson")
             || db["time_discretization"].is("backward_euler") ) )
      return;
    
    if( db["afc_fct_scheme"].is("linear")
          || db["afc_fct_scheme"].is("non-linear") )
      if( db["afc_limiter"].is("zalesak") )
        if( db["time_discretization"].is("backward_euler")
            || db["time_discretization"].is("crank_nicolson") )
          return;
    
    ErrThrow("Unsupported combination of afc_fct_scheme,"
             " afc_limiter and time_discretization.\n",
             "Following combinations are supported:\n",
             "afc_fct_scheme = explicit and afc_limiter = monolithic "
             "and time_discretization = Runge_Kutta_Heun,\n",
             "afc_fct_scheme = non-linear and afc_limiter = monolithic "
             "and time_discretization = crank_nicolson,\n",
             "afc_fct_scheme = linear/non-linear and afc_limiter = zalesak "
             "and time_discretization = backward_euler/crank_nicolson");
    
    if( db["projection_onto_admissible_interval"].is("yes") )
      if( (double)db["admissible_minimum"] >= (double)db["admissible_maximum"] )
        ErrThrow("Projection onto admissible interval with endpoints that are not valid.");
  }
};


/** @brief Records active and nonactive master rows
 * and indices of diagonal entries in schemes
 * for algebraic flux correction for TCD.*/
template<int d>
struct RowsColsRecord
{
  public:
    
    using Params = AFC_TCD_params<d>;
    
    RowsColsRecord(const Params& params);
    
    RowsColsRecord() = delete;
    RowsColsRecord(const RowsColsRecord&) = delete;
    RowsColsRecord(RowsColsRecord&&) = delete;
    RowsColsRecord& operator=(const RowsColsRecord&) = delete;
    RowsColsRecord& operator=(RowsColsRecord&&) = delete;
    
    
    std::vector<int> myMasterRows;
    std::vector<int> myActiveMasterRows;
    std::vector<int> myNonActiveMasterRows;
    std::vector<int> indicesOfDiagEntries;
    
    const int* rows;
    const int* cols;
  
  protected:
    
    const FEMatrix& massMat;
    
    void classifyRows();
    void recordMyMasterRows();
    void recordMyActiveMasterRows();
    void recordMyNonActiveMasterRows();
    
    void storeIndicesOfDiagEntries();
};


/** @brief Sets boundary values of BlockVectors
 * in algebraic flux correction for TCD.*/
template<int d>
class TCD_BC_setter
{
  public:
    
    using FEFunction = typename Template_names<d>::FEFunction;
    using BoundaryValuesFunction
      = typename Template_names<d>::BoundaryValuesFunction;
    using BoundaryConditionFunction
      = typename Template_names<d>::BoundaryConditionFunction;
    using Params = AFC_TCD_params<d>;
    
    
    TCD_BC_setter(const Params& params);
    
    TCD_BC_setter() = delete;
    TCD_BC_setter(const TCD_BC_setter&) = delete;
    TCD_BC_setter(TCD_BC_setter&&) = delete;
    TCD_BC_setter& operator=(const TCD_BC_setter&) = delete;
    TCD_BC_setter& operator=(TCD_BC_setter&&) = delete;
    
    void setBoundaryValuesVectorForTime(double time);
    void setBCsIn(BlockVector& vec)
    { vec.copy_nonactive(boundaryValuesVector); }
    
    void setBoundaryTimeDerivativeVectorForTime(double time);
    void setBTimeDerivativesIn(BlockVector& vec)
    { vec.copy_nonactive(boundaryTimeDerivativeVector); }
  
  private:
    
    BlockVector boundaryValuesVector;
    const BoundaryValuesFunction& boundaryValuesFunction;
    const BoundaryConditionFunction& boundaryConditionFunction;
    std::vector<double> boundaryValuesFEFunctionEntries;
    FEFunction boundaryValuesFEFunction;
    
    const BoundaryValuesFunction* boundaryTimeDerivativeFct;
    BlockVector boundaryTimeDerivativeVector;
    std::vector<double> boundaryTimeDerivativeFEFctEntries;
    FEFunction boundaryTimeDerivativeFEFct;
    
    double timeStep;
    
};


/** @brief Base class that carries out limiting of fluxes
 * in algebraic flux correction for TCD.*/
template<int d>
class TCD_limiter
{
  public:
    
    using Params = AFC_TCD_params<d>;
    
    TCD_limiter(const Params& params);
      
    TCD_limiter() = delete;
    TCD_limiter(const TCD_limiter&) = delete;
    TCD_limiter(TCD_limiter&&) = delete;
    TCD_limiter& operator=(const TCD_limiter&) = delete;
    TCD_limiter& operator=(TCD_limiter&&) = delete;
  
  protected:
    
    const RowsColsRecord<d>& rowsColsRecord;
    
    std::vector<double>& fluxes;
    const FEMatrix& diffMat;
    const std::vector<double>& lumpMatDiagonal;
    const BlockVector& sol;
    
    TCD_BC_setter<d>& bcSetter;
    
    #ifdef _MPI
    void updateConsistency(double* data, size_t level)
    {
      const TParFECommunicator3D& comm = diffMat.GetFESpace3D()->get_communicator();
      comm.consistency_update(data, level);
    }
    #endif
};


/** @brief Carries out Zalesak's limiting of fluxes
 * in algebraic flux correction for TCD.*/
template<int d>
class TCD_limiter_zalesak : public TCD_limiter<d>
{
  public:
    
    using Params = AFC_TCD_params<d>;
    
    TCD_limiter_zalesak(const Params& params);
      
    TCD_limiter_zalesak() = delete;
    TCD_limiter_zalesak(const TCD_limiter_zalesak&) = delete;
    TCD_limiter_zalesak(TCD_limiter_zalesak&&) = delete;
    TCD_limiter_zalesak& operator=(const TCD_limiter_zalesak&) = delete;
    TCD_limiter_zalesak& operator=(TCD_limiter_zalesak&&) = delete;
  
    void limitFluxes();
    void limitFluxesExplicit(const BlockVector& solDot, const BlockVector& sol);
    const BlockVector& getSolInterm()
    { return solInterm; }
    
    void calculateSolInterm();
    void calculateSolIntermExplicit(const BlockVector& solDot, const BlockVector& sol);
    
    void setSolutionTime(double time)
    { solutionTime = time; }
  
  protected:
    
    double solutionTime;
    
    std::vector<double> P_plus;
    std::vector<double> P_minus;
    std::vector<double> Q_plus;
    std::vector<double> Q_minus;
    std::vector<double> R_plus;
    std::vector<double> R_minus;
    double timeStep;
    
    void prelimitFluxes();
    void calculateRs();
    void updateFluxes();
  
    BlockVector solInterm;
    const BlockVector& solOld;
    const BlockVector& rhsOld;
    const FEMatrix& stiffMat;
    
    
    double toleranceForZero;
    
    static double minMod(double a, double b)
    {
      if(a*b < 0)
        return 0.0;
      else
        return std::abs(a) <  std::abs(b) ? a : b;
    }
    
    bool makesSenceToCompareSignOf(double number)
    {
      if( std::abs(number) >= toleranceForZero )
        return true;
      else
        return false;
    }
    
    const AFC_TCD_TimeSteppingScheme timeSteppingScheme;
    const AFC_TCD_Prelimiter prelimiter;
    const double oneMinusTheta;
    
    
    void calculateSolIntermNonOneTheta();
    void calculateSolIntermNonOneThetaActiveRows();
    void calculateSolIntermNonOneThetaNonActiveRows();
    void calculateSolIntermBackwardEuler()
    { solInterm = solOld; }
};


/** @brief Carries out monolithic limiting of fluxes
 * in algebraic flux correction for TCD.*/
template<int d>
class TCD_limiter_monolithic : public TCD_limiter<d>
{
  public:
    
    using Params = AFC_TCD_params<d>;
    
    TCD_limiter_monolithic(const Params& params);
      
    TCD_limiter_monolithic() = delete;
    TCD_limiter_monolithic(const TCD_limiter_monolithic&) = delete;
    TCD_limiter_monolithic(TCD_limiter_monolithic&&) = delete;
    TCD_limiter_monolithic& operator=(const TCD_limiter_monolithic&) = delete;
    TCD_limiter_monolithic& operator=(TCD_limiter_monolithic&&) = delete;
  
    void limitFluxesUsing(const BlockVector& sol);
    
    FEMatrix& getOrigStiffMat()
    { return origStiffMat; }
    
    void setOrigStiffMat(const FEMatrix& stiffMat)
    { origStiffMat = stiffMat; }
    
    FEMatrix& getDiffPartMatOrStiffWithoutDiffMat()
    { return diffPartMatOrStiffWithoutDiffMat; }
  
  protected:
    
    BlockVector solMax;
    BlockVector solMin;
    
    FEMatrix origStiffMat;
    FEMatrix diffPartMatOrStiffWithoutDiffMat;
    
    void determineExtremesOf(const BlockVector& sol);
    
    
    static double minFluxMin(double flux, double number1, double number2)
    {
      double innerMin = std::min({number1, number2});
      return std::min({flux, innerMin});
    }
    
    static double maxFluxMax(double flux, double number1, double number2)
    {
      double innerMax = std::max({number1, number2});
      return std::max({flux, innerMax});
    }
  
};

/** @brief Base class that carries out algebraic flux correction
 * for TCD. */
template<int d>
class AFC_TCD
{
  using Limiter = TCD_limiter<d>;
  using Params = AFC_TCD_params<d>;
  using MonolithicLimiter = TCD_limiter_monolithic<d>;
  
  public:
    
    AFC_TCD(Params& params);
      
    AFC_TCD() = delete;
    AFC_TCD(const AFC_TCD&) = delete;
    AFC_TCD( AFC_TCD&&) = delete;
    AFC_TCD& operator=(const AFC_TCD&) = delete;
    AFC_TCD& operator=( AFC_TCD&&) = delete;
  
    virtual ~AFC_TCD()
    { delete limiterPtr; }
    
    FEMatrix& getDiffPartMat()
    {
      MonolithicLimiter* monolithicLimiterPtr;
      monolithicLimiterPtr = static_cast<MonolithicLimiter*>(limiterPtr);
      return monolithicLimiterPtr->getDiffPartMatOrStiffWithoutDiffMat();
    }
    
    const bool projectionOntoAdmissibleValues;
    void projectOntoAdmissibleValues(BlockVector& solution);
    
    AFC_TCD_Complexity getComplexity()
    { return complexity; }
    
    AFC_TCD_limiter getLimiter()
    { return limiter; }
    
  protected:
    
    #ifdef _MPI
    void updateConsistency(double* data, size_t level)
    {
      const TParFECommunicator3D& comm = stiffMat.GetFESpace3D()->get_communicator();
      comm.consistency_update(data, level);
    }
    #endif
    
    const RowsColsRecord<d>& rowsColsRecord;
    TCD_BC_setter<d>& bcSetter;
    
    std::vector<double> lumpMatDiagonal;
    FEMatrix diffMat;
    
    FEMatrix& stiffMat;
    const FEMatrix& massMat;
    const BlockVector& solOld;
    const BlockVector& rhs;
    BlockVector solDot;
    
    std::vector<double> fluxes;
    
    AFC_TCD_limiter limiter;
    Limiter* limiterPtr;
    
    double timeStep;
    const AFC_TCD_TimeSteppingScheme timeSteppingScheme;
    AFC_TCD_Complexity complexity;
    
    void calculateDiffMatFrom(const FEMatrix& someStiffMat);
    void calculateStiffMatFromDiffMat(FEMatrix& someStiffMat)
    { someStiffMat += &diffMat; }
    
    void calculateSolDot(const BlockVector& sol, const BlockVector& rhs);
    
    void calculateFluxes(const BlockVector& solDot, const BlockVector& sol);
    
    virtual void checkCFLCondition() = 0;
    
    void calculateLumpMatDiagonal();
    
    double admissibleMinimum;
    double admissibleMaximum;
};


/** @brief Carries out algebraic flux correction
 * for TCD using an explicit scheme.*/
template<int d>
class AFC_TCD_explicit : public AFC_TCD<d>
{
  using Params = AFC_TCD_params<d>;

  public:
  
    AFC_TCD_explicit(Params& params);
    
    AFC_TCD_explicit() = delete;
    AFC_TCD_explicit(const AFC_TCD_explicit&) = delete;
    AFC_TCD_explicit(AFC_TCD_explicit&&) = delete;
    AFC_TCD_explicit& operator=(const AFC_TCD_explicit&) = delete;
    AFC_TCD_explicit& operator=(AFC_TCD_explicit&&) = delete;
  
    BlockVector& getRungeKuttaCoeff()
    {
      if( AFC_TCD<d>::limiter == ZALESAK)
        return getZalesakRungeKuttaCoeff();
      else// MONOLITHIC
        return getMonolithicRungeKuttaCoeff();
    }
    
    bool doAssembling;
    
    void readSolForRungeKutta(const BlockVector& passedSolution)
    { solRK = passedSolution; }
    
    void readTimeStep(double passedTimeStep)
    { 
        AFC_TCD<d>::timeStep = passedTimeStep; }
    
    void setZalesakSolutionTime(double time)
    { 
      solutionTime = time; 
      
      if( AFC_TCD<d>::limiter == ZALESAK)
      {
        TCD_limiter_zalesak<d>* limiterPtr;
        limiterPtr = static_cast<TCD_limiter_zalesak<d>*>( AFC_TCD<d>::limiterPtr);
        limiterPtr->setSolutionTime(time);
      }
    }
    
  protected:
    
    BlockVector& getMonolithicRungeKuttaCoeff();
    BlockVector& getZalesakRungeKuttaCoeff();
    
    double solutionTime;
    
    BlockVector solRK;
    BlockVector rungeKuttaCoeff;
    
    void addDiffContributionIntoDiffAndStiffMat(const FEMatrix& diffPartMat);
    void assembleRightHandSideRows();
    
    void checkCFLCondition();
};


/** @brief Carries out algebraic flux correction
 * for TCD using an implicit scheme.*/
template<int d>
class AFC_TCD_implicit : public AFC_TCD<d>
{
  using Params = AFC_TCD_params<d>;
  
  public:
    
    AFC_TCD_implicit(Params& params);
    
    AFC_TCD_implicit() = delete;
    AFC_TCD_implicit(const AFC_TCD_implicit&) = delete;
    AFC_TCD_implicit(AFC_TCD_implicit&&) = delete;
    AFC_TCD_implicit& operator=(const AFC_TCD_implicit&) = delete;
    AFC_TCD_implicit& operator=(AFC_TCD_implicit&&) = delete;
  
    void readSolPrevNonlin(const BlockVector& passedSolution)
    { solPrevNonlin = passedSolution; }
    
    void calculateMatAndRhs();
  
  protected:
    
    void calculateZalesakMatAndRhs();
    void calculateMonolithicMatAndRhs();
    
    
    BlockVector solPrevNonlin;
    
    const BlockVector& rhsOld;
    BlockVector& rhsFinal;
    
    BlockVector termWithMatrixForRhsFinal;// needed only for monolithic limiter
    BlockVector fluxOldPartOfRhsFinal;// needed only for monolithic limiter
    
    const unsigned int& iterationNumber;
    
    const double theta;
    const double oneMinusTheta;
    
    bool notInNonlinearLoopAfterFirstIteration();
    void setSolPrevNonlinForLinearScheme(const BlockVector& solInterm)
    {
      for(unsigned int row = 0; row < solPrevNonlin.length(); row++)
        solPrevNonlin[row] = 2.0 * solInterm[row] - AFC_TCD<d>::solOld[row];
    }
    
    void calculateImplicitZalesakFluxes();
    void calculateImplicitMonolithicFluxes(const BlockVector& sol);
    void checkCFLCondition();
    void calculateZalesakRhsFinal();
    void calculateZalesakRhsFinalActiveRowsWith (const BlockVector& solInterm);
    void calculateZalesakRhsFinalNonActiveRows();
    
    void calculateStiffMatFinal();
    void calculateStiffMatFinalActiveRows();
    void calculateStiffMatFinalNonActiveRows();
    
    void switchMats(FEMatrix& mat1, FEMatrix& mat2);
    
    void calculateTermWithMatrixForFinalRhsFinal();
    void calculateFluxesOldContributionToRhsFinal();
    
    void calculateMonolithicRhsFinalActiveRows();
    
};


ParameterDatabase default_afc_tcd_database();

