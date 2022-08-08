#ifndef __MAINUTILITIES__
#define __MAINUTILITIES__

#ifdef __2D__
  #include <FESpace2D.h>
  class TSquareMatrix2D;
#endif
#ifdef __3D__
  #include <FESpace3D.h>
  #include <FEVectFunct3D.h>
  #include <Variational_multiscale.h>
  #include <PointwiseAssemblyData.h>
  class TSquareMatrix3D;
#endif

#include <NonNewtonianViscosity.h>

#include <JointEqN.h>

#include <string.h>
#include <deque>
#include <queue>

double GetTime();
int GetMemory();
void display_mallinfo(const std::string& program_part);


/** @brief fixed size queue
 * 
 * This implements a queue with a fixed size given by the first template
 * parameter. The second parameter describes the type of objects in the 
 * queue. It is only possible to call three methods (besides the 
 * constructor). You can access the first and last entry in the queue and
 * add a new one.
 */
template<unsigned int queueSize, class T>
struct FixedSizeQueue : private std::queue<T>
{
  /// @brief construct a FixedSizeQueue with the given size
  ///
  /// this is done via a std::deque. However this is of no concern here.
  FixedSizeQueue() : std::queue<T>( std::deque<T>(queueSize) )
  {
  }
  /// @brief add an element to the queue
  ///
  /// This adds a new elements to the back and removes one from the front.
  void add(const T& t)
  {
    this->push(t);
    this->pop();
  }
  /// @brief get the first (oldest) element
  const T& front() const
  {
    return this->std::queue<T>::front();
  }
  /// @brief get the last (newest) element
  const T& back() const
  {
    return this->std::queue<T>::back();
  }
};


void SwapDoubleArray(double *doublearray, int length);
void SwapIntArray(int *intarray, int length);


#ifdef __2D__
// ====================================================================
// calculate the streamfunction from the (u1, u2) velocity field
// ====================================================================
void StreamFunction(const TFESpace2D *velo, double *u1, double *u2,
                    const TFESpace2D *stream, double *psi);

// determine L2 error
void L2Error(int N_Points, std::array<const double*, 2> xy,
              const double *AbsDetjk, const double *Weights, double hK, 
              const double *const* Der, const double *const* Exact,
              const double *const* coeffs, double *LocError);

// determine L2 and H1 error
void L2H1Errors(int N_Points, std::array<const double*, 2> xy,
                const double *AbsDetjk, const double *Weights, double hK, 
                const double *const* Der, const double *const* Exact,
                const double *const* coeffs, double *LocError);

// determine L2-error, divergence error and H1 error, 2D
void L2DivH1Errors(int N_Points, std::array<const double*, 2> xy,
                   const double *AbsDetjk, const double *Weights, double hK, 
                   const double *const* Der, const double *const* Exact,
                   const double *const* coeffs, double *LocError);

// determine L1 error, 2D
void L1Error(int N_Points, std::array<const double*, 2> xy,
             const double *AbsDetjk, const double *Weights, double hK, 
             const double *const* Der, const double *const* Exact,
             const double *const* coeffs, double *LocError);

// determine L2, H1 and SDFEM error
void SDFEMErrors(int N_Points, std::array<const double*, 2> xy,
                 const double *AbsDetjk, const double *Weights, double hK,
                 const double *const* Der, const double *const* Exact,
                 const double *const* coeffs, double *LocError);

// determine L2, H1 and SDFEM error, in (0,P6)^2
void SDFEMErrorsSmooth(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);

// determine L2, H1 and SDFEM error for smooth region in the
// example JohnMaubachTobiska1997 (x-0.5)^2+(y-0.5)^2 > r^2 
void SDFEMErrorsSmooth_JohnMaubachTobiska1997
(int N_Points, double *X, double *Y, double *AbsDetjk, 
 const double *Weights, double hK, double **Der, double **Exact,
 double **coeffs, double *LocError);

// determine errors to interpolant
// paper with Julia Novo
 void SDFEMErrorsInterpolant(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);
     
// determine deformation tensor error
void DeformationTensorError(int N_Points, double *X, double *Y, double *AbsDetjk, 
                const double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);
// compute the subgrid dissipation 
void SubGridDissipation(int N_Points, double *X, double *Y, double *AbsDetjk, 
                const double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);

// determine H1 norm
void H1Norm(int N_Points, double *X, double *Y, double *AbsDetjk, 
            const double *Weights, double hK, 
            double **Der, double **Exact,
            double **coeffs, double *LocError);

// compute the error in the divergence
void DivergenceError(int N_Points, double *X, double *Y,
		     double *AbsDetjk, const double *Weights, double hK, 
		     double **Der, double **Exact,
		     double **coeffs, double *LocError);

// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *X, double *Y, double *AbsDetjk, 
           const double *Weights, double hK, 
           double **Der, double **Exact,
           double **coeffs, double *LocError);

void DivergenceErrorGradDivOseen(int N_Points, double *X, double *Y,
         double *AbsDetjk, const double *Weights, double hK, 
         double **Der, double **Exact,
         double **coeffs, double *LocError);
         
void Parameters_Gradient_Residual(int N_Points, double *X, double *Y, double *AbsDetjk,
           const double *Weights, double hK,
           double **Der, double **Exact,
           double **coeffs, double *LocError);
#endif

double graddiv_parameterOseen(double hK, double nu, double b1, double b2);

#ifdef __3D__

// determine L2 error
void L2Error(int N_Points, std::array<const double*, 3> xyz,
                const double *AbsDetjk, 
                const double *Weights, double hK, 
                const double *const* Der, const double *const* Exact,
                const double *const* coeffs, double *LocError);

// determine L2 and H1 error
void L2H1Errors(int N_Points, std::array<const double*, 3> xyz,
                const double *AbsDetjk, 
                const double *Weights, double hK, 
                const double *const* Der, const double *const* Exact,
                const double *const* coeffs, double *LocError);

void L2H1ErrorsSmooth(int N_Points, double *X, double *Y, double *Z,
                double *AbsDetjk, 
                const double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);

// compute L2 error, L2 error of divergence, and H1 error for vector valued
// basis functions (Raviart-Thomas or Brezzi-Douglas-Marini)
void L2DivH1Errors(int N_Points, std::array<const double*, 3> xyz,
                   const double *AbsDetjk, const double *Weights, double hK,
                   const double *const* Der, const double *const* Exact,
                   const double *const* coeffs, double *LocError);

// determine L1 error
void L1Error(int N_Points, std::array<const double*, 3> xyz,
             const double *AbsDetjk, const double *Weights, double hK,
             const double *const* Der, const double *const* Exact,
             const double *const* coeffs, double *LocError);

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *X, double *Y, double *Z, 
                            double *AbsDetjk, 
                            const double *Weights, double hK, 
                            double **Der, double **Exact,
                            double **coeffs, double *LocError);
// compute the subgrid dissipation 
void SubGridDissipation(int N_Points, double *X, double *Y, double *Z, 
                        double *AbsDetjk, 
                        const double *Weights, double hK, 
                        double **Der, double **Exact,
                        double **coeffs, double *LocError);


// compute the error in the divergence
void DivergenceError(int N_Points, double *X, double *Y, double *Z,
		     double *AbsDetjk, const double *Weights, double hK, 
		     double **Der, double **Exact,
		     double **coeffs, double *LocError);

// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *X, double *Y, double *Z,
                      double *AbsDetjk, 
                      const double *Weights, double hK, 
                      double **Der, double **Exact,
                      double **coeffs, double *LocError);

void Q_criterion(TCollection *Coll,
                 TFEFunction3D *velocity1, TFEFunction3D *velocity2,
                 TFEFunction3D *velocity3, double *Qcrit);

// compute wall shear stress. Note that tau_w *has* to be a 3-component
// vector with 1:1 vertex/DOF correspondence (e.g. P1 on tetra meshes)
void ComputeWallShearStress(const double nu, const TFEVectFunct3D &u,
  TFEVectFunct3D &tau_w, const ViscositySettings& settings);
void ComputeWallShearStressFace(const double nu, const TFEVectFunct3D &u,
  TFEVectFunct3D &tau_w, const ViscositySettings& settings);

void ComputeTurbulentKineticEnergy(const TFEVectFunct3D &u,
  TFEFunction3D &tke, double scale);
void ComputeTurbulentKineticEnergyCell(const TFEVectFunct3D &u,
  TFEFunction3D &tke, double scale);

void ComputeEffectiveViscosity(const TFEVectFunct3D &u,
  TFEFunction3D &nu_eff, const ViscositySettings& settings);

double ComputeTurbulentKineticEnergySmagorinsky(const TFEVectFunct3D &u, double scale);

double ComputeTurbulentKineticEnergyRBVMS(const TFEVectFunct3D &u,
  const TFEVectFunct3D &u_m1, const TFEFunction3D &p, double nu,
  RBVMS_Settings settings);

double ComputeTurbulentKineticEnergyRBVMSTime(const TFEVectFunct3D &u,
  std::shared_ptr<PointwiseAssemblyData> persistent_data);

#endif

// ====================================================================
// auxiliary routines
// ====================================================================
void linfb(int N_Points, double **Coeffs, double ** Params,
           TBaseCell *cell);
           
void ave_l2b_quad_points(int N_Points, double **Coeffs, double **Params,
           TBaseCell *cell);

void CFPM2D_linfb(int N_Points, double **Coeffs, double ** Params,
           TBaseCell *cell);

void LInfU(int N_Points, double **Coeffs, double **Params, 
           TBaseCell *cell);


// set the exact_solution in examples to this, if it is unknown
void unknown_solution_2d(double x, double y, double *values);
void unknown_solution_3d(double x, double y, double z, double *values);
void BoundConditionNoBoundCondition(int BdComp, double t, BoundCond &cond);
void BoundConditionNoBoundCondition(int comp, double x, double y, double z, BoundCond &cond);
void BoundaryValueNoBoundaryValue(int BdComp, double Param, double &value);
void BoundaryValueHomogenous(int BdComp, double Param, double &value);
void BoundaryValueHomogenous(int comp, double x, double y, double z, double &value);
void BoundConditionVMM(int BdComp, double t, BoundCond &cond);
void BoundConditionNSE(int BdComp, double t, BoundCond &cond);
void BoundaryConditionPressSep(int i, double t, BoundCond &cond);
void BoundaryValuePressSep(int BdComp, double Param, double &value);
void BoundaryConditionPressSep3D(int comp, double x, double y, double z, BoundCond &cond);
void BoundaryValuePressSep3D(int comp, double x, double y, double z, double &value);
void BoundaryConditionNewton(int comp, double x, double y, double z, BoundCond &cond);
void BoundaryValueNewton(int comp, double x, double y, double z, double &value);
void BoundCondition_FEM_FCT(int i, double t, BoundCond &cond);
void BoundValue_FEM_FCT(int BdComp, double Param, double &value);
void BoundCondition_FEM_FCT(int comp, double x, double y, double z, BoundCond &cond);
void BoundValue_FEM_FCT(int comp, double x, double y, double z, double &value);
void BoundConditionAuxProblem(int i, double t, BoundCond &cond);
void BoundValueAuxProblem(int BdComp, double Param, double &value);
void ho_BoundCondition(int i, double t, BoundCond &cond);
void ho_BoundValue(int BdComp, double Param, double &value);

void SetPolynomialDegree();
void CheckMaximumPrinciple(TSquareMatrix *A, double *sol, int N_Active,
			   double *errors);
void SaveData(char *name, int N_Array, double **sol, int *N_Unknowns);
void ReadData(char *name, int N_Array, double **sol, int *N_Unknowns);

void SaveData(const std::string& basename, double *sol, int nDOF);
void ReadData(const std::string& filename, double *sol, int nDOF);

void ComputeVorticityDivergence(TFEFunction2D *u1, TFEFunction2D *u2,
                                const TFESpace2D *vorticity_space, 
                                double *vort,  double *div);

#endif // __MAINUTILITIES__









