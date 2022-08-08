//
// Created by Moritz Hoffmann on 05/07/15.
//

#ifndef CDERRORESTIMATOR2D_H
#define CDERRORESTIMATOR2D_H

#include "templateNames.h"
#include "ErrorEstimator.h"
#ifdef __2D__
#include "Example_CD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_CD3D.h"
#include "FEFunction3D.h"
#endif

class TJoint;

enum class CDErrorEstimatorType
{
  // 0 - gradient indicator
  GradientIndicator = 0,
  // 1 - H^1 estimator
  H1_ResidualEstimator,
  // 2 - L^2 estimator
  L2_ResidualEstimator,
  // 3 - energy norm + dual norm, Verf"urth 2005
  Energy_ResidualEstimatorQuasiRobust,
  // 4 - energy norm estimator without jumps
  Energy_ResidualEstimatorWithoutJumps,
  // 5 - supg estimator John/Novo, upper estimate
  SUPG_Upper,
  // 6 - supg estimator John/Novo, lower estimate
  SUPG_Lower
};

constexpr int N_CD2D_ESTIMATOR_TYPES = 7;

// enable output of CD error estimator type
std::ostream & operator<<(std::ostream &os, CDErrorEstimatorType type);

/**
 * @brief Error estimator specialization for convection-diffusion problems.
 * 
 * This (template) class essentially provides only one method 'estimate' whose
 * main purpose it is to fill the vector 'eta_K' along with the members 
 * 'maximal_local_error', 'estimated_global_error', and 'currentCollection' in
 * its base (template) class.
 * 
 * @note This is only a start, it is not thoroughly tested and currently only 
 *       works in 2D.
 */
template <int d>
class CDErrorEstimator : public ErrorEstimator<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_CD = typename Template_names<d>::Example_CD;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;

    // opaque struct holding relevant data of the joints
    struct JointData;
    // opaque struct holding relevant data of the ref joints in jump calculation
    struct JointRefData;
  
  protected:
    // needed derivatives
#ifdef __3D__
    const MultiIndex_vector derivatives{MultiIndex3D::D000, MultiIndex3D::D100,
                                        MultiIndex3D::D010, MultiIndex3D::D001,
                                        MultiIndex3D::D200, MultiIndex3D::D110,
                                        MultiIndex3D::D101, MultiIndex3D::D020,
                                        MultiIndex3D::D011, MultiIndex3D::D002};
#else
    const MultiIndex_vector derivatives{MultiIndex2D::D00, MultiIndex2D::D10,
                                        MultiIndex2D::D01, MultiIndex2D::D20,
                                        MultiIndex2D::D02};
#endif

    // internal function calculating jumps across joints at the boundary
    bool handleJump_BoundaryEdge(double *result, const Example_CD &example,
                                 const int estimatorType,
                                 const TQuadFormula& qf,
                                 const CDErrorEstimator::JointData &edgeData,
                                 const double meas, const double *coeff,
                                 double linfb, const std::vector<double> &alpha,
                                 int edgeIdx, const TJoint *joint) const;


    // calculates eta_K^2 for a single cell K
    double calculateEtaK(const TBaseCell *cell,
                         const FEFunction &fe_function,
                         const TCollection &coll, const Example_CD & example,
                         std::vector<double *> &derivativesPerQuadPoint,
                         const TQuadFormula& qf_cell,
                         std::vector<double *> &coefficientsPerQuadPoint,
                         const TQuadFormula& qf_line, const JointData &edgeData,
                         JointRefData &edgeRefData) const;

  public:
    // constructor
    CDErrorEstimator(const ParameterDatabase& param_db);

    void estimate(const Example_CD &ex, const FEFunction &fe_function);
    
    virtual void info() override;
};


#endif //CDERRORESTIMATOR2D_H
