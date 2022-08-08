//
// Created by Moritz Hoffmann on 05/07/15.
//

#ifndef CDERRORESTIMATOR2D_AFC_H
#define CDERRORESTIMATOR2D_AFC_H

#include "templateNames.h"
#include "ErrorEstimator.h"
#include "FEMatrix.h"
#ifdef __2D__
#include "Example_CD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_CD3D.h"
#include "FEFunction3D.h"
#endif
#include "Joint.h"

enum class CDErrorEstimatorType_AFC
{
  // 0- AFC 
  AFC,
  // 1- SUPG_AFC 
  SUPG_AFC
};

constexpr int N_CD2D_ESTIMATOR_TYPES_AFC = 2;

// enable output of CD error estimator type
std::ostream &operator<<(std::ostream &os, CDErrorEstimatorType_AFC type);

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
class CDErrorEstimator_AFC : public ErrorEstimator<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_CD = typename Template_names<d>::Example_CD;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;

    // opaque struct holding relevant data of the edges
    struct JointData;
    // opaque struct holding relevant data of the ref edges in jump calculation
    struct JointRefData;
    
    double estimated_afc_error;
    double estimated_edge_error;
    double estimated_element_error;
  
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
    
    // 13.6568*(Perimeter^3/(|K|^2)*(1-C_cos)) constant in eta_3
    double edge_max;

    // internal function calculating jumps across edges at the boundary
    bool handleJump_BoundaryEdge(double *result, const Example_CD &example,
                                 const TQuadFormula& qf,
                                 const CDErrorEstimator_AFC<d>::JointData &edgeData,
                                 const double *coeff,
                                 int edgeIdx, const TJoint *joint) const;
    
    // calculates eta_K^2 for a single cell K
    double calculateEtaK(const TBaseCell *cell,
                         const FEFunction &fe_function,
                         const TCollection &coll, const Example_CD & example,
                         std::vector<double *> &derivativesPerQuadPoint,
                         std::vector<double *> &initial_derivativesPerQuadPoint,
                         const TQuadFormula& qf_cell,
                         std::vector<double *> &coefficientsPerQuadPoint,
                         const TQuadFormula& qf_line, const JointData &edgeData,
                         JointRefData &edgeRefData,
                         const FEMatrix& afc_matrix_D_entries,
                         const std::vector<double>& alphas,
                         double& result_afc_error,
                         double& element_error,
                         double& edge_error ) const;

  public:
    // constructor
    CDErrorEstimator_AFC(const ParameterDatabase& param_db);

    void estimate(const Example_CD &ex, const FEFunction &fe_function,
                  const FEMatrix& afc_matrix_D_entries,
                  const std::vector<double>& alphas,
                  const FEFunction& initial_fe_function);
    
    double get_estimated_afc_error() const 
    { return estimated_afc_error; }
    
    double get_estimated_element_error() const 
    { return estimated_element_error; }
    
    double get_estimated_edge_error() const 
    { return estimated_edge_error; }
    
    virtual void info() override;
};


#endif //CDERRORESTIMATOR2D_H
