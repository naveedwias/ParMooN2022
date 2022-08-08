#ifndef NSEERRORESTIMATOR2D_H
#define NSEERRORESTIMATOR2D_H

#include "templateNames.h"
#include "ErrorEstimator.h"
#ifdef __2D__
#include "Example_NSE2D.h"
#include "FEFunction2D.h"
#else
#include "Example_NSE3D.h"
#include "FEFunction3D.h"
#endif

class TAuxParam2D; // todo: remove this

enum class NSE2DErrorEstimatorType {
    // 0 - gradient indicator
    gradient_indicator = 0,
    // 1 - residual estimator h1
    residual_estimator_h1,
    // 2 - residual estimator l2
    residual_estimator_l2,
    // 3 - residual estimator energy norm quasi robust
    residual_estimator_energy_quasi_robust,
    // 4 - gradient recovery
    gradient_recovery,
    // 5 - implicit estimator neumann
    implicit_estimator_neumann
};

constexpr int N_NSE2D_ESTIMATOR_TYPES = 6;

// enable output of NSE error estimator type
std::ostream & operator<<(std::ostream &os, NSE2DErrorEstimatorType &type);

template <int d>
class NSEErrorEstimator : public ErrorEstimator<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_NSE = typename Template_names<d>::Example_NSE;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;


    NSEErrorEstimator(const ParameterDatabase& param_db);

    void estimate(const Example_NSE & ex, const FEVectFunct & fe_function_u,
                  const FEFunction & fe_function2D_p, TAuxParam2D &Aux);
    
    virtual void info() override;

  protected:
    // the selected estimator type
    NSE2DErrorEstimatorType estimatorType;
    
    // opaque struct holding relevant data of the joints
    struct EdgeData;
    // opaque struct holding relevant data of the ref joints in jump calculation
    struct EdgeRefData;
    
    // parameter indicating whether the grid is conforming or not
    int conform_grid;
    // indicating if navier stokes or not
    bool is_nse;

    void calculateEtaK(const TFEVectFunct2D &fe_function2D_u,
                       const TFEFunction2D &fe_function2D_p, const TBaseCell *cell,
                       unsigned int N_Points, unsigned int N_Points1D,
                       const double *weights, double **Derivatives, double **coeffs,
                       const Example_NSE &example, EdgeData &edgeData,
                       EdgeRefData &edgeRefData, double *estimated_local_error);
};

#endif //NSEERRORESTIMATOR2D_H
