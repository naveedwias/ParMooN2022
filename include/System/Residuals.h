/**
 * @file Residuals.h
 * Contains a small struct which contains residuals
 * as appearing in approximate solutions to the Navier--Stokes Equations.
 *
 * @date 2016/04/01
 * @author Clemens Bartsch (moved it here from NSE2D)
 */

#ifndef INCLUDE_SYSTEM_RESIDUALS_H_
#define INCLUDE_SYSTEM_RESIDUALS_H_

#include <ostream>

/**
 * @brief a simple struct storing one set of residuals
 *
 * The full residual is the \f$\ell^2\f$-norm of the vector \f$Ax-b\f$
 * where \f$A\f$ is the current matrix and \f$b\f$ the right hand side. It
 * is composed of two parts, the momentum and the mass.
 *
 * If not default constructed it holds
 *     fullResidual*fullResidual = momentumResidual*momentumResidual
 *                                 +massResidual*massResidual
 */
struct Residuals
{
    /// @brief the impulse residual
    double momentumResidual;
    /// @brief the mass residual
    double massResidual;

    /// @brief the full residual
    double fullResidual;

    double momentumResidualMax;
    double massResidualMax;
    double fullResidualMax;

    ///@brief standard constructor, initialize with large numbers (1e10).
    Residuals();

    /// @brief constructor given the \e square of the momentum and mass
    /// residuals
    Residuals(double moR2, double maR2, double moRM, double maRM);

    /// @brief cast to double.
    /// This is needed for this struct to be used in the `LoopInfo` template 
    /// class.
    operator double() { return fullResidual; }
};

/// @brief Write out the three numbers to a stream.
std::ostream& operator<<(std::ostream& s, const Residuals& n);



#endif /* INCLUDE_SYSTEM_RESIDUALS_H_ */
