#ifndef PARMOON__UTILITIES
#define PARMOON__UTILITIES

#include <string>
#include <vector>

namespace utilities
{

/// @brief return the host name of the computer ParMooN is running on
std::string get_host_name();

/// @brief return the date and time as a string
std::string get_date_and_time();

/** @brief returns minmod of a vector of numbers
 *
 * @details This function is the implementation of \f$\mathrm{minmod}(a_0, a_1,
 * \ldots, a_{N-1})\f$, \f$N\ge 1\f$, for \f$N\f$ input values where
 * \f$\mathrm{minmod}\f$ is defined by
 * \f[
 *     \mathrm{minmod}(a_0, a_1, \ldots, a_{N-1}):=
 *       \begin{cases}
 *          s \min_{i = 1, 2, \ldots, N-1}{\vert a_i \vert} &\text{ if } s :=
 *          \mathrm{sign}(a_0) = \mathrm{sign}(a_1) = \ldots =
 *          \mathrm{sign}(a_{N-1}) \\
 *          0 & \text{ else.}
 *        \end{cases}
 * \f]
 * @param[in] input_vec vector of numbers \f$a_0, a_1,
 * ..., a_{N-1}\f$
 * @return \f$\mathrm{minmod}(a_0, a_1,\ldots, a_{N-1})\f$
 */
double minmod(const std::vector< double >& input_vec);

/**
 * @brief checks if two numbers are equal up to given tolerance
 *
 * @details If the first number \f$|a|\ge 1\f$, the relative error
 *\f[
     |a-b| < \mathrm{tol} \cdot |a|
 \f]
 * between \p a and \p b is investigated. If \f$|a| < 1\f$ the absolute error
 *\f[
   |a-b| < tol
 \f]
 * is checked.
 *
 * The default value for \p tol is 1e-12.
 */
bool are_equal(double a, double b, double tol = 1e-12);
}

#endif// PARMOON__UTILITIES
