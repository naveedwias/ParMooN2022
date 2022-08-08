/** ****************************************************************************
*
* @date   30.01.2021
* @brief  
* @author Baptiste Moreau
*
***************************************************************************** */

#ifndef __PERIODIC_DOMAIN__
#define __PERIODIC_DOMAIN__

#include <Constants.h>
#include <Domain.h>
#include <Point.h>
#include <templateNames.h>

#include <vector>


/** @brief set the norm of the vector to 1 */
void normalized_vector(std::vector<double>& v, double tolerance=1e-10);

/** @brief the norm of the vectors should be larger than tol */
bool are_vectors_collinear(const std::vector<double>& a,
                           const std::vector<double>& b,
                           double tol=1e-10);

bool are_vectors_equal(const std::vector<double>& a,
                       const std::vector<double>& b,
                       double tol=1e-10);

/** @brief check if two sorted vectors<int> contain at least
 *         one common element */
bool have_common_element(std::vector<int>& V1, std::vector<int>& V2);

/** @brief sort elements of the vector and keep only one of each */
void sort_and_select(std::vector<int>& V);

/** @brief search periodic boundary joints
 *  @param[in]     BoundaryCondition
 *  @param[in]     Domain
 *  @param[in,out] joint_list cell-joint identifier for periodic boundary
 *  @param[in,out] barycenter barycenter of the periodic joint
 *  @return        true if periodic boundaries are found
 */
#ifdef __2D__
bool find_PeriodicJoint(BoundCondFunct2D* BoundaryCondition,
#else
bool find_PeriodicJoint(BoundCondFunct3D* BoundaryCondition,
#endif
                        TDomain *Domain,
                        std::vector<std::pair<int, int>>& joint_list,
                        std::vector<parmoon::Point>& barycenter);

/** @brief change the boundary joints to periodic joints */
#ifdef __2D__
void set_PeriodicJoint(BoundCondFunct2D* BoundaryCondition,
#else
void set_PeriodicJoint(BoundCondFunct3D* BoundaryCondition,
#endif
                       TDomain *Domain);

#endif
