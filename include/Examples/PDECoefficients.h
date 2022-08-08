#ifndef __PDECOFFICIENTS__
#define __PDECOFFICIENTS__

#include <functional> // std::function
#include <ostream>
#include <Point.h>

using namespace parmoon;

enum class Problem_type
{
  ConvDiffReac,
  Stokes,
  NavierStokes,
  Darcy,
  dummy
};
std::ostream& operator<<(std::ostream& out, const Problem_type value);

/** ************************************************************************
*
* @class      PDECoefficients
* @brief      handle coefficient functions of a pde
*
* Essentially this stores a function which can be evaluated at coordinates in
* the computational domain.
*
* @author     Ulrich Wilbrandt
* @date       03.11.2014
*
****************************************************************************/
class PDECoefficients
{
 public:
  /** @brief construct a PDECoefficients object for time dependent problems
   * 
   * The function f is set to this->coeff. It is recommended to set the second
   * parameter as well, but it is not strictly needed, see the member variable
   * type in this class.
   * 
   * The double parameter is the point in time at which the coefficients can 
   * be evaluated. For stationary problems this is simply left unused.
   */
  PDECoefficients(
      std::function<void(const Point&, double, std::vector<double>&)> f,
      Problem_type t = Problem_type::dummy, bool lhs_depends_on_time = true,
      bool rhs_depends_on_time = true);

  /** @brief evaluate all coefficients of the pde at a given point
   *
   * The minimum size of the vector depends on the problem type, see coeff.
   * 
   * This is a convenience method for stationary problems. It simply calls the
   * other `get` method with a fixed time (which is then probably not used).
   */
  void get(const Point& p, std::vector<double>& v) const
  { this->get(p, 0.0, v); }

  /** @brief evaluate all coefficients of the pde at a given point in space and
   * time.
   *
   * The minimum size of the vector depends on the problem type, see coeff.
   */
  void get(const parmoon::Point&, double, std::vector<double>&) const;

 private:
  /** @brief type of the problem for this example.
   *
   * The number of coefficients and their order depends on the type of problem.
   * If it is not used, it will be set to Problem_type::dummy. Then in the
   * functions this->get(...) no check will be done to ensure the vector is
   * long enough.
   */
  Problem_type type;

  /** @brief the function to evaluate the coefficients of the pde2D
   *
   * The coefficients are written into the vector. The size of the vector is
   * checked outside of this function but within this class. It depends on the
   * problem type at hand and the space dimension. The order of the entries in
   * this vector is as follows, depending on the
   * Problem_type in 2D:
   *
   *  <TABLE>
   *    <TR>
   *      <TD>Problem_type</TD>
   *      <TD>size</TD>
   *      <TD>entries (ordered)</TD>
   *    </TR>
   *    <TR>
   *      <TD>ConvDiffReac</TD>
   *      <TD>5</TD>
   *      <TD>diffusion, convection in x-direction, <BR>
   *          convection in y-direction, reaction, <BR>
   *          right hand side</TD>
   *    </TR>
   *    <TR>
   *      <TD>Stokes</TD>
   *      <TD>4</TD>
   *      <TD>viscosity, right hand side (first component),<BR>
   *          right hand side (second component),<BR>
   *          right hand side (third component, divergence)</TD>
   *    </TR>
   *    <TR>
   *      <TD>NavierStokes</TD>
   *      <TD>4</TD>
   *      <TD>viscosity, right hand side (first component),<BR>
   *          right hand side (second component),<BR>
   *          right hand side (third component, divergence)</TD>
   *    </TR>
   *    <TR>
   *      <TD>Darcy</TD>
   *      <TD>7</TD>
   *      <TD>permeability (top left entry), permeability (top right entry),<BR>
   *          permeability (lower left entry),
   *          permeability (lower right entry),<BR>
   *          right hand side (first component),
   *          right hand side (second component),<BR>
   *          right hand side (third component, divergence)</TD>
   *    </TR>
   *  </TABLE>
   *
   *  And in 3D:
   *  <TABLE>
   *    <TR>
   *      <TD>Problem_type</TD>
   *      <TD>size</TD>
   *      <TD>entries (ordered)</TD>
   *    </TR>
   *    <TR>
   *      <TD>ConvDiffReac</TD>
   *      <TD>6</TD>
   *      <TD>diffusion, convection in x-direction, <BR>
   *          convection in y-direction, convection in z-direction <BR>
   *          reaction, right hand side</TD>
   *    </TR>
   *    <TR>
   *      <TD>Stokes</TD>
   *      <TD>5</TD>
   *      <TD>viscosity, right hand side (first component),<BR>
   *          right hand side (second component),
   *          right hand side (third component) <BR>
   *          right hand side (fourth component, divergence)</TD>
   *    </TR>
   *    <TR>
   *      <TD>NavierStokes</TD>
   *      <TD>5</TD>
   *      <TD>viscosity, right hand side (first component),<BR>
   *          right hand side (second component),
   *          right hand side (third component) <BR>
   *          right hand side (fourth component, divergence)</TD>
   *    </TR>
   *    <TR>
   *      <TD>Darcy</TD>
   *      <TD>13</TD>
   *      <TD>permeability (top left entry),
   *          permeability (top center entry),<BR>
   *          permeability (top right entry),
   *          permeability (middle left entry), <BR>
   *          permeability (middle center entry),
   *          permeability (middle right entry), <BR>
   *          permeability (lower left entry),
   *          permeability (lower center entry), <BR>
   *          permeability (lower right entry),
   *          right hand side (first component), <BR>
   *          right hand side (second component),
   *          right hand side (third component),<BR>
   *          right hand side (fourth component, divergence)</TD>
   *    </TR>
   *  </TABLE>
   *
   */
  std::function<void(const Point&, double, std::vector<double>&)> coeff;
};

#endif // __PDECOFFICIENTS__
