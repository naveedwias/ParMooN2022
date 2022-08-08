#ifndef __WINDKESSELBOUNDARY__
#define __WINDKESSELBOUNDARY__

#include <set>
#include <string>
#include <vector>

/** @brief store the variables for lumped parameter bc models
 * 
 * Lumped parameters models are used for boundary conditions when the
 * flow downstream is unknown, but it is expected to have an influence on 
 * the upstream dynamics. These models impose a special Neumann boundary condition in
 * which the Neumann value P(t) depends on the outgoing flow rate via an ODE:
 * P(t) = W(Q,d_t Q, ....)
 *
 */

class WindkesselBoundary
{
  public:

    /// @brief Create object for given parameters
    WindkesselBoundary(double _Rp, double _C, double _Rd, double _pi = 0.0)
    {
        Rp = _Rp;
        C = _C;
        Rd = _Rd;

        pi = _pi;
        last_pi = _pi;
        pressure = 0.0;
    };

    /// @brief destructor
    ~WindkesselBoundary() = default;

    /// @brief set parameters
    void set_proximal_resistance(double _Rp) { Rp = _Rp; }
    void set_distal_resistance(double _Rd) { Rd = _Rd; }
    void set_capacitance(double _C) { C = _C; }

    double get_proximal_resistance() const { return Rp; }
    double get_distal_resistance() const { return Rd; }
    double get_capacitance() const { return C; }

    /// @brief compute Neumann value
    double solve(double q, double dt);

    /// @brief advance to next time step
    void advance();

    /// @brief return last computed Neumann value
    double get_pressure() const;

  protected:

    double Rp;
    double Rd;
    double C;

    double last_pi;

    double tau;
    double pi;
    double pressure;
};

#endif // __WINDKESSELBOUNDARY__
