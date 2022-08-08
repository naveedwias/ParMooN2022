#include "Database.h"
#include <WindkesselBoundary.h>
#include <algorithm>
#include <iomanip>      // std::setw
#include <sstream>
#include <iterator>
#include <vector>
#include <MooNMD_Io.h>

double WindkesselBoundary::solve(double q, double delta_t)
{
  /* Backwards Euler solver for the RCR model's distal pressure, which evolves as:
   *
   * d/dt pi = 1/C Q - 1/(R_d C) pi
   *
   * so for one backwards Euler step we have
   *
   * (pi_new - pi) / tau = 1/C Q - 1/(R_d C) pi_new
   *
   * (1 + tau / R_d C) pi_new = tau/C Q + pi
   *
   * (1 + beta) pi_new = R_d Q + beta pi
   *
   * with
   *
   * beta = R_d C / tau.
   *
   * The overall pressure is p = pi + R_p Q.
   *
   * Notice that if C = 0, this correctly resolves to the resistive model
   * p = (R_d + R_p) Q.
   */

  double beta = this->C * this->Rd / delta_t;
  double pi_new = 1.0 / (1.0 + beta) * (beta * last_pi + this->Rd * q);
  double p = pi_new + q * this->Rp;

  tau = delta_t;
  pi = pi_new;
  pressure = p;

  return get_pressure();
}

double WindkesselBoundary::get_pressure() const
{
  return pressure;
}

void WindkesselBoundary::advance()
{
  last_pi = pi;
}