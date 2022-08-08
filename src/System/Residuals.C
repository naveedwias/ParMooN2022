#include <Residuals.h>
#include <cmath>
#include <iomanip>


Residuals::Residuals()
  : momentumResidual(1e10), massResidual(1e10), fullResidual(1e10),
    momentumResidualMax(1e10), massResidualMax(1e10), fullResidualMax(1e10)
{}

Residuals::Residuals(double moR2, double maR2, double moRM, double maRM)
 : momentumResidual(std::sqrt(moR2)), massResidual(std::sqrt(maR2)),
   fullResidual(std::sqrt(moR2 + maR2)),
   momentumResidualMax(moRM), massResidualMax(maRM),
   fullResidualMax(std::max(moRM, maRM))
{
}

std::ostream& operator<<(std::ostream& s, const Residuals& n)
{
  using namespace std;
  s << "\n\t2-residual: " << setprecision(10) << setw(15) << n.fullResidual
    << " (momentum: " << setprecision(6) << setw(10) << n.momentumResidual
    << "  mass: " << setprecision(6) << setw(10) << n.massResidual << ")";

  s << "\n\tmax residual: " << setprecision(10) << setw(15) << n.fullResidualMax
    << " (momentum: " << setprecision(6) << setw(10) << n.momentumResidualMax
    << "  mass: " << setprecision(6) << setw(10) << n.massResidualMax << ")";

  return s;
}
