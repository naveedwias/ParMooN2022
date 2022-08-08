#include "NSE_2D/backward_facing_step.h"

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double y, double *values)
{
  values[0] = y >= 0. ? 4*y*(1-y) : 0.;
}

void InitialU2(double, double, double *values)
{
  values[0] = 0;
}

void InitialP(double, double, double *values)
{
  values[0] = 0;
}
