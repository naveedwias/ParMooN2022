#ifndef COMMONROUTINETNSE3D_H
#define COMMONROUTINETNSE3D_H

// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double hK, double* gradU, double* u,
			    double* uConv, double* x, double* y, double* z,
                            double proj_space);
// ========================================================================
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff,
                                          double* params);

// ======================================================================
// compute stabilization for div--div term
// ======================================================================
double DivDivStab3D(double u1, double u2, double u3, double hK, double eps);

// ======================================================================
/// compute the squared Frobenius norm of the turbulent viscosity tensor
// ======================================================================
double frobeniusNormTensor(const double* u, const double* gradu,
                           const double *uConv, int proj_space = 0);

double turbulentViscosityDelta(double hK);

// ======================================================================
/// compute turbulence viscosity using Smagorinsky model
// ======================================================================
double turbulentViscosity3D(double hK, const double* u, const double* gradu,
                            const double* uConv, const double* x,
                            const double* y, const double* z,
                            double proj_space);

/// @brief Estimate SGS kinetic energy k using
/// \nu_T ~= 2 C \delta \sqrt k
double sgsKineticEnergy3D(double hK, const double* u, const double* gradu,
                          const double* uConv, const double* x,
                          const double* y, const double* z,
                          double proj_space, double scale);

#endif // COMMONROUTINETNSE3D_H
