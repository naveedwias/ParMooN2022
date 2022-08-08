#ifndef __NON_NEWTONIAN__
#define __NON_NEWTONIAN__

#include <Database.h>

typedef enum
{
  Newtonian,
  CarreauYasuda,
  PowerLaw,
  Cross,
  HerschelBulkley
} ViscosityMode;

struct ViscositySettings
{
  public:
    ViscositySettings();
    ViscositySettings(const ParameterDatabase &db);

    ViscosityMode mode;

    double nu_0;
    double nu_infty;

    double tau_0;

    double k;
    double n;

    double lambda;
    double a;

    double g_0;
};

template<int d>
double NonNewtonianViscosity(const double* grad_u,
  const ViscositySettings& settings);

double CarreauYasudaViscosity(double shear_rate, const ViscositySettings& settings);
double PowerLawViscosity(double shear_rate, const ViscositySettings& settings);
double CrossViscosity(double shear_rate, const ViscositySettings& settings);
double HerschelBulkleyViscosity(double shear_rate, const ViscositySettings& settings);

template<int d>
void NSLaplaceGradGradSingleNonNewtonian(double Mult, const double *,
                                        const double *param, double,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **,
                                        const ViscositySettings& settings);

template<int d>
void NSLaplaceGradGradNonNewtonian(double Mult, const double *,
                                        const double *param, double,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **,
                                        const ViscositySettings& settings);

template<int d>
void NSLaplaceDeformationNonNewtonian(double Mult, const double *,
                                        const double *param, double,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **,
                                        const ViscositySettings& settings);

#endif