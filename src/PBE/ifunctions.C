#include "ifunctions.h"
#include "constants.h"
#include "attrition.h"
#include "basictools.h"

double nu_frag(double x, double* r, function_3D3D velocity, double t)
{
	return V_attr(x, r, velocity, t) * (std::pow(L_min, -2.25) - std::pow(L_max(x, r, velocity, t), -2.25)) /
		   (3 * k_V * (std::pow(L_max(x, r, velocity, t), 0.75) - std::pow(L_min, 0.75)));
}

double V_attr(double x, double* r, function_3D3D velocity, double t)
{
	return 2 * std::pow(H_V, 2/3) * E_kin(x, r, velocity, t) / (3 * mu * Gamma_Kr);
}

double E_kin(double x, double* r, function_3D3D velocity, double t)
{
	return 0.5 * ro_d * k_V * x * x * x * std::pow(spectral_norm(velocity(r, t), 2), 2);
}

double L_max(double x, double* r, function_3D3D velocity, double t)
{
	return 0.5 * std::pow(H_V, 2/9) * E_kin(x, r, velocity, t) / (std::pow(mu,1/3) * std::pow(Gamma_Kr, 1/3));
}

double L_star(double x, double* r, function_3D3D velocity, double t)
{
	return std::pow(x * x * x - V_attr(x, r, velocity, t) / k_V, 1/3);
}

double b(double* r, function_3D3D velocity, double t)
{
	if(is_on_boundary(r))
		return 0;
	return 1; // TODO
}