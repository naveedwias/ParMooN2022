// =======================================================================
// @(#)ConvDiff.C
//
// Purpose: contains routines which are relevant for ConvDiff2D.C and
//          ConvDiff3D.C
//
// Author: Volker John
//
// History: start of implementation 09.08.2011
//
// =======================================================================

#include <ConvDiff.h>

#include <Database.h>
#include <MainUtilities.h>
#include <LinAlg.h>
#ifdef __2D__
  #include <FEFunction2D.h>
#endif
#include "BaseCell.h"

#include <cmath>

template<>
double Mesh_size_in_convection_direction_without_storing<2>(
  double hK, std::array<double, 2> convection)
{
  int i;
  double x[4], y[4], sx, sy, a[16], b[16], den, val, norm_b;

  // compute numerator
  norm_b = std::sqrt(convection[0]*convection[0] + convection[1]*convection[1]);
  
  // triangles
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3] == -4711)
  {
    for (i=0;i<3;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
    }
    // initialize rhs
    memset(b,0,9*sizeof(double));

    // set matrices for computation of the coefficients
    // of the linear function
    for (i=0;i<3;i++)
    {
      a[3*i] = 1;
      a[3*i+1] = x[i];
      a[3*i+2] = y[i];
      b[4*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,3,3,3,3);

    // compute denominator
    den = 0;
    for (i=0;i<3;i++)
    {
      // value of gradient basis fct. in bary centre is a constant
      den += std::abs(convection[0]*b[3*i+1] + convection[1]*b[3*i+2]);
    }
    // return the mesh size in convection direction
    if (den<1e-10)
      return(hK);
    else
      return(2*norm_b/den);
  }
  else
  {                                               // quadrilateral
    sx = sy = 0;
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      //Output::print(x[i], " ", y[i]);
      sx += x[i];
      sy += y[i];
    }
    //bary centre
    sx /= 4;
    sy /= 4;
    // initialize rhs
    memset(b,0,16*sizeof(double));

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = x[i]*y[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,4,4,4,4);

    // compute denominator
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      val = convection[0]*(b[4*i+1] + b[4*i+3] * sy);
      val += convection[1]*(b[4*i+2] + b[4*i+3] * sx);
      den += std::abs(val);
    }
    // return the mesh size in convection direction
    if (den<1e-10)
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
}

template<>
double Mesh_size_in_convection_direction_without_storing<3>(
  double hK, std::array<double, 3> convection)
{
  int i;
  double x[8], y[8], z[8], sx, sy, sz, a[64], b[64], den, val, norm_b;

  // compute numerator
  norm_b = std::sqrt(convection[0]*convection[0] + convection[1]*convection[1]
                     + convection[2]*convection[2]);
  
  // tetrahedra
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
  {
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
    }
    // initialize rhs
    memset(b,0,16*sizeof(double));

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = z[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    // which is faster ?
    // SolveMultipleSystemsLapack(a,b,4,4,4,4);
    SolveMultipleSystems(a,b,4,4,4,4);
    
    // compute denominator 
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      // is a constant
      den += std::abs(convection[0]*b[4*i+1] + convection[1]*b[4*i+2]
                      +convection[2]*b[4*i+3]);
    } 
    // return the mesh size in convection direction
    if (den<1e-10)
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
  else
  { // hexahedron 
    sx = sy = sz = 0;
    for (i=0;i<8;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
      //Output::print(x[i], " ", y[i]);
      sx += x[i];
      sy += y[i];
      sz += z[i];
    }
    //bary centre
    sx /= 8;
    sy /= 8;
    sz /= 8;
    // initialize rhs
    memset(b,0,64*sizeof(double));
    for (i=0;i<8;i++)
    {
      b[9*i] = 1;
    }
    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<8;i++)
    {
      a[8*i] = 1;
      a[8*i+1] = x[i];
      a[8*i+2] = y[i];
      a[8*i+3] = z[i];
      a[8*i+4] = x[i]*y[i];
      a[8*i+5] = x[i]*z[i];
      a[8*i+6] = y[i]*z[i];
      a[8*i+7] = x[i]*y[i]*z[i];
    }

    // solve system for the coefficients of the bilinear function
    // which is faster ??
    //SolveMultipleSystemsLapack(a,b,8,8,8,8);
    SolveMultipleSystems(a,b,8,8,8,8);

    // compute denominator 
    den = 0;
    for (i=0;i<8;i++)
    {
      // value of gradient basis fct. in bary centre
      val = convection[0]*(b[8*i+1]+ b[8*i+4]*sy + b[8*i+5]*sz + b[8*i+7]*sy*sz);
      val += convection[1]*(b[8*i+2]+ b[8*i+4]*sx + b[8*i+6]*sz + b[8*i+7]*sx*sz);
      val += convection[2]*(b[8*i+3]+ b[8*i+5]*sx + b[8*i+6]*sy + b[8*i+7]*sx*sy);
      den += std::abs(val);
    } 
    // return the mesh size in convection direction
    if (den<1e-10)
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
}

template<int d>
double Mesh_size_in_convection_direction(double hK, std::array<double, d> b)
{
  if(TDatabase::ParamDB->INTERNAL_HK_CONVECTION < 0)
  {
    // not yet computed for this mesh cell
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = 
      Mesh_size_in_convection_direction_without_storing<d>(hK, b);
  }
  // else: already computed for this mesh cell
  return TDatabase::ParamDB->INTERNAL_HK_CONVECTION;
}
template 
double Mesh_size_in_convection_direction<3>(double hK, std::array<double, 3> b);
template 
double Mesh_size_in_convection_direction<2>(double hK, std::array<double, 2> b);


/******************************************************************************/
//
// definitions of parameters for the SUPG method
// steady-state and time-dependent convection-diffusion-reaction equation
//
/******************************************************************************/

template <int d>
double Compute_SDFEM_delta(double hK, double eps, 
                           std::array<double, d> convection,
                           double react, double linfb)
{
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double alpha, alpha2, delta, nu, reaction, h_K;
  double norm_b = 0;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  int i;
  
  if(std::abs(eps)<1e-20)
    eps = 1e-20;

  // compute cell diameter in convection direction

  if(TDatabase::ParamDB->CELL_MEASURE==4)
    h_K = Mesh_size_in_convection_direction<d>(hK, convection);
  else
    h_K = hK;

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:
    case 1:
    case 2:
    case 5:
    case 6:
    case 7:
    case 8:
    case 11:
      norm_b = convection[0]*convection[0] + convection[1]*convection[1];
      if(d == 3)
        norm_b += convection[2]*convection[2];
      norm_b = std::sqrt(norm_b);
      break;
  }

  // just for safety
  reaction = std::abs(react);

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:
      // version from Roos, Stynes, Tobiska 2006
      if(eps < h_K*norm_b)
        delta = delta0 * h_K/norm_b;
      else
        delta = delta1 *h_K*h_K/eps ;
      break;
    case 1:
    case 2:
      // delta based on 1d Green's formula
      // case 2 is the the standard choice, in this choice
      // - h_K is mesh size in convection direction
      // - delta0 = 1
      // are forced
      nu = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      if (norm_b > 0)
      {
        alpha = nu*norm_b*h_K/(2*eps);
        delta = nu*h_K*(1/std::tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
        delta = 0;
      // scaling
      delta *= delta0;
      // for computation of reference solutions for parameter optimization
      //if (TDatabase::ParamDB->INTERNAL_MOMENT == 10)
    	//  delta = 0;
      break;

      // THESE ARE THE PARAMETERS WHERE REACTION BECOMES IMPORTANT
      // can be used for steady-state or time-dependent case
    case 4:
      // Lube/Rapin (2006)
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 0 (e.g. linear and bilinear elements)
      // polynomial degree = 1

      delta = -1;

      if (linfb > 0)
        delta = h_K/linfb;

      if (reaction > 0)
      {
        alpha = 1.0/reaction;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      if (eps>0)
      {
        alpha = h_K*h_K/eps;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      delta *= delta0;
      if (delta<0)
        delta = 0;
      break;
    case 5:
    case 8:
      // Lube/Rapin (2006), slightly modified
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 0 (e.g. linear and bilinear elements)
      // polynomial degree = 1
      // !!! case 8: reaction without term 1
      delta = -1;
      if (norm_b > 0)
        delta = h_K/(2*norm_b);

      if (reaction > 0)
      {
        alpha = 1.0/reaction;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      if (eps>0)
      {
        alpha = h_K*h_K/eps;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      delta *= delta0;
      if (delta<0)
        delta = 0;
      break;

    case 6:
      // Franca/Valentin (2000)
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 1/3 (e.g. linear and bilinear elements)
      // polynomial degree = 1
      if (reaction<1e-20)
        reaction = 1e20;

      alpha  = 6*eps/(h_K*h_K*reaction);
      if(alpha <= 1)
        alpha = 1;

      alpha2 = h_K*norm_b/(3*eps);
      if(alpha2 <= 1)
        alpha2 = 1;

      delta=1/(reaction*h_K*h_K*alpha+6*eps*alpha2);
      delta *= delta0*h_K*h_K;
      break;

    case 7:
      // Codina (2000)
      // equations with reaction
      delta = delta0 * h_K * h_K/(4*eps+2 * h_K * norm_b + h_K*h_K * reaction);
      break;

      // THESE ARE PARAMETERS FOR THE TIME-DEPENDENT CASE
    case 9:
      // first estimate in paper John, Novo, SIMUM 2011, formula (3.3)
      delta = delta0*time_step;
      break;
    case 10:
      // second estimate in paper John, Novo, SINUM 2011, asymptotic of formula (5.1)
      // get unscaled diffusion
      eps /= (theta1*time_step);
      if(eps <= h_K)
        delta = delta0 * h_K;
      else
        delta = delta0 * h_K*h_K/eps;
      break;
    case 11:
      // for estimate in paper John, Novo, SINUM 2011, formula (3.10)
      // parameter depending on square root of time step
      norm_b /= time_step;
      if (norm_b >0)
        delta = delta0 * h_K * std::sqrt(time_step)/norm_b;
      else
        delta = 0;
      break;
    case 12:
    {
      // for estimate in paper John, Novo, 2021, formula (35)
      // norm_b := max(abs(b_i))          // norm L_inf
      // norm_b = std::max(std::abs(convection[0]), std::abs(convection[1]));
      // if(d == 3)
      //   norm_b = std::max(norm_b, std::abs(convection[2]));
      // norm_c := abs(react) = reaction  // norm L_inf
      // double min_1 = std::min(0.5, mu_0/(4.norm_c));
      // double min_2 = std::min(sqrt(mu_0/norm_c), norm_b*h_K/(4.*eps*c_inv));
      // delta = std::min(h_K/(4.*c_inv*norm_b)*std::min(min_1, min_2),
      //                  1./mu_0, 1.norm_c);
      // Values for exemple traveling_wave_sqrt.h are:
      //   norm_b = max(cos(Pi/3), sin(Pi/3)) = sqrt(3)/2
      //   norm_c = 1
      //   mu_0   = 1 (because b is constant)
      norm_b = std::sqrt(3)/2;
      double c_inv = delta0;

      delta = 100*std::min(h_K/(2.*sqrt(3)*c_inv)*std::min(0.25,
                                                     sqrt(3)*h_K/(8*eps*c_inv)),
                       1.);
      break;
    }
    case 100:
      // take value from piecewise constant field
      // only on finest level available
      i = TDatabase::ParamDB->INTERNAL_LOCAL_DOF;
      delta =  TDatabase::ParamDB->INTERNAL_P1_Array[i];
      break;
    default :
      ErrThrow("SDFEM_TYPE ", TDatabase::ParamDB->SDFEM_TYPE, 
        " not implemented !!!");
  }

  Output::print<5>("delta ", delta);
  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 120814)
    TDatabase::ParamDB->INTERNAL_P1_Array[TDatabase::ParamDB->INTERNAL_MOMENT] = delta;

  return(delta);
}


// ========================================================================
// compute parameters for SOLD schemes
// 2D stationary cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = ||u^h||_{H^1,K}
//        param[6] = ||R(u^h)||_{L^2,K}
// 2D time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = u_old
//        param[6] = u_old_x
//        param[7] = u_old_y
//        param[8] = u_old_xx
//        param[9] = u_old_yy
// 3D stationary cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = ||u^h||_{H^1,K}
//        param[5] = ||R(u^h)||_{L^2,K}
// 3D time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = u_old
//        param[5] = u_old_x
//        param[6] = u_old_y
//        param[7] = u_old_z
//
// The implementation is for the 3D case. The unnecessary flops for 2D
// should be neglibible in the computing time
// ========================================================================

template <int d>
double Compute_SOLD_sigma(double hK, double eps, 
                          std::array<double, d> convection, double c, double f,
                          double, double deltaK, double *param,
                          double residual, int residual_computed,
                          int time_dependent_problem)
{
  if(param == nullptr)
  {
    ErrThrow("Compute_SOLD_sigma needs param to be valid but it's nullptr");
  }
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE, i, N;
  double u_x, u_y, u_z, sigma, res=0.0, norm_b2, value;
  double b1_orth, b2_orth,  b3_orth, norm_der_u2, linfb_orth, z1, z2, z3, linfz, normz;
  double alpha, beta, gamma, lambda, kappa, omega, rho, norm_b_orth;
  double epsilon= 1e-10, hK_project, y, z, u_xx=0., u_yy=0.;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  std::array<double, 3> b;
  if(d == 3) b = {{convection[0], convection[1], convection[2]}};
  if(d == 2) b = {{convection[0], convection[1], 0.}};
  u_z = d == 2 ? 0.0 : param[3];
  // square of norm of convection
  norm_b2 =  b[0]*b[0] + b[1]*b[1] + b[2]*b[2];

  // the residual is already computed
  // if neither residual nor params are provided: initialization is res=0
  if (residual_computed)
    res = residual;

  u_x = param[1];
  u_y = param[2];
  // square of the norm of the derivative of the current solution
  norm_der_u2 = u_x*u_x + u_y*u_y + u_z*u_z;

  // compute the residual if this is not already done
  // in 3D, the Laplacian is not provided and the corresponding factors
  // were set to be zero
  if (!residual_computed)
    // only for stationary problems
    res = - eps*(u_xx + u_yy) + b[0]*u_x + b[1]*u_y + +b[2] *u_z + c*param[0] - f;

  // compute the parameter for the SOLD scheme
  switch (sold_parameter_type)
  {
    case JSW87:                                   // Johnson,Schatz,Wahlbin (1987) (linear)
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction<2>(hK, {{b[0], b[1]}});
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction<3>(hK, b);
#endif
      sigma = std::sqrt(norm_b2)*hK_project*std::sqrt(hK_project) - eps;
      if (sigma < 0)
        sigma = 0;
      break;

    case HMM86:                                   // Hughes, Mallet, Mizukami (1986)
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      value = b[0]*u_x + b[1]*u_y + b[2]* u_z;       // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
      // ||b_orth||_\infty
      if (std::abs(b2_orth)< std::abs (b1_orth))
        linfb_orth = std::abs(b1_orth);
      else
        linfb_orth = std::abs(b2_orth);
      if (linfb_orth < std::abs(b3_orth))
        linfb_orth = std::abs(b3_orth);
      // \tau(\b_orth)
#ifdef __2D__
      value = Compute_SDFEM_delta<2>(hK, eps, {{b1_orth, b2_orth}}, c, linfb_orth);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta<3>(hK, eps, {{b1_orth, b2_orth, b3_orth}}, c, linfb_orth);
#endif
      // sigma = max ( \tau(b_orth) - \tau(b))
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      //Output::print("sigma ", sigma, " orth ", value, " deltaK ", deltaK);
      //Output::print("u_x ", u_x, " u_y ", u_y);
      value = b[0]*u_x + b[1]*u_y + b[2]*u_z;           // b \cdot nabla u^h
      if (norm_der_u2>0)
        // sigma = sigma * residual * (b\cdot u^h)/||\nabla u^h||^2
        sigma *= res*value/norm_der_u2;
      else
        sigma = 0;
      break;

    case TP86_1:                                  // Tezduyar, Park (1986), first parameter
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b[0]*u_x + b[1]*u_y + b[2]*u_z;           // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;

      rho =  std::sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/std::sqrt(norm_b2);
      value = 2*value*(1-value);
#ifdef __2D__
      kappa = Mesh_size_in_convection_direction_without_storing<2>(hK,{{b1_orth,b2_orth}});
#endif
#ifdef __3D__
      kappa = Mesh_size_in_convection_direction_without_storing<3>(hK,{{b1_orth,b2_orth,b3_orth}});
#endif
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * value/(2*rho)*res*alpha/norm_der_u2;
      break;

    case TP86_2:                                  // Tezduyar, Park (1986), second parameter
      if (TDatabase::ParamDB->SOLD_U0==0)
      {
        ErrThrow("Paramter SOLD_U0 is zero ");
      }
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b[0]*u_x + b[1]*u_y + b[2]*u_z;           // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;

      rho =  std::sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/std::sqrt(norm_b2);
      value = 2*value*(1-value);
#ifdef __2D__
      kappa = Mesh_size_in_convection_direction_without_storing<2>(hK,{{b1_orth,b2_orth}});
#endif
#ifdef __3D__
      kappa = Mesh_size_in_convection_direction_without_storing<3>(hK,{{b1_orth,b2_orth,b3_orth}});
#endif
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * kappa* value/(2*rho)*res*alpha/std::sqrt(norm_der_u2);
      sigma /= TDatabase::ParamDB->SOLD_U0;
      break;

    case GdC88:                                   // Galeao, do Carmo (1988)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
        z3 *= theta1*time_step;
      }
      if (std::abs(z2)< std::abs (z1))
        linfz = std::abs(z1);
      else
        linfz = std::abs(z2);
      if (linfz < std::abs(z3))
        linfz = std::abs(z3);
#ifdef __2D__
      value = Compute_SDFEM_delta<2>(hK, eps, {{z1, z2}}, c, linfz);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta<3>(hK, eps, {{z1, z2, z3}}, c, linfz);
#endif
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;

    case dCG91:                                   // do Carmo, Galeao (1991)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= std::sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)
        normz *= theta1*time_step;
      sigma = deltaK*(std::sqrt(norm_b2)/normz-1);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      //Output::print(deltaK, " ", std::sqrt(norm_b2), " ", normz, " ", sigma);
      break;

    case dCA03:                                   // do Carmo, Alvarez (2003)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
        z3 *= theta1*time_step;
      }
      if (std::abs(z2)< std::abs (z1))
        linfz = std::abs(z1);
      else
        linfz = std::abs(z2);
      if (linfz < std::abs(z3))
        linfz = std::abs(z3);
#ifdef __2D__
      value = Compute_SDFEM_delta<2>(hK, eps, {{z1, z2}}, c, linfz);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta<3>(hK, eps, {{z1, z2, z3}}, c, linfz);
#endif
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }

      // sigma of [GdC88] is computed, now compute rho
      normz = std::sqrt(z1*z1+z2*z2+z3*z3);

      alpha = normz/std::sqrt(norm_b2);
      if (alpha >= 1)
      {
        rho = 1;
      }

#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction<2>(hK, {{b[0], b[1]}});
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction<3>(hK, b);
#endif
      beta = std::pow(hK_project,1-alpha*alpha);
      if (beta > 1)
        beta = 1;

      gamma = (alpha+beta)/2.0;
      if (gamma > beta)
        gamma = beta;

      if (alpha > std::abs (res))
        value = alpha;
      else
        value = std::abs(res);
      lambda = std::pow(value,3+alpha/2+alpha*alpha);
      if (0.25+alpha>0.5)
        value = 0.25+alpha;
      else
        value = 0.5;
      lambda /= std::pow(gamma, value);
      if (lambda >= 1)
      {
        rho = 1;
      }

      value = (1-lambda)/(1+lambda);
      kappa = std::pow(std::abs(2-lambda),value);
      kappa -= 1;

      value = std::pow(gamma,2-alpha*alpha);
      omega = alpha*alpha * value/deltaK;

      if ((alpha < 1) && (lambda < 1))
      {
        rho = std::pow(omega*sigma,kappa);
      }
      sigma *= rho;
      sigma *= res*res/norm_der_u2;
      break;

    case AS97:                                    // Almeida, Silva (1997)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= std::sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)
      {
        normz *= theta1*time_step;
      }
      value = b[0]*u_x + b[1]*u_y + b[2]*u_z;

      if (res==0)
        value = 1;
      else
        value /= res;

      if (value < 1)
        value = 1;
      sigma = deltaK*(std::sqrt(norm_b2)/normz-value);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;

    case C93:                                     // Codina (1993)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      value = b[0]*u_x + b[1]*u_y + b[2]*u_z;
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
      norm_b_orth = std::sqrt(b1_orth*b1_orth+b2_orth*b2_orth+b3_orth*b3_orth);
      if (norm_b_orth == 0)
      {
        sigma = 0;
        break;
      }
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(hK * norm_b_orth);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * std::abs(res)/(2*std::sqrt(norm_der_u2));
      break;

    case KLR02_1:                                 // Knopp, Lube, Rapin (2002)
      //case KLR02_3:                      // same as KLR02_1 with SOLD_S = 0, new version see below
#ifdef __2D__
      i = 5;
#endif
#ifdef __3D__
      i = 4;
#endif
      value = TDatabase::ParamDB->SOLD_S + param[i];
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      alpha = param[i+1] / value;
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;

    case KLR02_3:                                 //pointwise evaluation of residual and norm of u_h
      if (std::sqrt(norm_der_u2)<epsilon)
      {
        sigma = 0;
        break;
      }
      sigma =  TDatabase::ParamDB->SOLD_CONST * hK *std::abs(res)/(2*std::sqrt(norm_der_u2)) - eps;
      if (sigma<0)
        sigma = 0;
      break;

    case KLR02_4:                                 //pointwise evaluation of residual and norm of u_h with additive constant
      // THIS IS NOT KLR02_4 FROM John/Knobloch 2007 !!!
      if (std::sqrt(norm_der_u2)<epsilon)
      {
        sigma = 0;
        break;
      }
      sigma =  TDatabase::ParamDB->SOLD_CONST * hK *std::abs(res)/(2*(TDatabase::ParamDB->SOLD_S +
        std::sqrt(norm_der_u2))) - eps;
      if (sigma<0)
        sigma = 0;
      break;

    case KLR02_2:                                 // similar to KLR02_3
    case CS99:
      // norm of convection
      value = std::sqrt(norm_der_u2);
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      // Q = res/|b|
      alpha = std::abs(res) / value;
      // C - (2 eps)/(h Q)
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;

    case J90:                                     // Johnson (1990)
      alpha = TDatabase::ParamDB->SOLD_CONST;
      kappa = TDatabase::ParamDB->SOLD_POWER;
      sigma = alpha * std::pow(hK,kappa)*std::abs(res)-eps;
      if (sigma < 0)
        sigma = 0;
      break;

    case BE02_1:                                  // Burman, Ern 2002
      alpha = M_PI/6;
      res = res*std::tanh(res/2.0);
      value = std::sqrt(norm_b2)*std::sqrt(norm_der_u2)+std::abs(res);
      if (value > 0)
        sigma = deltaK * norm_b2*std::abs(res)/value;
      else
      {
        sigma = 0;
        break;
      }
      if (norm_b2>0)
      {
        z1 = (1-b[0]*b[0]/norm_b2)*u_x - b[0]*(b[1]*u_y+ b[2]*u_z)/norm_b2;
        z2 = -b[1]*(b[0]*u_x + b[2]*u_z)/norm_b2 + (1-b[1]*b[1]/norm_b2)*u_y;
        z3 = -b[2]*(b[0]*u_x + b[1]*u_y)/norm_b2 + (1-b[2]*b[2]/norm_b2)*u_z;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = std::abs(res)+std::tan(alpha)*std::sqrt(norm_b2)*std::sqrt(z1*z1+z2*z2+z3*z3);
      if (value > 0)
        sigma *= (std::sqrt(norm_b2)*std::sqrt(norm_der_u2) + value)/value;
      else
      {
        sigma = 0;
      }
      break;

    case BE02_2:                                  // modified Burman, Ern 2002
      value = std::sqrt(norm_b2)*std::sqrt(norm_der_u2)+std::abs(res);
      if (value > 0)
        sigma = deltaK * norm_b2*std::abs(res)/value;
      else
        sigma = 0;
      break;

    case BE02_3:                                  // Burman, Ern 2002, formula (29)
      alpha = M_PI/6;
      res = res*std::tanh(res/2.0);
      if (norm_b2>0)
      {
        z1 = (1-b[0]*b[0]/norm_b2)*u_x - b[0]*(b[1]*u_y+ b[2]*u_z)/norm_b2;
        z2 = -b[1]*(b[0]*u_x + b[2]*u_z)/norm_b2 + (1-b[1]*b[1]/norm_b2)*u_y;
        z3 = -b[2]*(b[0]*u_x + b[1]*u_y)/norm_b2 + (1-b[2]*b[2]/norm_b2)*u_z;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = std::sqrt(res*res + std::tan(alpha)*std::tan(alpha)*norm_b2 * (z1*z1+z2*z2+z3*z3));
      if (value > 0)
        sigma = deltaK * norm_b2*std::abs(res)/value;
      else
        sigma = 0;
      break;

    case Y_Z_BETA:                                // Tezduyar
      y = std::abs(TDatabase::ParamDB->SOLD_U0);
      beta = TDatabase::ParamDB->SOLD_POWER;
      if (y==0)
      {
        sigma = 0;
        break;
      }
      if (norm_der_u2==0)
      {
        sigma = 0;
        break;
      }
      z = std::abs(res)/y;
      norm_der_u2 /= (y*y);
      norm_der_u2 = std::pow(norm_der_u2,beta/2.0-1);
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction_without_storing<2>(hK, {{u_x, u_y}})/2.0;
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction_without_storing<3>(hK, {{u_x, u_y, u_z}})/2.0;
#endif
      sigma = TDatabase::ParamDB->SOLD_CONST * z * norm_der_u2 * std::pow(hK_project,beta);
      break;

    case JSW87_1:                                 // Johnson,Schatz,Wahlbin (1987) (linear)
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction<2>(hK, {{b[0], b[1]}});
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction<3>(hK, b);
#endif
      sigma = std::sqrt(norm_b2)*hK_project*std::sqrt(hK_project) - eps;
      sigma /= theta1*time_step;
      if (sigma < 0)
        sigma = 0;
      break;

    case BH04:                                    // Burman, Hansbo 2004, edge stabilization
    case BE05_1:                                  // Burman, Hansbo 2004, edge stabilization
    case BE05_2:                                  // Burman, Hansbo 2004, edge stabilization
    case LP96:                                    // Layton, Polman 1996
    case MH_Kno06:                                // improved Mizukami-Hughes, by Knobloch 2006
    case FEM_TVD:                                 // algebraic flux correction
      sigma = 0;
      break;

    case GENERAL_SOLD:
      // TDatabase::ParamDB->SOLD_CONST is the eta from the SOLD2-paper
      if (norm_der_u2>0)
        sigma = TDatabase::ParamDB->SOLD_CONST*hK*std::abs(res)/(2*std::sqrt(norm_der_u2));
      else
        sigma = 0.0;
      break;

    case 100:
      // take value from piecewise constant field
      // only on finest level available
      i = TDatabase::ParamDB->INTERNAL_LOCAL_DOF;
      N = TDatabase::ParamDB->INTERNAL_ARRAY_LENGTH;
      sigma =  TDatabase::ParamDB->INTERNAL_P1_Array[N+i];
      return(sigma);
      break;

      // for lower levels
    case 101:
      return(0.0);
      break;

    default :
      ErrThrow("SOLD type ", sold_parameter_type, " not available");
  }
#ifdef __2D__
  i = 5;
#endif
#ifdef __3D__
  i = 4;
#endif
  // scaling of sigma accordingly to Knopp, Lube, Rapin (2002)
  if (TDatabase::ParamDB->SOLD_PARAMETER_SCALING)
  {
    value = TDatabase::ParamDB->SOLD_S + param[i];
    if (value == 0)
      sigma = 0;
    else
    {
      alpha = param[i+1] / value;
      sigma *= alpha*alpha;
    }
  }
  else
    sigma *= TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR;

  return (sigma);
}

template
double Compute_SOLD_sigma<2>(double hK, double eps, 
                             std::array<double, 2> convection, double c,
                             double f, double linfb, double deltaK,
                             double *param, double residual,
                             int residual_computed, int time_dependent_problem);
template
double Compute_SOLD_sigma<3>(double hK, double eps, 
                             std::array<double, 3> convection, double c,
                             double f, double linfb, double deltaK,
                             double *param, double residual,
                             int residual_computed, int time_dependent_problem);

/*************************************************************************/
// estimate of coercivity constant 
// used eg for residual based estimator of Verf"uhrt 2005
/*************************************************************************/
double EstimateCoercivityConstant(TCollection *Coll,
#ifdef __2D__
                                  const CoeffFct2D& Coeffs
#else // 3D
                                  const CoeffFct3D& Coeffs
#endif
                                  )
{
  int i, j, N_Cells, N_V;
  double coerc = 4711.0, x, y, xs, ys, *coeff;
  #ifdef __3D__
  double z, zs;
  #endif
  TBaseCell *cell;
  TVertex *vertex;

  coeff = new double[13];
  // index of the reaction term in coeff array
  #ifdef __2D__
  int reac_index = 3;
  #else
  int reac_index = 4;
  #endif
  
    
  N_Cells = Coll->GetN_Cells();                   // number of mesh cells
  // loop over the cells
  for (i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    // loop over the vertices
    xs = ys = 0;
#ifdef __3D__
    zs = 0;
#endif
    for (j=0;j<N_V;j++)
    {
      vertex = cell->GetVertex(j);
      #ifdef __2D__
      vertex->GetCoords(x, y);
      Coeffs(1, &x, &y, nullptr, &coeff);
      #else
      vertex->GetCoords(x, y, z);
      Coeffs(1, &x, &y, &z, nullptr, &coeff);
      #endif
      // assuming that convection is divergence-free
      if (coeff[reac_index] < coerc)
        coerc = coeff[reac_index];
      xs += x;
      ys += y;
      #ifdef __3D__
      zs += z;
      #endif
    }
    xs /= N_V;
    ys /= N_V;
    #ifdef __3D__
    zs /= N_V;
    #endif
    
    #ifdef __2D__
    Coeffs(1, &xs, &ys, nullptr, &coeff);
    #else
    Coeffs(1, &xs, &ys, &zs, nullptr, &coeff);
    #endif
    // assuming that convection is divergence-free
    if (coeff[reac_index] < coerc)
      coerc = coeff[reac_index];
  }

  delete coeff;
  Output::print("coercivity constant (assuming div-free convection): ", coerc);
  return(coerc);
}


void SetSoldParameters(int i)
{
  const int max = 16;
  TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;
  TDatabase::ParamDB->SOLD_TYPE = 1;
  
  if ((i>max-1) && (i<2*max))
  {
    i -= max;
    TDatabase::ParamDB->SOLD_TYPE = 2;
  }

  switch(i)
  {
    case 0:
      TDatabase::ParamDB->SOLD_TYPE = 0;
      break;
    case 1:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
      break;
    case 2:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 1;
      break;
    case 3:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 3;
      break;
    case 4:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 4;
      break;
    case 5:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 5;
      break;
    case 6:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 6;
      break;
    case 7:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      break;
    case 8:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      break;
    case 9:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      break;
    case 10:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 11:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 12:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 13:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 1.0;
      break;
    case 14:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 1.0;
      break;
    case 15:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 2.0;
      break;
    default:
      ErrThrow("Wrong input in SetSoldParameters");
  }
}


#ifdef __2D__
void EdgeStabilization(TFESpace2D *fespace, TFEFunction2D *u, 
                       const CoeffFct2D& Coeffs, double *rhs,
                       int time_dependent,
                       double *time_step, TFEFunction2D *old_u)
{
  int i, j, ii, N_Cells;
  int ActiveBound, N_Edges, locdof; //boundedge;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double sx, sy, meas, area, rho = 2.0; //tmp
  TBaseCell *cell, *neigh;
  FE_type CurrentElement;
  TJoint *joint;
  const TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  const int *TmpEdVer;

  coeff = new double[6];

  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->get_n_active();
  // get collection and number of cells
  auto coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the edge stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    auto DOF = fespace->GetGlobalDOF(i);

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->get_fe_type(i);
    if (CurrentElement!=C_P1_2D_T_A)
    {
      if (sold_parameter_type!=BE05_2)
      {
        ErrThrow("Edge stabilization for element ", CurrentElement,
                 " not implemented !!!");
      }
      if ((CurrentElement!=C_Q1_2D_Q_A)&&(CurrentElement!=C_Q1_2D_Q_M))
      {
        ErrThrow("Edge stabilization for element ", CurrentElement,
                 " not implemented !!!");
      }
    }
    // # of edges
    N_Edges = cell->GetN_Edges();

    sx = sy = 0;
    // compute derivatives for basis functions
    for (j=0;j<N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
      sx += x[j];
      sy += y[j];
      u->FindGradientLocal(cell, i, x[j], y[j], val);
    }
    sx /= N_Edges;
    sy /= N_Edges;
    // compute twice area of triangle
    if (N_Edges==3)
    {
      area = 2*meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }
    else
    {
      // Q1
      ErrThrow("Implementation not complete");
      area = meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }

    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    // compute gradient of current solution (constant)
    u->FindGradientLocal(cell, i, sx, sy, val);

    // compute maximal normal jump
    if (sold_parameter_type==BH04)
    {
      maxjump = 0;
      for(j=0;j<N_Edges;j++)                      // loop over all edges of cell
      {
        joint=cell->GetJoint(j);
        ver0=cell->GetVertex(TmpEdVer[2*j]);      // get vertices of face j
        ver1=cell->GetVertex(TmpEdVer[2*j+1]);
        x0 = ver0->GetX();                        // coordinates of face j
        y0 = ver0->GetY();
        x1 = ver1->GetX();
        y1 = ver1->GetY();

        // compute tangential
        t1 = x1 - x0;
        t2 = y1 - y0;
        norm_t = std::sqrt(t1*t1+t2*t2);
        t1 /= norm_t;
        t2 /= norm_t;
        // compute normal
        n1 = -t2;
        n2 = t1;
        // compute solution (including derivative) in midpoint of tangential
        // from point of view of this mesh cell
        xs = (x1+x0)/2;
        ys = (y1+y0)/2;

        // compute solution (including derivative) in midpoint of tangential
        // from point of view of neighbour mesh cell
        // NO ADAPTIVE MESHES ALLOWED
        neigh=joint->GetNeighbour(cell);          // neighbour cell
        if (neigh!=nullptr)
        {
          ii =  neigh->GetClipBoard();
          u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        }
        else
        {
          // boundary edge
          // continue;
          val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
        }
        jump = (n1 * val[1] + n2 * val[2]) - (n1 * val_neigh[1] + n2 * val_neigh[2]);
        jump = std::abs(jump);
        if (jump > maxjump)
          maxjump = jump;
      }
    }

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //Output::print("ed ", j, " ", x0, " ", y0, " ; ", x1, " ", y1);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = std::sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      //Output::print(t1, " ", t2, " ", t1*t1+t2*t2);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      //u->FindGradientLocal(cell, i, xs, ys, val);
      //Output::print("grad_i ", val[1], " ", val[2]);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of neighbour mesh cell
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=nullptr)
      {
        ii =  neigh->GetClipBoard();
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
//        boundedge = 0;
      }
      else
      {
        // boundary edge
        val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
//        boundedge = 1;
      }
      //Output::print("grad_ii ", val_neigh[1], " ", val_neigh[2]);

      // compute factor which defines the sign
      fac0 = t1 * val[1] + t2 * val[2];
      //Output::print("fac 0 ", fac0);
      /*if (fac0 > 1e-8)
          fac0 = 1;
      else
      {
          if (fac0<-1e-8)
        fac0 = -1;
          else
        fac0 = 0;
        }*/
      //Output::print("fac0 ", fac0);
//      tmp = std::tanh(fac0/1.0)/fac0;
      fac0 = std::tanh(fac0/1.0);
      //Output::print(fac0);
      // compute nonlinear factor depending on u (Psi_K(u))
      switch (sold_parameter_type)
      {
        case BH04:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, nullptr, &coeff);
          fac1 = coeff[0] * TDatabase::ParamDB->SOLD_CONST;
          if (time_dependent)
            fac1 *= time_step[0];
          fac1 += h * TDatabase::ParamDB->SOLD_S;
          fac1 *= h *maxjump/meas;
          break;
        case BE05_1:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, nullptr, &coeff);
          // norm of convection
          fac1 = std::sqrt(coeff[1]*coeff[1]+coeff[2]*coeff[2]);
          if (time_dependent)
            fac1 *= time_step[0];
          //Output::print("conv ", fac1, " ", h);
          fac1 *= h*h;
          // reaction term
          fac2 = rho * std::abs(coeff[3]);
          if (time_dependent)
          {
            fac2 *= time_step[0];
            fac2 += 1.0;
            //Output::print(fac2);
          }
          fac2 *=h*h*h;
          fac1 += fac2;
          // jump of gradient
          jump = (val[1]-val_neigh[1]) * (val[1]-val_neigh[1]);
          jump += (val[2]-val_neigh[2]) * (val[2]-val_neigh[2]);
          fac1 *= std::sqrt(jump);
          //Output::print("jump ", std::sqrt(jump));
          fac1 *= TDatabase::ParamDB->SOLD_CONST/meas;
          break;

        case BE05_2:
          // this case does not depend from formerly computed values
          // Simpson rule
          // compute coefficients in (x0,y0)
          Coeffs(1, &x0, &y0, nullptr, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x0<0.25)||(x0>0.75)||(y0<0.25)||(y0>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 = coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x0, y0, oldval);
            fac1 = val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }
          Coeffs(1, &xs, &ys, nullptr, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((sx<0.25)||(sx>0.75)||(sy<0.25)||(sy>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += 4*(coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4]);
          }
          else
          {
            old_u->FindGradientLocal(cell, i, xs, ys, oldval);
            fac1 += 4 * (val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4]);
          }
          // compute coefficients in (x1,y1)
          Coeffs(1, &x1, &y1, nullptr, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x1<0.25)||(x1>0.75)||(y1<0.25)||(y1>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x1, y1, oldval);
            fac1 += val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }

          fac1 = TDatabase::ParamDB->SOLD_CONST * std::abs(fac1)/6.0;
          break;
        default :
          ErrThrow("Edge stabilization ", sold_parameter_type,
                   " not implemented !!!");
      }
      // norm_t is the length of the edge
      fac1 = fac0*meas*fac1*norm_t;
      // update the rhs
      switch(j)
      {
        // edge zero, active dof are local 0 and 1
        case 0:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge one, active dof are local 1 and 2
        case 1:
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge two, active dof are local 0 and 2
        case 2:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
      }
    }
  }                                               // loop over cells
  /*
    // loop over all cells for computing the edge stabilization
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      // pointer to global indices of dof connected with this cell
      DOF = GlobalNumbers + BeginIndex[i];

      // local dofs are arranged as follows
      // local dof 0 on vertex 0 opposite to edge 1
  // local dof 1 on vertex 1 opposite to edge 2
  // local dof 2 on vertex 2 opposite to edge 0

  // # of edges
  N_Edges = cell->GetN_Edges();

  // compute derivatives for basis functions
  for (j=0;j<N_Edges; j++)
  {
  x[j] = cell->GetVertex(j)->GetX();
  y[j] = cell->GetVertex(j)->GetY();
  Output::print(x[j], " ", y[j], " ori ", rhsori[DOF[j]], " update ",
                rhs[DOF[j]], " diff ", rhsori[DOF[j]] - rhs[DOF[j]]);
  }
  }
  */
  delete coeff;

}

#endif // 2D


double ComputeAlpha(double hK) // copy from TCD2D.C
{
  double alpha;
  
  alpha = TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT*
     std::pow(hK, TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_POWER);
  return(alpha);

  // this is just for the response to the referee and the special example 
  // in [JKL05]
//   double b, eps, Pe, t;
// 
//   b = std::sqrt(5.0);
//   eps = 1/TDatabase::ParamDB->RE_NR;
//   Pe = b*hK/(2*eps);
//   t = 1/std::tanh(Pe) - 1/Pe;
//   alpha = t*hK/(2*b);
//   return(alpha);
}

template<>
void BilinearAssembleGalerkin<2>(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs)
{
  double val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  int j;
  
  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  const double *u = OrigValues[0];
  const double *ux = OrigValues[1];
  const double *uy = OrigValues[2];

  // coefficients
  const double c0 = coeff[0]; // eps
  const double c1 = coeff[1]; // b_1
  const double c2 = coeff[2]; // b_2
  const double c3 = coeff[3]; // c
  const double c4 = coeff[4]; // f

  for(int i=0;i<N_;i++)
  {
    // test function
    test10 = ux[i]; // xi derivative
    test01 = uy[i]; // eta derivative
    test00 = u[i]; // function
    
    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = ux[j]; // xi derivative
      ansatz01 = uy[j]; // eta derivative
      ansatz00 = u[j]; // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assemble reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;
      // lumping
      //
      //Matrix[i][i]  += Mult*c3*ansatz00 * test00;
      // quad weigth
      val *= Mult;

      // update matrix entry
      Matrix[i][j] += val;
    } // endfor j
  } // endfor i
}

template<>
void BilinearAssembleGalerkin<3>(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs)
{
  double val;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double ** Matrix = LocMatrices[0];
  double * Rhs = LocRhs[0];
  int N_ = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uz = OrigValues[3];
  double nu = coeff[0];
  double b_1 = coeff[1];
  double b_2 = coeff[2];
  double b_3 = coeff[3];
  double c = coeff[4];
  double f = coeff[5];
  for(int i = 0; i < N_; i++)
  {
    test100 = ux[i];
    test010 = uy[i];
    test001 = uz[i];
    test000 = u[i];
    Rhs[i] += Mult*test000*f;
    for(int j = 0; j < N_; j++)
    {
      ansatz100 = ux[j];
      ansatz010 = uy[j];
      ansatz001 = uz[j];
      ansatz000 = u[j];
      val = nu *(test100*ansatz100 + test010*ansatz010 + test001*ansatz001);
      val += (b_1*ansatz100 + b_2*ansatz010 + b_3*ansatz001)*test000;
      val += c*ansatz000*test000;
      Matrix[i][j] += val*Mult;
    }
  }
}

template<>
void BilinearAssemble_SD<2>(double Mult, const double *coeff, const double*,
                            double hK, const double **OrigValues,
                            const int *N_BaseFuncts, double ***LocMatrices,
                            double **LocRhs)
{
  double val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test10, test01;
  int j;
  double bgradv;

  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  const double *u = OrigValues[0];
  const double *ux = OrigValues[1];
  const double *uy = OrigValues[2];
  const double *uxx = OrigValues[3];
  const double *uyy = OrigValues[5];

  const double c0 = coeff[0]; // nu
  const double c1 = coeff[1]; // b_1
  const double c2 = coeff[2]; // b_2
  const double c3 = coeff[3]; // c
  const double c4 = coeff[4]; // f
  const double c5 = coeff[5]; // \|b\|_infty (in the entire cell)

  const double delta = Compute_SDFEM_delta<2>(hK, c0, {{c1, c2}}, c3, c5);

  for(int i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = ux[i];
    test01 = uy[i];
    //double test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += Mult*delta*bgradv*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = ux[j];
      ansatz01 = uy[j];
      ansatz00 = u[j];
      ansatz20 = uxx[j];
      ansatz02 = uyy[j];
      
      // assemble stabilization
      val = delta * ( -c0*(ansatz20 + ansatz02)
                       +c1*ansatz10 + c2*ansatz01
                       +c3*ansatz00 ) * bgradv;
      
      // quad weigth
      val *=Mult;
      
      // update matrix entry
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

template<>
void BilinearAssemble_SD<3>(double Mult, const double *coeff, const double*,
                            double hK, const double **OrigValues,
                            const int *N_BaseFuncts, double ***LocMatrices,
                            double **LocRhs)
{
  double val;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  double test100, test010, test001;
  double ** Matrix = LocMatrices[0];
  double * Rhs = LocRhs[0];
  int N_ = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uz = OrigValues[3];
  const double * uxx = OrigValues[4];
  const double * uyy = OrigValues[7];
  const double * uzz = OrigValues[9];
  double nu = coeff[0];
  double b_1 = coeff[1];
  double b_2 = coeff[2];
  double b_3 = coeff[3];
  double c = coeff[4];
  double f = coeff[5];
  double b_norm = coeff[6];
  double delta = Compute_SDFEM_delta<3>(hK, nu, {{b_1, b_2, b_3}}, c, b_norm);
  for(int i = 0; i < N_; i++)
  {
    test100 = ux[i];
    test010 = uy[i];
    test001 = uz[i];
    double bgradv = b_1*test100 + b_2*test010 + b_3*test001;
    Rhs[i] += Mult*(delta*bgradv)*f;
    for(int j = 0; j < N_; j++)
    {
      ansatz100 = ux[j];
      ansatz010 = uy[j];
      ansatz001 = uz[j];
      ansatz000 = u[j];
      ansatz200 = uxx[j];
      ansatz020 = uyy[j];
      ansatz002 = uzz[j];
      val = -nu * (ansatz200 + ansatz020 + ansatz002);
      val += b_1*ansatz100 + b_2*ansatz010 + b_3*ansatz001 + c*ansatz000;
      Matrix[i][j] += val * delta * bgradv * Mult;
    }
  }
}

template<>
void BilinearAssemble_GLS<2>(double Mult, const double *coeff, const double*,
                             double hK, const double **OrigValues,
                             const int *N_BaseFuncts, double ***LocMatrices,
                             double **LocRhs)
{
  double val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double Lu;

  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  const double * u = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uxx = OrigValues[3];
  const double * uyy = OrigValues[5];

  const double c0 = coeff[0];                                  // nu
  const double c1 = coeff[1];                                  // b_1
  const double c2 = coeff[2];                                  // b_2
  const double c3 = coeff[3];                                  // c
  const double c4 = coeff[4];                                  // f
  const double c5 = coeff[5];                                  // \|b\|_infty

  const double delta = Compute_SDFEM_delta<2>(hK, c0, {{c1, c2}}, c3, c5);

  for(int i = 0; i < N_; i++)
  {
    test10 = ux[i];
    test01 = uy[i];
    test00 = u[i];
    test20 = uxx[i];
    test02 = uyy[i];
      
    Lu = (-c0*(test20+test02)+c1*test10+c2*test01 +c3*test00);
    Rhs[i] += Mult*delta*Lu*c4;

    for(int j = 0; j < N_; j++)
    {
      ansatz10 = ux[j];
      ansatz01 = uy[j];
      ansatz00 = u[j];
      ansatz20 = uxx[j];
      ansatz02 = uyy[j];
      val = delta * (-c0*(ansatz20+ansatz02)
            +c1*ansatz10+c2*ansatz01
            +c3*ansatz00) * Lu;
      val *=Mult;
      Matrix[i][j] += val;
    }
  }
}

template<>
void BilinearAssemble_GLS<3>(double Mult, const double *coeff, const double*,
                             double hK, const double **OrigValues,
                             const int *N_BaseFuncts, double ***LocMatrices,
                             double **LocRhs)
{
  double val;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001, test200, test020, test002;
  double Lu;

  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  const double * u = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uz = OrigValues[3];
  const double * uxx = OrigValues[4];
  const double * uyy = OrigValues[7];
  const double * uzz = OrigValues[9];

  double nu = coeff[0];
  double b_1 = coeff[1];
  double b_2 = coeff[2];
  double b_3 = coeff[3];
  double c = coeff[4];
  double f = coeff[5];
  double b_norm = coeff[6];

  const double delta = Compute_SDFEM_delta<3>(hK, nu, {{b_1, b_2, b_3}}, c,
                                              b_norm);
  for(int i = 0; i < N_; i++)
  {
    test100 = ux[i];
    test010 = uy[i];
    test001 = uz[i];
    test000 = u[i];
    test200 = uxx[i];
    test020 = uyy[i];
    test002 = uzz[i];
    Lu = -nu*(test200 + test020 + test002);
    Lu += + b_1*test100 + b_2*test010 + b_3*test001 + c*test000;
    Rhs[i] += Mult*delta*Lu*f;
    for(int j = 0; j < N_; j++)
    {
      ansatz100 = ux[j];
      ansatz010 = uy[j];
      ansatz001 = uz[j];
      ansatz000 = u[j];
      ansatz200 = uxx[j];
      ansatz020 = uyy[j];
      ansatz002 = uzz[j];
      val = -nu*(ansatz200 + ansatz020 + ansatz002);
      val += b_1*ansatz100 + b_2*ansatz010 + b_3*ansatz001 + c*ansatz000;
      val *=Mult * delta * Lu;
      Matrix[i][j] += val;
    }
  }
}

template<int d>
void conv_diff_l2_h1_linf_error(int N_Points, std::array<const double*, d>,
                                const double *AbsDetjk, const double *Weights,
                                double hK, const double *const* Der,
                                const double *const* Exact,
                                const double *const* coeffs, double *LocError)
{
  LocError[0] = 0.0;
  LocError[1] = 0.0;
  LocError[2] = 0.0;
  LocError[3] = 0.0;
  LocError[4] = 0.0;
  for(int i=0;i<N_Points;i++)
  {
    double nu = coeffs[i][0];
    std::array<double, d> convection;
    for(int j = 0; j < d; ++j)
      convection[j] = coeffs[i][1+j];
    double c = coeffs[i][d+1];
    double b_norm = std::max(std::abs(convection[0]), std::abs(convection[1]));
    if(d == 3)
      b_norm = std::max(b_norm, std::abs(convection[2]));
    // SUPG parameter
    double delta = Compute_SDFEM_delta<d>(hK, nu, convection, c, b_norm);

    const double * deriv = Der[i];
    const double * exactval = Exact[i];
    double w = Weights[i]*AbsDetjk[i];

    // error in solution
    double e0 = deriv[0]-exactval[0];
    // for l2 norm
    LocError[0] += w*e0*e0;
    // for l_inf norm
    if(std::abs(e0) > LocError[4])
      LocError[4] = std::abs(e0);

    // error in x-derivative of solution
    double e1 = deriv[1]-exactval[1];
    LocError[1] += w*e1*e1;
    // error in y-derivative of solution
    double e2 = deriv[2]-exactval[2];
    LocError[1] += w*e2*e2;
    // error in z-derivative of solution
    double e3 = d == 2 ? 0. : (deriv[3] - exactval[3]);
    LocError[1] += w*e3*e3;

    // sd error
    // THIS IS ONLY CORRECT IF DIV(b) = 0
    double e4 = convection[0]*e1 + convection[1]*e2;
    if(d == 3)
      e4 += convection[2]*e3;
    LocError[2] += w*(nu*(e1*e1 + e2*e2 + e3*e3) + c*e0*e0 + delta*e4*e4);

    // DG error
    // THIS IS ONLY THE VOLUME PART
      LocError[3] += w*(nu*(e1*e1 + e2*e2 + e3*e3)  // diffusion_coef x H1-semi
          + (e4  // convection * nabla(eh) x eh
            + c*e0)*e0);  // reaction x L2
  } // endfor i
}

template<>
void PoissonAssemble<2>(double Mult, const double *coeff,
                        const double*, double,
                        const double **OrigValues,
                        const int *N_BaseFuncts, double ***LocMatrices,
                        double **LocRhs)
{
  double **MatrixA, *Rhs, *MatrixRowA;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  const double *Orig0, *Orig1, *Orig2;
  
  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  int N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; //u
  Orig1 = OrigValues[1]; //ux
  Orig2 = OrigValues[2]; //uy
  
  double c4 = coeff[4]; // f
  
  for(int i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test00 = Orig0[i]; //u
    test10 = Orig1[i]; //ux
    test01 = Orig2[i]; //uy

    Rhs[i] += Mult*test00*c4;

    for(int j=0;j<N_;j++)
    {
      ansatz10 = Orig1[j];
      ansatz01 = Orig2[j];      

      MatrixRowA[j] += Mult * (test10*ansatz10+test01*ansatz01);                
    } // endfor j
  } // endfor i
}

template<>
void PoissonAssemble<3>(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs)
{
  double val;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double ** Matrix = LocMatrices[0];
  double * Rhs = LocRhs[0];
  int N_ = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uz = OrigValues[3];
  double nu = coeff[0];
  double b_1 = coeff[1];
  double b_2 = coeff[2];
  double b_3 = coeff[3];
  double c = coeff[4];
  double f = coeff[5];
  for(int i = 0; i < N_; i++)
  {
    test100 = ux[i];
    test010 = uy[i];
    test001 = uz[i];
    test000 = u[i];
    Rhs[i] += Mult*test000*f;
    for(int j = 0; j < N_; j++)
    {
      ansatz100 = ux[j];
      ansatz010 = uy[j];
      ansatz001 = uz[j];
      ansatz000 = u[j];
      val = test100*ansatz100 + test010*ansatz010 + test001*ansatz001;
      Matrix[i][j] += val*Mult;
    }
  }
}
template<>
void AssembleCDRHomogeneous<2>(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs)
{
  double **MatrixA, *Rhs, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  const double *Orig0, *Orig1, *Orig2;
  double c0, c1, c2, c3;
  
  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c

  for(int i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += 0;

    for(int j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      MatrixRowA[j] += Mult * val;                
    } // endfor j
  } // endfor i
}

template<>
void AssembleCDRHomogeneous<3>(double , const double *,
                                 const double*, double,
                                 const double **,
                                 const int *, double ***,
                                 double **)
{
  ErrThrow("AssembleCDRHomogeneous 3D case is not supported so far");
}

#ifdef __2D__
template void conv_diff_l2_h1_linf_error<2>(
  int N_Points, std::array<const double*, 2> xyz, const double *AbsDetjk,
  const double *Weights, double hK, const double *const*Der,
  const double *const*Exact, const double *const*coeffs, double *LocError);
#endif // 2D
#ifdef __3D__
template void conv_diff_l2_h1_linf_error<3>(
  int N_Points, std::array<const double*, 3> xyz, const double *AbsDetjk,
  const double *Weights, double hK, const double *const*Der,
  const double *const*Exact, const double *const*coeffs, double *LocError);
#endif // 3D
