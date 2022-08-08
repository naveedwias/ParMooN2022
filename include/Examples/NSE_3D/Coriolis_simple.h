/**
 * @file A Navier--Stokes test problem with coriolis force.
 *
 * The boundary data is adapted to the [0,1]^3 unit cube example, which is
 * availabe as default geometry in ParMooN. It will throw an error if you
 * try running it on any other domain - just to make you aware of that fact.

 */
#include <cmath>

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;
double omega; // coriolis force is the vector (0,0,omega)
// switch between Stokes (false) and Navier-Stokes (true)
bool include_nonlinear_term = true;
// turn on or off the coriolis force
bool include_coriolis_term = true;

// example specific parameters:
double alpha; // angle, interesting values are 0, 1/20, pi/2-1/20, pi/2
double v_0; // velocity scaling
double K = 1.;
double R = 4.;
// switch between different analytic solutions
int example = 2;


void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","Simple example with coriolis force");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  switch(example)
  {
    case 0:
      values[0] = -v_0 * std::cos(alpha) * y;
      values[1] = 0.;
      values[2] = -v_0 * std::cos(alpha);
      values[3] = 0.;
      values[4] = 0.; //Laplacien
      break;
    case 1:
    {
      values[0] = y*(K*std::pow(x, 4) - 6*K*std::pow(x, 2)*std::pow(y, 2)
                     + 12*K*std::pow(x, 2)*std::pow(z, 2) + K*std::pow(y, 4) 
                     - 4*K*std::pow(y, 2)*std::pow(z, 2) - omega);
      values[1] = 4*K*x*y*(std::pow(x, 2) - 3*std::pow(y, 2) 
                  + 6*std::pow(z, 2));
      values[2] = K*std::pow(x, 4) - 18*K*std::pow(x, 2)*std::pow(y, 2)
                  + 12*K*std::pow(x, 2)*std::pow(z, 2) + 5*K*std::pow(y, 4)
                  - 12*K*std::pow(y, 2)*std::pow(z, 2) - omega;
      values[3] = 8*K*y*z*(3*std::pow(x, 2) - std::pow(y, 2));
      values[4] = 0.;
      break;
    }
    case 2:
    {
      double xx = x*x;
      double yy = y*y;
      double zz = z*z;
      values[0] = std::pow(xx + yy + zz, -0.5*R)
                  *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                    - K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + K*xx*y*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + K*y*y*y*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    - omega*xx*y*std::pow(xx + yy + zz, 0.5*R)
                    - omega*y*y*y*std::pow(xx + yy + zz, 0.5*R))
                  /(xx + yy);
      if(values[0] != values[0])
        ErrThrow("evaluating exact solution at (", x, ",", y, ",", z, ")");
      values[1] = K*R*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, -0.5*R)
                  *(-R*std::pow(x, 3)*y*zz*std::cos(R*std::atan2(y, x))
                    - R*xx*yy*zz*std::sin(R*std::atan2(y, x))
                    + R*xx*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    - R*x*std::pow(y, 3)*zz*std::cos(R*std::atan2(y, x))
                    - 2*R*x*y*std::pow(z, 4)*std::cos(R*std::atan2(y, x))
                    - R*std::pow(y, 4)*zz*std::sin(R*std::atan2(y, x))
                    - R*yy*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    + std::pow(x, 4)*yy*std::sin(R*std::atan2(y, x))
                    - std::pow(x, 4)*zz*std::sin(R*std::atan2(y, x))
                    + 3*std::pow(x, 3)*y*zz*std::cos(R*std::atan2(y, x))
                    + 2*xx*std::pow(y, 4)*std::sin(R*std::atan2(y, x))
                    + xx*yy*zz*std::sin(R*std::atan2(y, x))
                    - xx*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    + 3*x*std::pow(y, 3)*zz*std::cos(R*std::atan2(y, x))
                    + 2*x*y*std::pow(z, 4)*std::cos(R*std::atan2(y, x))
                    + std::pow(y, 6)*std::sin(R*std::atan2(y, x))
                    + 2*std::pow(y, 4)*zz*std::sin(R*std::atan2(y, x))
                    + yy*std::pow(z, 4)*std::sin(R*std::atan2(y, x)))
                  /(std::pow(x, 6) + 3*std::pow(x, 4)*yy + std::pow(x, 4)*zz
                    + 3*xx*std::pow(y, 4) + 2*xx*yy*zz + std::pow(y, 6)
                    + std::pow(y, 4)*zz);
      values[2] = -std::pow(xx + yy + zz, -1.5L*R - 2)
                  *(R*y*(xx + yy)*std::pow(xx + yy + zz, R + 1)
                    *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                      - K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                      + K*xx*y*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                      + K*std::pow(y, 3)*std::pow(xx + yy, 0.5*R)
                        *std::cos(R*std::atan2(y, x))
                      - omega*xx*y*std::pow(xx + yy + zz, 0.5*R)
                      - omega*std::pow(y, 3)*std::pow(xx + yy + zz, 0.5*R))
                    + 2*y*std::pow(xx + yy + zz, R + 2)
                      *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                        - K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                        + K*xx*y*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                        + K*std::pow(y, 3)*std::pow(xx + yy, 0.5*R)
                          *std::cos(R*std::atan2(y, x))
                        - omega*xx*y*std::pow(xx + yy + zz, 0.5*R)
                        - omega*std::pow(y, 3)*std::pow(xx + yy + zz, 0.5*R))
                      + std::pow(xx + yy + zz, R + 1)
                        *(-K*R*std::pow(xx + yy, 0.5*R)*(xx + yy + zz)
                          *(R*xx*zz*std::cos(R*std::atan2(y, x))
                          + 2*R*x*y*zz*std::sin(R*std::atan2(y, x))
                          - R*yy*zz*std::cos(R*std::atan2(y, x))
                          - std::pow(x, 3)*y*std::sin(R*std::atan2(y, x))
                          + xx*yy*std::cos(R*std::atan2(y, x))
                          - x*std::pow(y, 3)*std::sin(R*std::atan2(y, x))
                          + std::pow(y, 4)*std::cos(R*std::atan2(y, x)))
                        + R*omega*yy*std::pow(xx + yy, 2)
                          *std::pow(xx + yy + zz, 0.5*R)
                        + (xx + yy)*(xx + yy + zz)
                          *(K*R*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            - K*xx*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            - 3*K*yy*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            + omega*xx*std::pow(xx + yy + zz, 0.5*R)
                            + 3*omega*yy*std::pow(xx + yy + zz, 0.5*R))))
                    /std::pow(xx + yy, 2);
      values[3] = K*R*z*std::pow(xx + yy, 0.5*R - 1)
                  *std::pow(xx + yy + zz, -0.5L*R - 1)
                  *(-R*x*zz*std::sin(R*std::atan2(y, x))
                  + R*y*zz*std::cos(R*std::atan2(y, x)) - xx*y*std::cos(R*std::atan2(y, x))
                  - std::pow(y, 3)*std::cos(R*std::atan2(y, x))
                  - 2*(-x*std::sin(R*std::atan2(y, x)) + y*std::cos(R*std::atan2(y, x)))
                    *(xx + yy + zz));
      values[4] = -K*R*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, -0.5*R)
                   *(std::pow(R, 2)*x*zz*std::sin(R*std::atan2(y, x))
                    - std::pow(R, 2)*y*zz*std::cos(R*std::atan2(y, x))
                    + R*xx*y*std::cos(R*std::atan2(y, x))
                    + 3*R*x*zz*std::sin(R*std::atan2(y, x))
                    + R*std::pow(y, 3)*std::cos(R*std::atan2(y, x))
                    - 3*R*y*zz*std::cos(R*std::atan2(y, x))
                    + 3*xx*y*std::cos(R*std::atan2(y, x))
                    + 3*std::pow(y, 3)*std::cos(R*std::atan2(y, x)))
                   /(std::pow(x, 4) + 2*xx*yy + xx*zz + std::pow(y, 4) + yy*zz);
      break;
    }
    default:
      ErrThrow("unknown analytic solution");
  }
}
void ExactU2(double x, double y,  double z, double *values)
{
  switch(example)
  {
    case 0:
      values[0] = v_0 * (std::cos(alpha) * x + std::sin(alpha) * z);
      values[1] = v_0 * std::cos(alpha);
      values[2] = 0.;
      values[3] = v_0 * std::sin(alpha);
      values[4] = 0.; //Laplacien
      break;
    case 1:
    {
      values[0] = x*(-K*std::pow(x, 4) + 6*K*std::pow(x, 2)*std::pow(y, 2) 
                     + 4*K*std::pow(x, 2)*std::pow(z, 2) - K*std::pow(y, 4) 
                     - 12*K*std::pow(y, 2)*std::pow(z, 2) + omega);
      values[1] = -5*K*std::pow(x, 4) + 18*K*std::pow(x, 2)*std::pow(y, 2) 
                  + 12*K*std::pow(x, 2)*std::pow(z, 2) - K*std::pow(y, 4) 
                  - 12*K*std::pow(y, 2)*std::pow(z, 2) + omega;
      values[2] = 4*K*x*y*(3*std::pow(x, 2) - std::pow(y, 2)
                  - 6*std::pow(z, 2));
      values[3] = 8*K*x*z*(std::pow(x, 2) - 3*std::pow(y, 2));
      values[4] = 0.;
      break;
    }
    case 2:
    {
      double xx = x*x;
      double yy = y*y;
      double zz = z*z;
      values[0] = std::pow(xx + yy + zz, -0.5*R)
                  *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                    - K*std::pow(x, 3)*std::pow(xx + yy, 0.5*R)
                      *std::cos(R*std::atan2(y, x))
                    - K*x*yy*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + omega*std::pow(x, 3)*std::pow(xx + yy + zz, 0.5*R)
                    + omega*x*yy*std::pow(xx + yy + zz, 0.5*R))
                  /(xx + yy);
      values[1] = std::pow(xx + yy + zz, -3.0L/2.0L*R - 2)
                  *(-R*x*(xx + yy)*std::pow(xx + yy + zz, R + 1)
                  *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                    - K*std::pow(x, 3)*std::pow(xx + yy, 0.5*R)
                      *std::cos(R*std::atan2(y, x))
                    - K*x*yy*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                    + omega*std::pow(x, 3)*std::pow(xx + yy + zz, 0.5*R)
                    + omega*x*yy*std::pow(xx + yy + zz, 0.5*R))
                  - 2*x*std::pow(xx + yy + zz, R + 2)
                    *(K*R*x*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                      + K*R*y*zz*std::pow(xx + yy, 0.5*R)*std::sin(R*std::atan2(y, x))
                      - K*std::pow(x, 3)*std::pow(xx + yy, 0.5*R)
                        *std::cos(R*std::atan2(y, x))
                      - K*x*yy*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                      + omega*std::pow(x, 3)*std::pow(xx + yy + zz, 0.5*R)
                      + omega*x*yy*std::pow(xx + yy + zz, 0.5*R))
                    + std::pow(xx + yy + zz, R + 1)
                      *(K*R*std::pow(xx + yy, 0.5*R)*(xx + yy + zz)
                        *(R*xx*zz*std::cos(R*std::atan2(y, x))
                          + 2*R*x*y*zz*std::sin(R*std::atan2(y, x))
                          - R*yy*zz*std::cos(R*std::atan2(y, x))
                          - std::pow(x, 4)*std::cos(R*std::atan2(y, x))
                          - std::pow(x, 3)*y*std::sin(R*std::atan2(y, x))
                          - xx*yy*std::cos(R*std::atan2(y, x))
                          - x*std::pow(y, 3)*std::sin(R*std::atan2(y, x)))
                        + R*omega*xx*std::pow(xx + yy, 2)
                          *std::pow(xx + yy + zz, 0.5*R)
                        + (xx + yy)*(xx + yy + zz)
                          *(K*R*zz*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            - 3*K*xx*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            - K*yy*std::pow(xx + yy, 0.5*R)*std::cos(R*std::atan2(y, x))
                            + 3*omega*xx*std::pow(xx + yy + zz, 0.5*R)
                            + omega*yy*std::pow(xx + yy + zz, 0.5*R))))
                  /std::pow(xx + yy, 2);
      values[2] = K*R*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, -0.5*R)
                  *(-R*std::pow(x, 4)*zz*std::sin(R*std::atan2(y, x))
                    + R*std::pow(x, 3)*y*zz*std::cos(R*std::atan2(y, x))
                    - R*xx*yy*zz*std::sin(R*std::atan2(y, x))
                    - R*xx*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    + R*x*std::pow(y, 3)*zz*std::cos(R*std::atan2(y, x))
                    + 2*R*x*y*std::pow(z, 4)*std::cos(R*std::atan2(y, x))
                    + R*yy*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    + std::pow(x, 6)*std::sin(R*std::atan2(y, x))
                    + 2*std::pow(x, 4)*yy*std::sin(R*std::atan2(y, x))
                    + 2*std::pow(x, 4)*zz*std::sin(R*std::atan2(y, x))
                    - 3*std::pow(x, 3)*y*zz*std::cos(R*std::atan2(y, x))
                    + xx*std::pow(y, 4)*std::sin(R*std::atan2(y, x))
                    + xx*yy*zz*std::sin(R*std::atan2(y, x))
                    + xx*std::pow(z, 4)*std::sin(R*std::atan2(y, x))
                    - 3*x*std::pow(y, 3)*zz*std::cos(R*std::atan2(y, x))
                    - 2*x*y*std::pow(z, 4)*std::cos(R*std::atan2(y, x))
                    - std::pow(y, 4)*zz*std::sin(R*std::atan2(y, x))
                    - yy*std::pow(z, 4)*std::sin(R*std::atan2(y, x)))
                  /(std::pow(x, 6) + 3*std::pow(x, 4)*yy + std::pow(x, 4)*zz
                    + 3*xx*std::pow(y, 4) + 2*xx*yy*zz + std::pow(y, 6)
                    + std::pow(y, 4)*zz);
      values[3] = K*R*z*std::pow(xx + yy, 0.5*R - 1)
                  *std::pow(xx + yy + zz, -0.5*R - 1)
                  *(-R*x*zz*std::cos(R*std::atan2(y, x))
                    - R*y*zz*std::sin(R*std::atan2(y, x))
                    + std::pow(x, 3)*std::cos(R*std::atan2(y, x))
                    + x*yy*std::cos(R*std::atan2(y, x)) + 2*(x*std::cos(R*std::atan2(y, x))
                    + y*std::sin(R*std::atan2(y, x)))*(xx + yy + zz));
      values[4] = K*R*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, -0.5*R)
                  *(-std::pow(R, 2)*x*zz*std::cos(R*std::atan2(y, x))
                    - std::pow(R, 2)*y*zz*std::sin(R*std::atan2(y, x))
                    + R*std::pow(x, 3)*std::cos(R*std::atan2(y, x))
                    + R*x*yy*std::cos(R*std::atan2(y, x)) - 3*R*x*zz*std::cos(R*std::atan2(y, x))
                    - 3*R*y*zz*std::sin(R*std::atan2(y, x))
                    + 3*std::pow(x, 3)*std::cos(R*std::atan2(y, x))
                    + 3*x*yy*std::cos(R*std::atan2(y, x)))
                  /(std::pow(x, 4) + 2*xx*yy + xx*zz + std::pow(y, 4) + yy*zz);
      break;
    }
    default:
      ErrThrow("unknown analytic solution");
  }
}
void ExactU3(double x, double y,  double z, double *values)
{
  switch(example)
  {
    case 0:
      values[0] = -v_0 * std::sin(alpha) * y;
      values[1] = 0.;
      values[2] = -v_0 * std::sin(alpha);
      values[3] = 0.;
      values[4] = 0.; //Laplacien
      break;
    case 1:
    {
      values[0] = 16*K*x*y*z*(-std::pow(x, 2) + std::pow(y, 2));
      values[1] = 16*K*y*z*(-3*std::pow(x, 2) + std::pow(y, 2));
      values[2] = 16*K*x*z*(-std::pow(x, 2) + 3*std::pow(y, 2));
      values[3] = 16*K*x*y*(-std::pow(x, 2) + std::pow(y, 2));
      values[4] = 0.;
      break;
    }
    case 2:
    {
      double xx = x*x;
      double yy = y*y;
      double zz = z*z;
      values[0] = -K*R*z*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, -0.5*R)
                   *std::sin(R*std::atan2(y, x));
      values[1] = K*std::pow(R, 2)*z
                  *(x*std::pow(xx + yy, 0.5*R + 1)*std::pow(xx + yy + zz, 0.5*R)
                    *std::sin(R*std::atan2(y, x))
                    + std::pow(xx + yy, 0.5*R)*(-x*std::sin(R*std::atan2(y, x))
                    + y*std::cos(R*std::atan2(y, x)))*std::pow(xx + yy + zz, 0.5*R
                    + 1))
                  *std::pow(xx + yy + zz, -R - 1)/(xx + yy);
      values[2] = K*std::pow(R, 2)*z
                  *(y*std::pow(xx + yy, 0.5*R + 1)*std::pow(xx + yy + zz, 0.5*R)
                    *std::sin(R*std::atan2(y, x))
                    - std::pow(xx + yy, 0.5*R)
                      *(x*std::cos(R*std::atan2(y, x)) + y*std::sin(R*std::atan2(y, x)))
                      *std::pow(xx + yy + zz, 0.5*R + 1)
                   )*std::pow(xx + yy + zz, -R - 1)/(xx + yy);
      values[3] = -K*R*std::pow(xx + yy, 0.5*R)
                   *std::pow(xx + yy + zz, -0.5*R - 1)
                   *(-R*zz + xx + yy + zz)*std::sin(R*std::atan2(y, x));
      values[4] = K*std::pow(R, 2)*z*std::pow(xx + yy + zz, -R)
                  *(R*xx*std::pow(xx + yy, 0.5*R)*std::pow(xx + yy + zz, 0.5*R)
                    + R*yy*std::pow(xx + yy, 0.5*R)
                      *std::pow(xx + yy + zz, 0.5*R)
                    + 2*R*zz*std::pow(xx + yy, 0.5*R)
                      *std::pow(xx + yy + zz, 0.5*R)
                    - R*zz*std::pow(std::pow(x, 4) + 2*xx*yy + xx*zz 
                                    + std::pow(y, 4) + yy*zz, 0.5*R)
                    + 3*xx*std::pow(std::pow(x, 4) + 2*xx*yy + xx*zz
                                    + std::pow(y, 4) + yy*zz, 0.5*R)
                    + 3*yy*std::pow(std::pow(x, 4) + 2*xx*yy + xx*zz
                                    + std::pow(y, 4) + yy*zz, 0.5*R)
                    + 2*zz*std::pow(xx + yy, 0.5*R)
                      *std::pow(xx + yy + zz, 0.5*R)
                    + zz*std::pow(std::pow(x, 4) + 2*xx*yy + xx*zz
                                  + std::pow(y, 4) + yy*zz, 0.5*R))
                  *std::sin(R*std::atan2(y, x))/(std::pow(x, 4)
                  + 2*xx*yy + 2*xx*zz + std::pow(y, 4) + 2*yy*zz
                  + std::pow(z, 4));
      break;
    }
    default:
      ErrThrow("unknown analytic solution");
  }
}

void ExactP(double x, double y,  double z, double *values)
{
  switch(example)
  {
    case 0:
    {
      double cosa = std::cos(alpha);
      double sina = std::sin(alpha);
      values[0] = 0.5 * v_0 * v_0 * (y*y + cosa*cosa*x*x + 2*x*z*cosa*sina 
                  + z*z*sina*sina) + omega * v_0 * (x*x + y*y) * cosa 
                  + 0.5 * omega*omega * (x*x + y*y);
      values[1] = v_0 * v_0 * (cosa*cosa*x + z*cosa*sina)
                  + omega * v_0 * 2*x * cosa + omega*omega * x;
      values[2] = v_0 * v_0 * y + omega * v_0 * 2*y * cosa + omega*omega * y;
      values[3] = v_0 * v_0 * (x*cosa*sina + z*sina*sina);
      values[4] = v_0 * v_0 * (1 + cosa*cosa + sina*sina)
                  + 4 * omega * v_0 * cosa + 2 * omega*omega;
      break;
    }
    case 1:
    case 2:
    {
      values[0] = 0.;
      values[1] = 0.;
      values[2] = 0.;
      values[3] = 0.;
      values[4] = 0.;
      break;
    }
    default:
      ErrThrow("unknown analytic solution");
  }
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}
// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU1(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U2BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU2(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U3BoundValue(int, double x, double y, double z, double &value)
{
  double diri[5];
  ExactU3(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}

void LinCoeffs(int n_points, const double * X, const double * Y,
               const double * Z, const double *const*, double **coeffs)
{
  const double eps = DIMENSIONLESS_VISCOSITY;
  double u1[5], u2[5], u3[5], p[5];
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = eps;

    ExactU1(X[i], Y[i], Z[i], u1);
    ExactU2(X[i], Y[i], Z[i], u2);
    ExactU3(X[i], Y[i], Z[i], u3);
    ExactP( X[i], Y[i], Z[i], p);

    coeffs[i][1] = -eps * u1[4] + p[1]; //Stokes: diffusion and pressure gradient
    coeffs[i][2] = -eps * u2[4] + p[2];
    coeffs[i][3] = -eps * u3[4] + p[3];
    if(include_nonlinear_term) // Navier--Stokes:
    {
      // add convective terms
      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3];
      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3];
      coeffs[i][3] += u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3];
    }
    if(include_coriolis_term)
    {
      coeffs[i][1] -= omega * u2[0];
      coeffs[i][2] += omega * u1[0];
    }
    coeffs[i][4] = u1[1] + u2[2] + u3[3]; // div(u) = g
    coeffs[i][5] = 0.; //sigma
    coeffs[i][6] = 0.; // first entry in coriolis force
    coeffs[i][7] = 0.; // second entry in coriolis force
    coeffs[i][8] = omega; // third entry in coriolis force
  }
}

