// Brinkman problem, solution in ansatz space
// Poiseuille (exact solution in P2/P1)
// 
// u(x,y,z) = ( 3*(1-r^2), 0 , 0)^T=(u1,u2,u3)^T
// r^2:= y^2+z^2
// p(x,y,z) = -x+1/2

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double, double, double *values)
{
    values[0] = 0;           // u1
    values[1] = 0;           // u1_x
    values[2] = 0;           // u1_y
    values[3] = 0;           // u1_z
    values[4] = 0;           // Delta u1
}

void ExactU2(double, double, double, double *values)
{
  values[0] = 0;            // u2
  values[1] = 0;            // u2_x
  values[2] = 0;            // u2_y
  values[3] = 0;            // u2_z
  values[4] = 0;            // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU3(double x, double y, double, double *values)
{
    values[0] = 1-(x*x+y*y);//3*(1-(y*y+x*x));      // u3
    values[1] = -2*x; //-6*x*x;        // u3_x
    values[2] = -2*y;        // u3_y
    values[3] = 0;        // u3_z
    values[4] = -4; //-12;        // Delta u3=u3_xx+u3_yy+u3_zz

}

void ExactP(double, double, double z, double *values)
{
  values[0] = z-8;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double, double, double, double &value)
{
    value = 0;
}

// value of boundary condition
void U2BoundValue(int, double, double, double, double &value)
{
    value = 0;
}

// value of boundary condition
void U3BoundValue(int, double x, double y, double, double &value)
{
    value = 1-(x*x+y*y);//3*(1-(y*y+x*x))
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
    double *coeff;
    
    for(int i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        
        coeff[0] = 1;//eps;
        coeff[1] = 0; // (1/coeff[0])*(-12)+(-1)+(coeff[5])*(3*(1-(Y[i]*Y[i]+Z[i]*Z[i]))); // f1
        coeff[2] = 0; // f2
        coeff[3] = 0.;//-(1/coeff[0])*(-4)+(1)+(coeff[5])*(1-(Y[i]*Y[i]+X[i]*X[i])); // 0; // f3
        coeff[4] = 0; // g
        coeff[5] = 0; // sigma
    }
}
