// ***********************************************************************
// Q1 element with space dependent exponential bubble, conforming, 2D
// *value : IN = space data, OUT = basis values
// Author : Sashi
// History: 26.05.2011 
// ***********************************************************************

// base function values
static void C_Q_EL1_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t6, t9,t12, t15, t18, t21, t41, t48, t57, t64, t65, t67, t71 ;
  double P, Q;
  
  P = values[0];
  Q = values[1];
  
  t1 = 1.0-xi;
  t2 = 1.0-eta;
  t6 = std::exp(-0.1577350269e1*P);
  t9 = std::exp(-2.0*P);
  t12 = std::exp(-0.1577350269e1*Q);
  t15 = std::exp(-2.0*Q);
  t18 = std::exp(-0.4226497308*P);
  t21 = std::exp(-0.4226497308*Q);
  t41 = -4.0*t6+0.4000000001e1*t9-4.0*t12+0.4000000001e1*t15
            -4.0*t18-4.0*t21+0.4000000001e1+0.4000000001e1*t9*t15
            -4.0*t9*t12-4.0*t6*t15+4.0*t6*t12-4.0*t18*t15+4.0*t18*t12
            -4.0*t9*t21+4.0*t6*t21+4.0*t18*t21;
  t48 = std::exp(-P*t1);
  t57 = std::exp(-Q*t2);
  t64 = 1/t41*(-1.0+t15)*(-1.0+t9)*(1.0/2.0+xi/2.0
            -(t48-t9)/(1.0-t9))*(1.0/2.0+eta/2.0-(t57-t15)/(1.0-t15));
  t65 = 0.4e1*t64;
  t67 = 1.0+xi;
  t71 = 1.0+eta;
      
   values[0] = t1*t2/4.0-t65;
   values[1] = t67*t2/4.0-t65;
   values[2] = t1*t71/4.0-t65;
   values[3] = t67*t71/4.0-t65;
   values[4] = 4.0*t64;
}

// values of the derivatives in xi direction
static void C_Q_EL1_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t3, t6, t9, t12, t15, t18, t38, t45, t55, t62, t63;
  double P, Q;
  
  P = values[0];
  Q = values[1];
  
  t1 = eta/4.0;
  t3 = std::exp(-0.1577350269e1*P);
  t6 = std::exp(-2.0*P);
  t9 = std::exp(-0.1577350269e1*Q);
  t12 = std::exp(-2.0*Q);
  t15 = std::exp(-0.4226497308*P);
  t18 = std::exp(-0.4226497308*Q);
  t38 = -4.0*t3+0.4000000001e1*t6-4.0*t9+0.4000000001e1*t12
             -4.0*t15-4.0*t18+0.4000000001e1+0.4000000001e1*t6*t12
             -4.0*t6*t9-4.0*t3*t12+4.0*t3*t9-4.0*t15*t12+4.0*t9*t15
             -4.0*t6*t18+4.0*t3*t18+4.0*t18*t15;
  t45 = std::exp(-P*(1.0-xi));
  t55 = std::exp(-Q*(1.0-eta));
  t62 = 1/t38*(-1.0+t12)*(-1.0+t6)*(1.0/2.0-P*t45/(1.0-t6))*(1.0/2.0+eta/2.0-(t55-t12)/(1.0-t12));
  t63 = 0.4e1*t62;
      
  values[0] = -1.0/4.0+t1-t63;
  values[1] = 1.0/4.0-t1-t63;
  values[2] = -1.0/4.0-t1-t63;
  values[3] = 1.0/4.0+t1-t63;
  values[4] = 4.0*t62;
}

// values of the derivatives in eta direction
static void C_Q_EL1_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t3, t6, t9, t12, t15, t18, t38, t46, t55, t62, t63;
  double P, Q;
  
  P = values[0];
  Q = values[1];
  
  t1 = xi/4.0;
  t3 = std::exp(-0.1577350269e1*P);
  t6 = std::exp(-2.0*P);
  t9 = std::exp(-0.1577350269e1*Q);
  t12 = std::exp(-2.0*Q);
  t15 = std::exp(-0.4226497308*P);
  t18 = std::exp(-0.4226497308*Q);
  t38 = -4.0*t3+0.4000000001e1*t6-4.0*t9+0.4000000001e1*t12
        -4.0*t15-4.0*t18+0.4000000001e1+0.4000000001e1*t6*t12
        -4.0*t6*t9-4.0*t3*t12+4.0*t3*t9-4.0*t15*t12+4.0*t9*t15
        -4.0*t6*t18+4.0*t3*t18+4.0*t18*t15;
  t46 = std::exp(-P*(1.0-xi));
  t55 = std::exp(-Q*(1.0-eta));
  t62 = 1/t38*(-1.0+t12)*(-1.0+t6)*(1.0/2.0+xi/2.0-(t46-t6)/(1.0-t6))*(1.0/2.0-Q*t55/(1.0-t12));
  t63 = 0.4e1*t62;
  
  values[0] = -1.0/4.0+t1-t63;
  values[1] = -1.0/4.0-t1-t63;
  values[2] = 1.0/4.0-t1-t63;
  values[3] = 1.0/4.0+t1-t63;
  values[4] = 4.0*t62;
}

// values of the derivatives in xi-xi  direction
static void C_Q_EL1_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t2, t5,t8, t11, t14, t17, t37, t43, t46, t53, t61, t62;
  double P, Q;
  
  P = values[0];
  Q = values[1];
  
  t2 = std::exp(-0.1577350269e1*P);
  t5 = std::exp(-2.0*P);
  t8 = std::exp(-0.1577350269e1*Q);
  t11 = std::exp(-2.0*Q);
  t14 = std::exp(-0.4226497308*P);
  t17 = std::exp(-0.4226497308*Q);
  t37 = -4.0*t2+0.4000000001e1*t5-4.0*t8+0.4000000001e1*t11
        -4.0*t14-4.0*t17+0.4000000001e1+0.4000000001e1*t5*t11
        -4.0*t5*t8-4.0*t2*t11+4.0*t2*t8-4.0*t14*t11+4.0*t14*t8
        -4.0*t5*t17+4.0*t2*t17+4.0*t14*t17;
  t43 = P*P;
  t46 = std::exp(-P*(1.0-xi));
  t53 = std::exp(-Q*(1.0-eta));
  t61 = 1/t37*(-1.0+t11)*(-1.0+t5)*t43*t46/(1.0-t5)*(1.0/2.0+eta/2.0-(t53-t11)/(1.0-t11));
  t62 = 0.4e1*t61;
  
  values[0] = t62;
  values[1] = t62;
  values[2] = t62;
  values[3] = t62;
  values[4] = -4.0*t61; 
}

// values of the derivatives in xi-eta direction
static void C_Q_EL1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t2, t5, t8, t11, t14, t17, t37, t44, t53, t60, t61, t62, t63;
  double P, Q;
  
  P = values[0];
  Q = values[1];
  
  t2 = std::exp(-0.1577350269e1*P);
  t5 = std::exp(-2.0*P);
  t8 = std::exp(-0.1577350269e1*Q);
  t11 = std::exp(-2.0*Q);
  t14 = std::exp(-0.4226497308*P);
  t17 = std::exp(-0.4226497308*Q);
  t37 = -4.0*t2+0.4000000001e1*t5-4.0*t8+0.4000000001e1*t11
        -4.0*t14-4.0*t17+0.4000000001e1+0.4000000001e1*t5*t11
        -4.0*t5*t8-4.0*t2*t11+4.0*t2*t8-4.0*t14*t11+4.0*t14*t8
        -4.0*t5*t17+4.0*t2*t17+4.0*t14*t17;
  t44 = std::exp(-P*(1.0-xi));
  t53 = std::exp(-Q*(1.0-eta));
  t60 = 1/t37*(-1.0+t11)*(-1.0+t5)*(1.0/2.0-P*t44/(1.0-t5))*(1.0/2.0-Q*t53/(1.0-t11));
  t61 = 0.4e1*t60;
  t62 = 1.0/4.0-t61;
  t63 = -1.0/4.0-t61;
      
  values[0] = t62;
  values[1] = t63;
  values[2] = t63;
  values[3] = t62;
  values[4] = 4.0*t60; 
}

// values of the derivatives in eta-eta direction
static void C_Q_EL1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t2, t5, t8, t11, t14, t17, t37, t46, t52, t56, t61, t62;
  double P, Q;
  
  P = values[0];
  Q = values[1];

  t2 = std::exp(-0.1577350269e1*P);
  t5 = std::exp(-2.0*P);
  t8 = std::exp(-0.1577350269e1*Q);
  t11 = std::exp(-2.0*Q);
  t14 = std::exp(-0.4226497308*P);
  t17 = std::exp(-0.4226497308*Q);
  t37 = -4.0*t2+0.4000000001e1*t5-4.0*t8+0.4000000001e1*t11-4.0*t14
        -4.0*t17+0.4000000001e1+0.4000000001e1*t5*t11-4.0*t5*t8
        -4.0*t2*t11+4.0*t2*t8-4.0*t14*t11+4.0*t14*t8-4.0*t5*t17
        +4.0*t2*t17+4.0*t14*t17;
  t46 = std::exp(-P*(1.0-xi));
  t52 = Q*Q;
  t56 = std::exp(-Q*(1.0-eta));
  t61 = 1/t37*(-1.0+t11)*(-1.0+t5)*(1.0/2.0+xi/2.0-(t46-t5)/(1.0-t5))*t52*t56/(1.0-t11);
  t62 = 0.4e1*t61;
  
  values[0] = t62;
  values[1] = t62;
  values[2] = t62;
  values[3] = t62;
  values[4] = -4.0*t61; 
}
