// First order Raviart-Thomas vector element on quads, nonconforming, 2D

static void N_Q_RT0_2D_Funct(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.25*(xi+1.);
  values[2]= 0.;
  values[3]= 0.25*(xi-1.);

  // second component
  values[4]= 0.25*(eta-1.);
  values[5]= 0.;
  values[6]= 0.25*(eta+1.);
  values[7]= 0.; 
}

// values of the derivatives in xi direction
static void N_Q_RT0_2D_DeriveXi(double, double, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= .25;
  values[2]= 0.;
  values[3]= .25;

   // second component
  values[4]= 0;
  values[5]= 0.;
  values[6]= 0;
  values[7]= 0.; 
}

// values of the derivatives in eta direction
static void N_Q_RT0_2D_DeriveEta(double, double, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.;
  values[2]= 0.;
  values[3]= 0.;
  
   // second component
  values[4]= .25;
  values[5]= 0;
  values[6]= .25;
  values[7]= 0; 
}

// values of derivatives in xi-xi direction
static void N_Q_RT0_2D_DeriveXiXi(double, double, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*sizeof(double)); // 2 is the space dimension
}

// values of derivatives in eta-eta direction
static void N_Q_RT0_2D_DeriveEtaEta(double, double, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*sizeof(double)); // 2 is the space dimension
}

// values of derivatives in xi-eta direction
static void N_Q_RT0_2D_DeriveXiEta(double, double, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*sizeof(double)); // 2 is the space dimension
}

static int N_Q_RT0_2D_ChangeJ0[1] = { 0 };
static int N_Q_RT0_2D_ChangeJ1[1] = { 1 };
static int N_Q_RT0_2D_ChangeJ2[1] = { 2 };
static int N_Q_RT0_2D_ChangeJ3[1] = { 3 };

static int *N_Q_RT0_2D_Change1[4] = { N_Q_RT0_2D_ChangeJ0, N_Q_RT0_2D_ChangeJ1,
                                     N_Q_RT0_2D_ChangeJ2, N_Q_RT0_2D_ChangeJ3 };
static int **N_Q_RT0_2D_Change[] = { N_Q_RT0_2D_Change1 };
