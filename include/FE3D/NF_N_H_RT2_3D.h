// ***********************************************************************
// Raviart-Thomas element of second order on tetrahedra, 3D
// ***********************************************************************

/* for all functionals */

//tschebyscheff point
static double tscheb_point_3 = 0.866025403784439;

static double NF_N_H_RT2_3D_Xi[] = {
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1,
  tscheb_point_3, tscheb_point_3, tscheb_point_3, 0, 0, 0, -tscheb_point_3, -tscheb_point_3, -tscheb_point_3,
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3
};
static double NF_N_H_RT2_3D_Eta[]  = {
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3
};
static double NF_N_H_RT2_3D_Zeta[] = {
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1
};

// face 0
static double NF_N_H_RT2_3D_F0_Xi[]   = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F0_Eta[]  = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F0_Zeta[] = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};

// face 1                               1
static double NF_N_H_RT2_3D_F1_Xi[]   = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F1_Eta[]  = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};
static double NF_N_H_RT2_3D_F1_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

// face 2                               2
static double NF_N_H_RT2_3D_F2_Xi[]   = { 1, 1, 1,  1, 1, 1,  1, 1, 1};
static double NF_N_H_RT2_3D_F2_Eta[]  = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F2_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

 //face 3                               3
static double NF_N_H_RT2_3D_F3_Xi[]   = { tscheb_point_3, tscheb_point_3, tscheb_point_3, 0, 0, 0, -tscheb_point_3, -tscheb_point_3, -tscheb_point_3 };
static double NF_N_H_RT2_3D_F3_Eta[]  = { 1, 1, 1,  1, 1, 1,  1, 1, 1};
static double NF_N_H_RT2_3D_F3_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

// face 4                               4
static double NF_N_H_RT2_3D_F4_Xi[]   = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};
static double NF_N_H_RT2_3D_F4_Eta[]  = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F4_Zeta[] = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };

// face 5                               5
static double NF_N_H_RT2_3D_F5_Xi[]   = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F5_Eta[]  = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F5_Zeta[] = { 1, 1, 1,  1, 1, 1,  1, 1, 1};

static double *NF_N_H_RT2_3D_XiArray[6] = {
                        NF_N_H_RT2_3D_F0_Xi,
                        NF_N_H_RT2_3D_F1_Xi,
                        NF_N_H_RT2_3D_F2_Xi,
                        NF_N_H_RT2_3D_F3_Xi,
                        NF_N_H_RT2_3D_F4_Xi,
                        NF_N_H_RT2_3D_F5_Xi };

static double *NF_N_H_RT2_3D_EtaArray[6] = {
                        NF_N_H_RT2_3D_F0_Eta,
                        NF_N_H_RT2_3D_F1_Eta,
                        NF_N_H_RT2_3D_F2_Eta,
                        NF_N_H_RT2_3D_F3_Eta,
                        NF_N_H_RT2_3D_F4_Eta,
                        NF_N_H_RT2_3D_F5_Eta };

static double *NF_N_H_RT2_3D_ZetaArray[6] = {
                        NF_N_H_RT2_3D_F0_Zeta,
                        NF_N_H_RT2_3D_F1_Zeta,
                        NF_N_H_RT2_3D_F2_Zeta,
                        NF_N_H_RT2_3D_F3_Zeta,
                        NF_N_H_RT2_3D_F4_Zeta,
                        NF_N_H_RT2_3D_F5_Zeta };

static double NF_N_H_RT2_3D_T[] = {-100};//???
static double NF_N_H_RT2_3D_S[] = {-100};//???



void NF_N_H_RT2_3D_EvalAll(const TCollection *, const TBaseCell *,
                           const double *, double *Functionals)
{
	cout << "Raviart-Thomas elements of order 2 on hexaeder: "
	     << "Nodal functionals are not fully implemented properly!" << endl;
	  for(int i=0; i<108; i++)
	    Functionals[i] = 0;
}

void NF_N_H_RT2_3D_EvalFace(const TCollection *, const TBaseCell *Cell, int face,
                            const double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 4, length == {4,4,4,4}
  Cell->GetVertex(faceVertex[4*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[4*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[4*face + 2])->GetCoords(x2,y2,z2);
  // compute measure of this face
  s = std::sqrt( std::pow((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + std::pow((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + std::pow((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;  
  Functionals[3] = PointValues[3]*s;
  Functionals[4] = PointValues[4]*s;
  Functionals[5] = PointValues[5]*s;
  Functionals[6] = PointValues[6]*s;
  Functionals[7] = PointValues[7]*s;
  Functionals[8] = PointValues[8]*s;
}

static int NF_N_H_RT2_3D_N_AllFunctionals = 108;
static int NF_N_H_RT2_3D_N_PointsAll = 54;
static int NF_N_H_RT2_3D_N_FaceFunctionals[] = { 9, 9, 9, 9, 9, 9 };
static int NF_N_H_RT2_3D_N_PointsFace[] = { 9, 9, 9, 9, 9, 9 };
