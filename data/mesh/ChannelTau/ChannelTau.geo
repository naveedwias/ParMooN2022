/*******************************************************************************
*                                                                              *
*      Gmsh file generating geometry and mesh for the ChannelTau example       *
*                           (gmsh ChannelTau.geo)                              *
*                                                                              *
*******************************************************************************/
// consecutive numbering of the surfaces
Geometry.OldNewReg=0;

/**************************         PARAMETERS        **************************
*                        (values to be set by the user)                       */
// for which Reynolds number (180, 395 or 590):
//Re = 180;
Re = GetValue("Reynolds number (180, 395 or 590)", 180);

// for which mesh progression (0:TANH or 1:COS):
//progression = 1;
progression = GetValue("mesh progression (0:TANH or 1:COS)", 1);

// number of refinement steps:
//n_refinement = 0;
n_refinement = GetValue("number of refinement steps", 0);

// stretching factor for TANH
If(progression==0)
  gamma = 2.75;
//  gamma = GetValue("stretching factor", 2.75);
EndIf

/**************************      MESH DEFINITION      *************************/
// number of cells per direction before refinement:
n0_x = 2;
n0_y = 2;
n0_z = 4;

// number of cells per direction after refinement:
n_x = n0_x * 2^(n_refinement);
n_y = n0_y * 2^(n_refinement);
n_z = n0_z * 2^(n_refinement);

If(progression==0)
  // For TANH progression:
  Printf("Generating channel geometry for Re: %g using TANH progression.", Re);
  For t In {0:n_z}
    position[t] = 1 + Tanh(gamma*(2.0*t/(n_z) -1))/Tanh(gamma);
    extrusion[t] = position[t];
    If(t>0)
      extrusion[t] = position[t] - position[t-1];
    EndIf
    Printf("Position z = %g ", position[t]);
  EndFor
EndIf
If(progression==1)
  // For COS progression:
  Printf("Generating channel geometry for Re: %g using COS progression.", Re);
  For t In {0:n_z}
    position[t] = 1-Cos(t*Pi/(n_z));
    extrusion[t] = position[t];
    If(t>0)
      extrusion[t] = position[t] - position[t-1];
    EndIf
    Printf("Position z = %g ", position[t]);
  EndFor
EndIf

Printf("%g refinement steps (%g cells in x-direction, %g cells in y-direction, %g cells in z-direction)", n_refinement, n_x, n_y, n_z);

/**************************    GEOMETRY DEFINITION    *************************/
If(Re==180)
  // For Re=180:
  x = 2*Pi;
  y = 2*Pi/3;
  z = 0;
  Point(1) = {-x, -y, z, 1.0};
  Point(2) = {-x, y, z, 1.0};
  Point(3) = {x, y, z, 1.0};
  Point(4) = {x, -y, z, 1.0};
EndIf
If(Re==395 || Re==590)
  // For Re=395 or 590
  x = Pi;
  y = Pi/2;
  z = 0;
  Point(1) = {-x, -y, z, 1.0};
  Point(2) = {-x, y, z, 1.0};
  Point(3) = {x, y, z, 1.0};
  Point(4) = {x, -y, z, 1.0};
EndIf

Line(1) = {2, 3};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {3, 4};
Transfinite Line {1, 2} = n_x+1 Using Progression 1;
Transfinite Line {3, 4} = n_y+1 Using Progression 1;

Line Loop(1) = {1, 4, 2, 3};
Plane Surface(1) = {1};
Transfinite Surface {1} = {2, 3, 4, 1};
Recombine Surface {1};

Z_down = 1;
Z_up = 1+6*n_z;

For t In {0:n_z-1}
  volumes[t] = t + 1;
  X_in[t] = 6+6*t;
  X_out[t] = 4+6*t;
  Y_period[2*t] = 3+6*t;
  Y_period[2*t+1] = 5+6*t;
  Extrude {0, 0, extrusion[t+1]} {
  Surface{1+6*t}; Layers{1}; Recombine;
  }
EndFor

Mesh 3;

Physical Surface(1) = {-X_in[]};       // inlet
Physical Surface(2) = {-X_out[]};      // outlet
Physical Surface(3) = {-Y_period[]};   // Y periodic faces
Physical Surface(4) = {Z_down, -Z_up}; // Z wall
Physical Volume(1) = {-volumes[]};

Mesh.Format = 30;
Mesh.SaveElementTagType = 2;

If(progression==0)
  Save Sprintf("ChannelTau%g_TANH%gx%gx%g.mesh", Re, n_x, n_y, n_z);
EndIf
If(progression==1)
  Save Sprintf("ChannelTau%g_COS%gx%gx%g.mesh", Re, n_x, n_y, n_z);
EndIf

Exit;

