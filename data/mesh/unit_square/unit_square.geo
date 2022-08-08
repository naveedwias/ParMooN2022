//for quad meshes
Mesh.RecombineAll=1;//recombine all defined surfaces
Mesh.RecombinationAlgorithm=1;
Mesh.Algorithm=6;//delquad mesher

h = 1;
lc = h;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// boundary
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Surface(1) = {6};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
