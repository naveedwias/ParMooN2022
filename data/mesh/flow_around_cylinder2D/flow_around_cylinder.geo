//for quad meshes
Mesh.RecombineAll=1;//recombine all defined surfaces
Mesh.RecombinationAlgorithm=1;
Mesh.Algorithm=6;//delquad mesher

h = 0.3; // 0.3, 0.1, 0.05
lc = h;
lc_at_circle = lc/10.;
Point(1) = {0, 0, 0, lc};
Point(2) = {2.2, 0, 0, lc};
Point(3) = {2.2, 0.41, 0, lc};
Point(4) = {0, 0.41, 0, lc};
Point(5) = {0.2, 0.2, 0, lc_at_circle};
Point(6) = {0.15, 0.2, 0, lc_at_circle};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 6};
Line Loop(8) = {1, 2, 3, 4, -5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};

Plane Surface(8) = {8};
Physical Surface(1) = {8};

