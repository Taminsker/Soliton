xm = 0.0;
xp = 10.0;
ym = 0.0;
yp = 10.0;
hsize = 2;
// point
Point(1) = {xm, ym, 0, hsize};
Point(2) = {xp, ym, 0, hsize};
Point(3) = {xp, yp, 0, hsize};
Point(4) = {xm, yp, 0, hsize};

// build
Line(1) = {4, 3};
Line(2) = {2, 3};
Line(3) = {4, 1};
Line(4) = {1, 2};
Curve Loop(1) = {1, -2, -4, -3};
Plane Surface(1) = {1};

// mesh generation;
//Mesh.SaveAll=1;
//Mesh.ElementOrder = 2;
Mesh.MshFileVersion = 2.2;
Mesh 2;
Save "mesh.vtk";

