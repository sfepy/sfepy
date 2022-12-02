SetFactory("OpenCASCADE");
// Parameters for geometry
mm = 1e-3; // conversion factor
plate_e = 10. * mm;
plate_l = 270. * mm;
plate_layers = 4;
cylinder_h = 15. * mm;
cylinder_r = 10. * mm;
cylinder_layers = 6;

// Geometry elements and plane mesh
Circle(1) = {0, 0, 0, cylinder_r, 0, 2*Pi};
Point(2) = {-plate_l/2., plate_l/2., 0};
Point(3) = {plate_l/2., plate_l/2., 0};
Point(4) = {plate_l/2., -plate_l/2., 0};
Point(5) = {-plate_l/2., -plate_l/2., 0};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 2};
Line Loop(6) = {5, 2, 3, 4};
Plane Surface(6) = {6};//+
BooleanFragments{ Surface{6}; Delete; }{ Curve{1}; Delete; }

// Extrusion towards bottom
Extrude {0, 0, -plate_e} {
  Surface{1}; Surface{2};
  Layers{plate_layers};
}
// Extrusion towards top
Extrude {0, 0, cylinder_h} {
  Surface{1}; Surface{2};
  Layers{cylinder_layers};
}

// Physical volumes
Physical Volume ("powder") = {3};
Physical Volume ("cylinder") = {4};
Physical Volume ("plate") = {1, 2};

// generate the mesh
// gmsh -3 -format mesh -o .\multi_material_cylinder_plate.mesh .\multi_material_cylinder_plate.geo

// convert to vtk format using sfepy
// sfepy-convert -d 3 .\multi_material_cylinder_plate.mesh .\multi_material_cylinder_plate.vtk

