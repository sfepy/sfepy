//General.BackgroundImageFileName = "appartement.png";
General.BackgroundImageFileName = "";
General.BackgroundImage3D = 0; // 2D bg image (if set to 1, will draw in 3D model coordinates)
General.BackgroundImagePositionX = 10; // in pixels (or in model units in 3D case)
General.BackgroundImagePositionY = 10;
General.BackgroundImageWidth = -1; // actual size
General.BackgroundImageHeight = -1; // actual size//+
SetFactory("OpenCASCADE");

// appartment floor with walls
// Rectangle(1) = {-1.73, -0.75, 0, 2.66, 1.77, 0};
Point(1) = {-1.73, -0.75, 0, 1};
Point(2) = {-1.73+2.66, -0.75, 0, 1};
Point(3) = {-1.73+2.66, -0.75+1.77, 0, 1};
Point(4) = {-1.73, -0.75+1.77, 0, 1};


// lots of points for the walls
Point(5) = {-1.65, 0.95, 0, 1.0};
Point(6) = {-0.42, 0.95, 0, 1.0};
Point(7) = {-0.42, -0.03, 0, 1.0};
Point(8) = {-0.05, -0.04, 0, 1.0};
Point(9) = {-0.05, 0.04, 0, 1.0};
Point(10) = {-0.34, 0.04, 0, 1.0};
Point(11) = {-0.34, 0.95, 0, 1.0};
Point(12) = {0.19, 0.96, 0, 1.0};
Point(13) = {0.2, -0.03, 0, 1.0};
Point(14) = {0.28, -0.03, 0, 1.0};
Point(15) = {0.28, 0.95, 0, 1.0};
Point(16) = {0.85, 0.95, 0, 1.0};
Point(17) = {0.85, 0.04, 0, 1.0};
Point(18) = {0.52, 0.04, 0, 1.0};
Point(19) = {0.52, -0.04, 0, 1.0};
Point(20) = {0.85, -0.04, 0, 1.0};
Point(21) = {0.85, -0.19, 0, 1.0};
Point(22) = {0.52, -0.19, 0, 1.0};
Point(23) = {0.52, -0.26, 0, 1.0};
Point(24) = {0.85, -0.26, 0, 1.0};
Point(25) = {0.85, -0.68, 0, 1.0};
Point(26) = {0.07, -0.68, 0, 1.0};
Point(27) = {0.07, -0.26, 0, 1.0};
Point(28) = {0.27, -0.26, 0, 1.0};
Point(29) = {0.27, -0.19, 0, 1.0};
Point(30) = {-0.01, -0.19, 0, 1.0};
Point(31) = {-0.01, -0.68, 0, 1.0};
Point(32) = {-0.34, -0.68, 0, 1.0};
Point(33) = {-0.34, -0.19, 0, 1.0};
Point(34) = {-0.42, -0.19, 0, 1.0};
Point(35) = {-0.42, -0.68, 0, 1.0};
Point(36) = {-1.65, -0.68, 0, 1.0};

// circular source
Circle(41) = {-1.475, -0.515, 0, 0.04, 0, 2*Pi};


// dilate the geometry since it was not drawn on the right scale
Dilate {{-0.4, 0.135, 0.}, {2.63, 2.63, 2.63}} {
	Point{4, 1, 2,3 };
	Point{6, 7, 36, 8, 11, 12, 13, 9, 10, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 5};
	Curve{41}; 
}

// add lines for the outer walls
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// add lines for the inner walls
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {9, 10};
Line(9) = {13, 14};
Line(10) = {14, 15};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {15, 16};
Line(15) = {16, 17};
Line(16) = {17, 18};
Line(17) = {18, 19};
Line(18) = {19, 20};
Line(19) = {20, 21};
Line(20) = {21, 22};
Line(21) = {22, 23};
Line(22) = {23, 24};
Line(23) = {24, 25};
Line(24) = {25, 26};
Line(25) = {26, 27};
Line(26) = {27, 28};
Line(27) = {28, 29};
Line(28) = {29, 30};
Line(29) = {30, 31};
Line(30) = {31, 32};
Line(31) = {32, 33};
Line(32) = {33, 34};
Line(33) = {34, 35};
Line(34) = {35, 36};
Line(35) = {36, 5};
Line(36) = {8, 9};


// close outer walls and make a surface
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

// close inner walls and make a surface
Curve Loop(2) = {6, 7, 36, 8, 11, 12, 13, 9, 10, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 5};
Plane Surface(2) = {2};

// surface for source
Curve Loop(3) = {41};
Plane Surface(3) = {3};

// walls = outer - inner
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; }

// air = inner - source
BooleanDifference{ Surface{2}; Delete; }{ Surface{3}; }

// physical surfaces
Physical Surface("wall", 41) = {1};
Physical Surface("air", 42) = {2};
Physical Surface("source", 43) = {3};

// mesh parameters
Field[1] = MathEval;
Field[1].F = "0.5";

// finer mesh for circle
Field[2] = MathEval;
Field[2].F = "0.05";
Field[3] = Restrict;
Field[3].InField = 2;
Field[3].CurvesList = {41};

// overall background field is the minimum of other fields
Field[4] = Min;
Field[4].FieldsList = {3, 1};
Background Field = 4;
Mesh.Algorithm = 5;

// to export
// gmsh -2 helmholtz_apartment.geo -format mesh -clscale 1

// to make fit for sfepy
// sfepy-convert -2 helmholtz_apartment.mesh helmholtz_apartment.vtk