// mesh properties
Mesh.CharacteristicLengthMax = 100;
w = 1920;
h = 918;
stepsize = 25.;

// exterior boundary points
Point(1) = {0, 0, 0, stepsize};
Point(2) = {w, 0, 0, stepsize};
Point(3) = {w, h, 0, stepsize};
Point(4) = {0, h, 0, stepsize};

// interior boundary points
Point(5) = {126.3, 376.1, 0, stepsize};
Point(6) = {168.1, 462.9, 0, stepsize};
Point(7) = {540.1, 544.4, 0, stepsize};
Point(8) = {776.3, 686.5, 0, stepsize};
Point(9) = {1173.4, 697, 0, stepsize};
Point(10) = {1481.7, 589.3, 0, stepsize};
Point(11) = {1582.1, 531.9, 0, stepsize};
Point(12) = {1599.8, 486.9, 0, stepsize};
Point(13) = {1591.5, 390.8, 0, stepsize};
Point(14) = {1610.3, 372, 0, stepsize};
Point(15) = {1604, 336.4, 0, stepsize};
Point(16) = {1502.6, 332.2, 0, stepsize};
Point(17) = {1470.2, 303, 0, stepsize};
Point(18) = {1444.1, 200.6, 0, stepsize};
Point(19) = {1371, 189.1, 0, stepsize};
Point(20) = {1290.5, 235.1, 0, stepsize};
Point(21) = {1270.6, 298.8, 0, stepsize};
Point(22) = {1239.3, 243.4, 0, stepsize};
Point(23) = {949.8, 240.3, 0, stepsize};
Point(24) = {546.4, 236.1, 0, stepsize};
Point(25) = {500.4, 190.1, 0, stepsize};
Point(26) = {417.9, 158.8, 0, stepsize};
Point(27) = {335.3, 199.5, 0, stepsize};
Point(28) = {315.5, 256, 0, stepsize};
Point(29) = {217.2, 284.2, 0, stepsize};
Point(30) = {130.5, 328.1, 0, stepsize};

// exterior boundary lines
Line(1) = {30, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 8};

// interior boundary lines
Line(5) = {8, 9};
Line(6) = {9, 10};
Line(7) = {10, 11};
Line(8) = {11, 12};
Line(9) = {12, 13};
Line(10) = {13, 14};
Line(11) = {14, 15};
Line(12) = {15, 16};
Line(13) = {16, 17};
Line(14) = {17, 18};
Line(15) = {18, 19};
Line(16) = {19, 20};
Line(17) = {20, 21};
Line(18) = {21, 22};
Line(19) = {22, 23};
Line(20) = {23, 24};
Line(22) = {24, 25};
Line(23) = {25, 26};
Line(24) = {26, 27};
Line(25) = {27, 28};
Line(26) = {28, 29};
Line(27) = {29, 30};
Line(28) = {1, 4};
Line(29) = {4, 3};
Line(30) = {3, 2};
Line(31) = {2, 1};

// contours
Line Loop(1) = {28, 29, 30, 31};
Line Loop(2) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 1, 2};

// define surface with hole
Plane Surface(1) = {1, 2};