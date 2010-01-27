a = 7.5;
Point(1) = {-a,-a,-a, 0.1};
Point(2) = {-a,a,-a, 0.1};
Point(3) = {a, a,-a, 0.1};
Point(4) = {a,-a,-a, 0.1};
Line(1) = {1,4};
Line(2) = {4,3};
Line(3) = {3,2};
Line(4) = {2,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Extrude {0,0,2*a} {
  Surface{6};
}
Physical Volume(1) = {1};
