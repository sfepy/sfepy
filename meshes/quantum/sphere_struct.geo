Print.GeoLabels=0;
lt=0.1;
R=3; //radius
n=30; //number of boxes along the perimeter
N=20; //number of boxes along the radius (counted from the surface)

phi=Pi/2;
nphi=n*phi/(2*Pi);
An=2*Pi/n;
P=An+1;
R1=R/(P^N);

Point(1) = {R1,0,0,lt};

Point(6) = {R,0,0,lt};
Line(1) = {1,6};
Transfinite Line {1} = N Using Progression P;

Extrude Line {1, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Line {2, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi} {
  Layers { {nphi}, {1} }; Recombine;
};

Extrude Surface {9, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {5, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {26, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {43, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {60, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {77, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {94, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};
Extrude Surface {111, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/2} {
  Layers { {nphi}, {1} }; Recombine;
};

Physical Volume(100) = {1,2,3,4,5,6,7,8};
Physical Surface(1) = {136,123,89,69,103,35,55,21};
