Print.GeoLabels=0;

lt_0=R/n/dens;
lt_h=R/n;
lt_R=R/n*dens;
phi=Pi/2;

Point(1) = {0,0,0,lt_0};
#Point(2) = {H,0,0,lt_h};
Point(2) = {R,0,0,lt_R};
Line(1) = {1,2};
#Line(3) = {2,3};
#Line Loop(1) = {2,3};

Extrude Line {1, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi};
Extrude Line {2, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi};

Extrude Surface {4, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {7, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};

Extrude Surface {19, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {31, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {43, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {55, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {67, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {79, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};

Physical Volume(100) = {1,2,3,4,5,6,7,8};
#Physical Surface(1) = {100,29,77,53,39,87,63,15};

