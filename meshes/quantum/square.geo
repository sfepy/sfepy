a = 50;
lt = a; // this needs to be bigger than the stuff below, so let's just use "a"
Point(1) = {-a,-a,-a, lt};
Point(2) = {-a,a,-a, lt};
Point(3) = {a, a,-a, lt};
Point(4) = {a,-a,-a, lt};
Line(1) = {1,4};
Line(2) = {4,3};
Line(3) = {3,2};
Line(4) = {2,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Physical Surface(1) = {6};

Field[1] = MathEval;
// This is generated using the formula a/5-(a/5-a/50)*(Cos(3.14*x/a)+..)/2
Field[1].F = "10-(10-0.05)*(Cos(3.14*x/50)+Cos(3.14*y/50))/2";

Background Field = 1;
