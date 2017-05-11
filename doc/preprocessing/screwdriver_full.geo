Merge "screwdriver_full.step";

//Periodic Surface 5 {7} = 26 {67};
//Periodic Surface 3 {6, 2, -6, 7} = 27 {68, 69, -68, 67};

Periodic Surface 5 {6} = 26 {66};
Periodic Surface 3 {5, 2, -5, 6} = 27 {67, 68, -67, 66};

Physical Volume(1) = {1};
Physical Volume(2) = {2};

Field[1] = MathEval;
Field[1].F = "0.0015";
Background Field = 1;
