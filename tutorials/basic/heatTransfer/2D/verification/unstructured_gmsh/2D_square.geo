//+
Point(1) = {0, 0, 0, 0.15};
//+
Point(2) = {5, 5, 0, 0.15};
//+
Point(3) = {0, 5, 0, 0.15};
//+
Point(4) = {5, 0, 0, 0.15};

Line(1) = {1, 4};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 4};
//+
Line Loop(1) = {3, 4, -1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers{1}; Recombine;
}

//+
Physical Surface("top") = {13};
//+
Physical Surface("bottom") = {21};
//+
Physical Surface("sides") = {25, 1, 17, 26};
//+
Physical Volume("body") = {1};
