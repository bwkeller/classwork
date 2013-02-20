Merge "NACA5420.stl" ;

Point(1) = {0.05, -0.05, 0.03, 0.01};
Point(2) = {-0.05, -0.05, 0.03, 0.01};
Point(3) = {-0.05, -0.05, -0.07, 0.01};
Point(4) = {0.05, -0.05, -0.07, 0.01};

Line(1)  = {1 , 2};
Line(2)  = {2 , 3};
Line(3)  = {3 , 4};
Line(4)  = {4 , 1};


Line Loop(1) = {1, 2, 3, 4};
Plane Surface(2) = {1};

Extrude{0, 0.10, 0 }{ Surface{2}; }

Surface Loop(1) = {1};
Surface Loop(2) = { 2, 13, 17, 21, 25, 26};

Volume(2) = {2,1};

Physical Volume("Fluid")  = {2};

Physical Surface("Airfoil") = {1};
Physical Surface("Floor") = {2};
Physical Surface("Inlet") = {13};
Physical Surface("Outlet") = {21};
Physical Surface("Walls") = {17, 25, 26};

