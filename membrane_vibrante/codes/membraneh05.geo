// definition du pas du maillage
h = 0.5;
// d�finition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {-2, 0, 0, h};
Point(4) = {0, 2, 0, h};
Point(5) = {0, -2, 0, h};
// d�finition des segments qui relient les points
Circle(1) = {2,1,4};
Circle(2) = {4,1,3};
Circle(3) = {3,1,5};
Circle(4) = {5,1,2};
// d�finition des contours ferm�s
Line Loop(1) = {1,2,3,4};
// d�finition des surfaces � partir contours ferm�s
Plane Surface(1) = {1};
// d�finition des �l�ments physiques : pour ces �l�ments, nous pourrons r�cup�rer
// les r�f�rences 
Physical Line(1) = {1,2,3,4};
Physical Surface(2) = {1};
