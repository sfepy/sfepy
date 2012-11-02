// The fill box we are making
DX=5 ;
DY=5 ;
DZ=1 ;
centerX=0;
centerY=0;
centerZ=0;
mesh_fill=0.5; // density of nodes in the fill box

// The cylinder
xcyl=2;
ycyl=2;
radius=0.2;
mesh_cyl=0.1; // density of nodes in the cylinder


//////////////////////////////// cylindrical heating element. //////////////////////////////////////////
// The circular cross-section
For j In {0:4}
  // Calculating the points on the circle
  omega=j*2*Pi/4;
  xm[j]=(xcyl+radius*Cos(omega));
  ym[j]=(ycyl+radius*Sin(omega));
  Printf("%g %g",xm[j],ym[j]);
  Point(j)={xm[j],ym[j],-DZ/2,mesh_cyl};
EndFor
Point(5)={xcyl,ycyl,-DZ/2};

// Creating a  "Circle"
Circle(6)={0,5,1};
Circle(7)={1,5,2};
Circle(8)={2,5,3};
Circle(9)={3,5,0};

// Defining a line loop, necessary for making a surface
CircleLoop_zm=newl; Line Loop(CircleLoop_zm) = {6,7,8,9};

// Defining the surface
CircleSurf_zm=news; Plane Surface(CircleSurf_zm)={CircleLoop_zm};

// When making an extruded line or surface it is automatically meshed identically to the parent surface.
out[]=Extrude{0,0,DZ}{ Line {6,7,8,9}; };

// out contains the extruded elements. index 0,4,8,12 are the extruded lines which we now use to create the
// other end of the cylinder.
CircleLoop_zp=newl; Line Loop(CircleLoop_zp) ={out[0],out[4],out[8],out[12]};
CircleSurf_zp=news; Plane Surface(CircleSurf_zp)={CircleLoop_zp};

// The two surfaces now need to be forced to get identical nodes. This is done with the Periodic Surface command
Periodic Surface CircleSurf_zm {6,7,8,9} = CircleSurf_zp {out[0],out[4],out[8],out[12]};

//out also contain the extruded sides of the surface at index 1,5,9,13. We use them to make the entire cylinder
CylSL=newsl;Surface Loop(CylSL) = {out[1],out[5],out[9],out[13],CircleSurf_zm,CircleSurf_zp} ;
Volume(444)={CylSL};


////////////////////////////////////////////////////Fill-volume: ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
p00=newp; Point(p00)={centerX - 0.5 * DX, centerY - 0.5 * DY, centerZ - 0.5 * DZ, mesh_fill};
p01=newp; Point(p01)={centerX + 0.5 * DX, centerY - 0.5 * DY, centerZ - 0.5 * DZ, mesh_fill};
p02=newp; Point(p02)={centerX + 0.5 * DX, centerY + 0.5 * DY, centerZ - 0.5 * DZ, mesh_fill};
p03=newp; Point(p03)={centerX - 0.5 * DX, centerY + 0.5 * DY, centerZ - 0.5 * DZ, mesh_fill};
p10=newp; Point(p10)={centerX - 0.5 * DX, centerY - 0.5 * DY, centerZ + 0.5 * DZ, mesh_fill};
p11=newp; Point(p11)={centerX + 0.5 * DX, centerY - 0.5 * DY, centerZ + 0.5 * DZ, mesh_fill};
p12=newp; Point(p12)={centerX + 0.5 * DX, centerY + 0.5 * DY, centerZ + 0.5 * DZ, mesh_fill};
p13=newp; Point(p13)={centerX - 0.5 * DX, centerY + 0.5 * DY, centerZ + 0.5 * DZ, mesh_fill};

// Periphery of Fill-volume

// -z, xy plane
l00=newl; Line(l00) = {p00,p01};
l01=newl; Line(l01) = {p01,p02};
l02=newl; Line(l02) = {p02,p03};
l03=newl; Line(l03) = {p03,p00};

// +z,  xy plane
l10=newl; Line(l10) = {p10,p11};
l11=newl; Line(l11) = {p11,p12};
l12=newl; Line(l12) = {p12,p13};
l13=newl; Line(l13) = {p13,p10};

// Corner lines along z
l20=newl; Line(l20) = {p00,p10};
l21=newl; Line(l21) = {p01,p11};
l22=newl; Line(l22) = {p02,p12};
l23=newl; Line(l23) = {p03,p13};



// z-, surface. This is defined as a rectangle with a hole cut into it (CircleLoop_zm)
loop_zm=newl; Line Loop(loop_zm) = {l00,l01,l02,l03};
surf_zm=news; Plane Surface(surf_zm)={loop_zm,CircleLoop_zm};
// z+, surface This is defined as a rectangle with a hole cut into it (CircleLoop_zp)
loop_zp=newl; Line Loop(loop_zp) = {l10,l11,l12,l13};
surf_zp=news; Plane Surface(surf_zp)={loop_zp,CircleLoop_zp};
// couple z- (slave) to z+ (master) surface for periodicity in z direction
Periodic Surface surf_zm { l00,l01,l02,l03,6,7,8,9} = surf_zp { l10,l11,l12,l13,out[0],out[4],out[8],out[12]};

// y-, surface.  Here we also define periodic surfaces.
loop_ym=newl; Line Loop(loop_ym) =  {l00,l21,-l10,-l20};
surf_ym=news; Plane Surface(surf_ym)={loop_ym};
// y+, surface
loop_yp=newl; Line Loop(loop_yp) ={-l02,l22,l12,-l23};
surf_yp=news; Plane Surface(surf_yp)={loop_yp};
// couple y- (slave) to y+ (master) surface for periodicity in z direction
Periodic Surface surf_ym {l00,l21,-l10,-l20} = surf_yp{-l02,l22,l12,-l23};

// x-, surface.  Here we do not define periodic surfaces.
loop_xm=newl; Line Loop(loop_xm) =  {l20,-l13,-l23,l03};
surf_xm=news; Plane Surface(surf_xm)={loop_xm};
// x+, surface
loop_xp=newl; Line Loop(loop_xp) = {l21,l11,-l22,-l01};
surf_xp=news; Plane Surface(surf_xp)={loop_xp};

// Now we assemble the fill-volume
FillSL=newsl;Surface Loop(FillSL) = {surf_xp,surf_xm,surf_yp,surf_ym,surf_zp,surf_zm,out[1],out[5],out[9],out[13]};
Volume(555)={FillSL};

Physical Volume ("cylinder") = 444;
Physical Volume ("fill") = 555;
