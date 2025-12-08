HAP_KOMPLEKS123:=fail;
DIMS:=0;;   #THIS IS NOT ALLOWED! CHANGE IT!!

############################################################
############################################################
VoronoiComplexGL:=function(n,d)
local file,file2,m,k,G;

if d<0 then
file:="~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/BasicData.m";
PrintTo(file, "n:=", n, "; d:=", d, "; steinitz:=1;\n");
fi;
if d>0 then
file:="~/pkg/Voronoi/VMH-TotallyRealNumberFields/Magma/BasicData.m";
PrintTo(file,"K:=QuadraticField(",d,"); n:=", n, "; V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)],j in [1..n]];"); 
fi;

file2:=DirectoriesPackageLibrary("HAP");
if d<0 then
file2:=Filename(file2,"Voronoi/voronoiGL.gi");
fi;
if d>0 then
file2:=Filename(file2,"Voronoi/voronoiRealGL.gi");
fi;

Read(file2);
G:=HAP_KOMPLEKS123!.group;
G!.bianchiInteger:=d;
SetName(G,Concatenation("GL(2,O",String(d),")"));
HAP_KOMPLEKS123!.group:=G;

for m in [0..Length(DIMS)-1] do
for k in [1..HAP_KOMPLEKS123!.dimension(m)] do
HAP_KOMPLEKS123!.boundary(m,k);
od;
od;
return HAP_KOMPLEKS123;

end;
#############################################################
#############################################################

############################################################
############################################################
VoronoiComplexSL:=function(n,d)
local file,file2,m,k,G;

if d<0 then
file:="~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/BasicData.m";
fi;
if d>0 then
Print("For a totally real quadratic number field we'd need the order of the unit group.\n");
return fail; 
file:="~/pkg/Voronoi/VMH-TotallyRealNumberFields/Magma/BasicData.m";
fi;
PrintTo(file, "n:=", n, "; d:=", d, "; steinitz:=1;\n");

file2:=DirectoriesPackageLibrary("HAP");
file2:=Filename(file2,"Voronoi/voronoiSL.gi");

Read(file2);
G:=HAP_KOMPLEKS123!.group;
G!.bianchiInteger:=d;
SetName(G,Concatenation("SL(2,O",String(d),")"));
HAP_KOMPLEKS123!.group:=G;
for m in [0..Length(DIMS)-1] do
for k in [1..HAP_KOMPLEKS123!.dimension(m)] do
HAP_KOMPLEKS123!.boundary(m,k);
od;
od;
return HAP_KOMPLEKS123;

end;
#############################################################
#############################################################

