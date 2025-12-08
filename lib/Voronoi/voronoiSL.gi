file:="~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/BasicData.m";
Read(file);

MAGMAHAP:=DirectoriesPackageLibrary("HAP");
if not d in [-1,-3] then
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/voronoiSL.m");
fi;
if d =-1 then
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/voronoiSL-1.m");
fi;
if  d=-3 then
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/voronoiSL-3.m");
fi;

Exec(Concatenation("magma ",MAGMAHAP));
Read("/tmp/HAP_voronoi.g");
RemoveFile("/tmp/HAP_voronoi.g");

MAGMAHAP:=DirectoriesPackageLibrary("HAP");
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/ReadWellRoundedComplex.gi");
Read(MAGMAHAP);

