
MAGMAHAP:=DirectoriesPackageLibrary("HAP");
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/voronoiGL.m");

Exec(Concatenation("magma ",MAGMAHAP));
Read("/tmp/HAP_voronoi.g");
RemoveFile("/tmp/HAP_voronoi.g");

MAGMAHAP:=DirectoriesPackageLibrary("HAP");
MAGMAHAP:=Filename(MAGMAHAP,"Voronoi/ReadWellRoundedComplex.gi");
Read(MAGMAHAP);

