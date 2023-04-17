gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"data1V2X.pdb");;
gap> K:=ReadPDBfileAsPurePermutahedralComplex(dir);
Reading chain containing 191 atoms.
Pure permutahedral complex of dimension 3.

gap> G:=KnotGroup(K);
<fp group of size infinity on the generators [ f1, f2 ]>
gap> RelatorsOfFpGroup(G);
[ f2^-1*f1*f2*f1*f2^-1*f1^-1 ]
gap> STOP_TEST( "tst.tst", 1000 );


