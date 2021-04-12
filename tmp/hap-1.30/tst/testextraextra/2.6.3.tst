gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/testall")[1],"25spheres.txt");;
gap> Read(dir);
gap> K:=PureCubicalComplexToCubicalComplex(M);;
gap> Size(K);
5999166
gap> L:=DVFReducedCubicalComplex(K);;
gap> D:=ChainComplex(L);;
gap> List([0..3],D!.dimension);
[ 1, 0, 25, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


