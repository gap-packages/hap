gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/25spheres.txt");;
gap> K:=PureCubicalComplexToCubicalComplex(M);;
gap> Size(K);
5999166
gap> L:=DVFReducedCubicalComplex(K);;
gap> D:=ChainComplex(L);;
gap> List([0..3],D!.dimension);
[ 1, 0, 25, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


