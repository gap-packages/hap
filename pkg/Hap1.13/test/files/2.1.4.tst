gap> START_TEST("HAP library");
gap> M:=ReadImageAsPureCubicalComplex(Concatenation(dir,"image2.1.1.eps"),600);
Pure cubical complex of dimension 2.

gap> C:=ChainComplex(M);
Chain complex of length 2 in characteristic 0 .

gap> EulerCharacteristic(C);
7
gap> STOP_TEST( "tst.tst", 1000 );


