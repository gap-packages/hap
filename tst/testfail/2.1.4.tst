#E
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"image2.1.1.eps");;
gap> t:=ReadImageAsPureCubicalComplex(dir,"matrix");;
gap> t:=Int((6/10)*Maximum(Flat(t)));;
gap> M:=ReadImageAsPureCubicalComplex(dir,t);
Pure cubical complex of dimension 2.

gap> C:=ChainComplex(M);
Chain complex of length 2 in characteristic 0 . 

gap> EulerCharacteristic(C);
7
gap> STOP_TEST( "tst.tst", 1000 );


