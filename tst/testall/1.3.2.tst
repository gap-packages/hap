#E
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"image1.3.2.png");;
gap> t:=ReadImageAsPureCubicalComplex(dir,"matrix");;
gap> t:=Int((1/10)*Maximum(Flat(t)));;
gap> M:=ReadImageAsPureCubicalComplex(dir,t);
Pure cubical complex of dimension 2.

gap> BettiNumber(M,0);
20
gap> STOP_TEST( "tst.tst", 1000 );


