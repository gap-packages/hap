gap> START_TEST("HAP library");
gap> M:=ReadImageAsPureCubicalComplex(Concatenation(dir,"image1.3.2.eps"),300);
Pure cubical complex of dimension 2.

gap> BettiNumber(M,0);
133
gap> STOP_TEST( "tst.tst", 1000 );


