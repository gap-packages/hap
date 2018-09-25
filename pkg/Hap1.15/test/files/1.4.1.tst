gap> START_TEST("HAP library");
gap> M:=ReadImageAsPureCubicalComplex(Concatenation(dir,"image1.3.2.eps"),300);;
gap> Y:=RegularCWComplex(M);
Regular CW-complex of dimension 2
gap> Size(Y);
3074967
gap> W:=ContractedComplex(Y);
Regular CW-complex of dimension 1
gap> Size(W);
gap> V:=SimplifiedComplex(W);
Regular CW-complex of dimension 1
gap> Size(V);
3507
gap> STOP_TEST( "tst.tst", 1000 );


