#E
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples")[1],"image1.3.2.png");;
gap> t:=ReadImageAsPureCubicalComplex(dir,"matrix");;
gap> t:=Int((3/10)*Maximum(Flat(t)));;
gap> M:=ReadImageAsPureCubicalComplex(dir,t);;
gap> Y:=RegularCWComplex(M);
Regular CW-complex of dimension 2

gap> Size(Y);
225206
gap> W:=ContractedComplex(Y);
Regular CW-complex of dimension 1

gap> Size(W);
1198
gap> V:=SimplifiedComplex(W);
Regular CW-complex of dimension 1

gap> Size(V);
58
gap> STOP_TEST( "tst.tst", 1000 );


