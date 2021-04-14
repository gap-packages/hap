gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"image2.1.1.eps");;
gap> x:=ReadImageAsWeightFunction(dir,3);;
gap> M:=x[1];
Regular CW-complex of dimension 2

gap> w:=x[2];
function( k, i ) ... end
gap> EulerIntegral(M,w);
32
gap> STOP_TEST( "tst.tst", 1000 );


