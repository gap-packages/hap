#D
gap> START_TEST("HAP library");
gap> M:=HenonOrbit([0,0],14/10,3/10,5*10^6,500,30);
Pure cubical complex of dimension 2.

gap> Size(M);
10279
gap> #Display(M);
gap> STOP_TEST( "tst.tst", 1000 );





