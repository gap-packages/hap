#
gap> START_TEST("HAP library");
gap> G:=AbelianGroup([4,4]);;
gap> F:=CocyclicHadamardMatrices(G);;
gap> Length(F);
192
gap> STOP_TEST( "tst.tst", 1000 );
