gap> START_TEST("HAP library");
gap> R:=ResolutionFiniteGroup(SymmetricGroup(5),5);;
gap> Homology(TensorWithIntegers(R),4);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


