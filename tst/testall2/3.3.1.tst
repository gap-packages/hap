gap> START_TEST("HAP library");
gap> R:=ResolutionFiniteGroup(SymmetricGroup(5),5);;
gap> S:=SimplifiedComplex(R);;
gap> Homology(TensorWithIntegers(S),4);
[ 2 ]
gap> Homology(TensorWithIntegers(R),4);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


