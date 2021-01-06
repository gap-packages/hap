gap> START_TEST("HAP library");
gap> F:=FreeGroup(3);;
gap> G:=NilpotentQuotient(F,2);;
gap> R:=ResolutionNilpotentGroup(G,4);;
gap> C:=TensorWithIntegers(R);;
gap> Homology(C,3);
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


