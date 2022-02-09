gap> START_TEST("HAP library");
gap> G:=SymmetricGroup(5);;
gap> H:=AlternatingGroup(5);;
gap> R:=ResolutionFiniteGroup(G,6);;
gap> S:=ResolutionSubgroup(R,H);;
gap> Homology(TensorWithIntegers(S),5);
[ 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


