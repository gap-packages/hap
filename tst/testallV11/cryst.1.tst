gap> START_TEST("HAP library");
gap> G:=SpaceGroupIT(3,227);;  ##
gap> R:=ResolutionSpaceGroup(G,2);;
gap> Homology(TensorWithIntegers(R),1);
[ 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
