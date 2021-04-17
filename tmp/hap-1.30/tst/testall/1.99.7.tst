#
gap> START_TEST("HAP library");
gap> D:=[[1,[2,3]],[2,[3,4]]];;
gap> R:=ResolutionCoxeterGroup(D,3);;
gap> Homology(TensorWithIntegers(R),2);
[ 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
