gap> START_TEST("HAP library");
gap> G:=SpaceGroup(4,1886);;
gap> R:=ResolutionBieberbachGroup(G);;
gap> C:=TensorWithIntegers(R);;
gap> Homology(C,0);
[ 0 ]
gap> Homology(C,1);
[ 2, 0 ]
gap> Homology(C,2);
[ 2 ]
gap> Homology(C,3);
[ 0 ]
gap> Homology(C,4);
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


