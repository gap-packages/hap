#
gap> START_TEST("HAP library");
gap> A:=AbelianPcpGroup([0]);; #infinite cyclic group
gap> K:=EilenbergMacLaneSimplicialFreeAbelianGroup(A,4,7);;
gap> Homology(K,6);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
