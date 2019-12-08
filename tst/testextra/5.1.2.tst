gap> START_TEST("HAP library");
gap> A:=AbelianPcpGroup([0]);;AbelianInvariants(A);
[ 0 ]
gap> K:=EilenbergMacLaneSimplicialGroup(A,3,6);
Simplicial group of length 6

gap> C:=ChainComplexOfSimplicialGroup(K);
Chain complex of length 6 in characteristic 0 . 

gap> Homology(C,5);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
