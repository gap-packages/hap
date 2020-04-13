gap> START_TEST("HAP library");
gap> A:=AbelianPcpGroup([0]);;AbelianInvariants(A);
[ 0 ]
gap> K:=EilenbergMacLaneSimplicialGroup(A,3,4);
Simplicial group of length 4

gap> C:=ChainComplexOfSimplicialGroup(K);
Chain complex of length 4 in characteristic 0 . 

gap> Homology(C,3);
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
