#
gap> START_TEST("HAP library");
gap> A:=AbelianGroup([3]);;
gap> K:=EilenbergMacLaneSimplicialGroup(A,3,5);;
gap> C:=ChainComplex(K);;
gap> Homology(C,4);
[  ]
gap> STOP_TEST( "tst.tst", 1000 );


