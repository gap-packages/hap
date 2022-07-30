gap> START_TEST("HAP library");
gap> K:=EilenbergMacLaneSimplicialGroup(AbelianPcpGroup([0]),2,5);;
gap> C:=ChainComplex(K);;
gap> Size(C);
20
gap> D:=CoreducedChainComplex(C);;
gap> Size(D);
10
gap> STOP_TEST( "tst.tst", 1000 );
