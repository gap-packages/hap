gap> START_TEST("HAP library");
gap> ap:=[[2,1],[2,1]];; #Arc presentation for the trivial knot
gap> L51:=ThreeManifoldViaDehnSurgery(ap,5,1);;
gap> D51:=DijkgraafWittenInvariant(L51,CyclicGroup(5));;
gap> L52:=ThreeManifoldViaDehnSurgery(ap,5,2);;
gap> D52:=DijkgraafWittenInvariant(L52,CyclicGroup(5));;
gap> STOP_TEST( "tst.tst", 1000 );

