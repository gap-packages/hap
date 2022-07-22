gap> START_TEST("HAP library");
gap> G:=Image(IsomorphismPcpGroup(SpaceGroup(2,10)));;
gap> R:=ResolutionAlmostCrystalQuotient(G,4,3);;
gap> Homology(TensorWithIntegers(R),3);
[ 2, 2, 4, 4 ]
gap> R:=ResolutionAlmostCrystalQuotient(G,4,3,true);;
gap> Homology(TensorWithIntegers(R),3);
[ 2, 2, 4, 4 ]
gap> R:=ResolutionAlmostCrystalGroup(G,4);;
gap> Homology(TensorWithIntegers(R),3);
[ 2, 4, 4 ]
gap> STOP_TEST( "tst.tst", 1000 );
