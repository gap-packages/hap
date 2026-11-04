gap> START_TEST("HAP library");
gap>  Q:=QuadraticNumberField(-2);;
gap> OQ:=RingOfIntegers(Q);;
gap> I:=QuadraticIdeal(OQ,4+5*Sqrt(-2));;
gap> G:=HAP_CongruenceSubgroupGamma0(I);;
gap> IndexInSL2O(G);
144
gap> R:=ResolutionSL2QuadraticIntegers(-2,3,true);;
gap> S:=ResolutionFiniteSubgroup(R,G);;
gap> STOP_TEST( "tst.tst", 1000 );

