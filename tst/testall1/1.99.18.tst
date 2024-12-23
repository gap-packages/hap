gap> START_TEST("HAP library");
gap> R:=ResolutionSL2QuadraticIntegers(-1,3);
Resolution of length 3 in characteristic 0 for matrix group . 
No contracting homotopy available. 

gap> R:=ResolutionPSL2QuadraticIntegers(-2,3);
Resolution of length 3 in characteristic 0 for PSL(2,O-2) . 
No contracting homotopy available. 

gap> R:=ResolutionSL2QuadraticIntegers(-2,2,true);;

gap> Q:=QuadraticNumberField(-2);;
gap> OQ:=RingOfIntegers(Q);;
gap> I:=QuadraticIdeal(OQ,3+1*Sqrt(-2));
ideal of norm 11 in O(Q(Sqrt(-2)))
gap> G:=HAP_CongruenceSubgroupGamma0(I);;
gap> S:=ResolutionFiniteSubgroup(R,G);;
gap> Homology(TensorWithIntegers(S),1);
[ 10, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
