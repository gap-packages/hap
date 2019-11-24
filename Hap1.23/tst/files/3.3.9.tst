gap> START_TEST("HAP library");
gap> S:=SylowSubgroup(SL(2,Integers mod 5^3),5);;
gap> M:=MaximalSubgroups(S);;
gap> G:=Image(IsomorphismPcGroup(M[1]));;
gap> L:=CompositionSeries(G);;
gap> R:=ResolutionSubnormalSeries(L,4);;
gap> Homology(TensorWithIntegers(R),3);
[ 5, 5, 5, 5, 5, 5, 125 ]
gap> STOP_TEST( "tst.tst", 1000 );


