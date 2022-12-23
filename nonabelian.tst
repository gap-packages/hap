gap> START_TEST("HAP library");
gap> G:=SmallGroup(8,4);;
gap> B:=BogomolovMultiplier(G);
[  ]
gap> T:=IdGroup(Source(NonabelianTensorSquare(G).homomorphism));
[ 64, 192 ]
gap> T:=IdGroup(Source(NonabelianExteriorProduct(G,G).homomorphism));
[ 2, 1 ]
gap> WeakCommutativityGroup(G);
Pcp-group with orders [ 2, 2, 2, 2, 2, 2, 2 ]
gap> Order(SymmetricCommutativityGroup(G));
4096
gap> G:=SmallGroup(8,5);;
gap> U:=IdGroup(UpperEpicentralSeries(G,1));
[ 1, 1 ]
gap> Size(List(StemGroups(G),IdGroup));
2
gap> G:=SmallGroup(16,13);;
gap> N:=Center(G);;
gap> RelativeSchurMultiplier(G,N);
[ 2, 2 ]
gap> G:=Image(NqEpimorphismNilpotentQuotient(FreeGroup(2),2));;
gap> T:=(Source(NonabelianTensorSquare(G).homomorphism));
Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
