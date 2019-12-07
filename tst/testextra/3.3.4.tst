gap> START_TEST("HAP library");
gap> G:=DihedralGroup(64);;
gap> Prank(G);
2
gap> R:=ResolutionFiniteGroup(G,9);;
gap> L:=List([0..8],R!.dimension);
[ 1, 2, 3, 4, 5, 6, 7, 8, 9 ]
gap> C:=TensorWithIntegersModP(R,2);;
gap> M:=List([0..8],i->Homology(C,i));;
gap> L-M;
[ 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


