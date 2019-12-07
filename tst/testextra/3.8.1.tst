gap> START_TEST("HAP library");
gap> G:=DihedralGroup(64);; N:=Center(G);;
gap> R:=ResolutionNormalSeries([G,N],6);;
gap> C:=FilteredTensorWithIntegersModP(R,2);;
gap> P:=PersistentHomologyOfFilteredChainComplex(C,5,2);
[ [ 1, 1, 0, 0, 0, 0, 0 ], [ 0, 3, 2, 2, 2, 2, 2 ], [ 0, 0, 5, 5, 2, 2, 2 ], 
  [ 0, 0, 0, 7, 4, 4, 4 ], [ 0, 0, 0, 0, 9, 9, 4 ], [ 0, 0, 0, 0, 0, 11, 6 ], 
  [ 0, 0, 0, 0, 0, 0, 6 ] ]
gap> #BarCodeDisplay(P);;
gap> STOP_TEST( "tst.tst", 1000 );

