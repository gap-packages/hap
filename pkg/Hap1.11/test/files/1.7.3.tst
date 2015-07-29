gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(9,1);
knot 1 with 9 crossings
gap> G:=KnotGroup(K);;
gap> L:=List(LowIndexSubgroupsFpGroup(G,6), AbelianInvariants);;
gap> SortedList(L);
[ [ 0 ], [ 0 ], [ 0, 0  ], [ 0, 0, 0, 3 ], [ 0, 0, 0, 3 ],
  [ 0, 0, 0, 3 ], [ 0,  0, 2, 3 ], [ 0, 0, 2, 3 ],
  [ 0, 0, 2, 3 ], [ 0,  0, 3 ], [ 0, 2, 2 ],
  [ 0, 2, 2, 2, 3 ], [  0, 2, 3 ], [ 0, 3, 9 ],
  [ 0, 3, 9, 9 ], [ 0,  9 ], [ 0, 9 ] ]
gap> STOP_TEST( "tst.tst", 1000 );


