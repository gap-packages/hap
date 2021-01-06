gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(3,1);
prime knot 1 with 3 crossings

gap> G:=KnotGroup(K);;
#I  there are 2 generators and 1 relator of total length 5
gap> L:=List(LowIndexSubgroupsFpGroup(G,6), AbelianInvariants);;
gap> SortedList(L);
[ [ 0 ], [ 0 ], [ 0, 0 ], [ 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ], 
  [ 0, 0, 2 ], [ 0, 0, 2 ], [ 0, 0, 2 ], [ 0, 2 ], [ 0, 2, 2 ], 
  [ 0, 2, 2, 2 ], [ 0, 3 ], [ 0, 3 ], [ 0, 3 ], [ 0, 3, 3 ] ]
gap> STOP_TEST( "tst.tst", 1000 );


