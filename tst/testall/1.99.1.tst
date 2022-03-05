#
gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(3,1);                                     
prime knot 1 with 3 crossings

gap> granny:=ArcPresentation(KnotSum(K,ReflectedCubicalKnot(K))); 
[ [ 6, 9 ], [ 8, 10 ], [ 7, 9 ], [ 6, 8 ], [ 2, 7 ], [ 5, 10 ], [ 1, 3 ], 
  [ 2, 4 ], [ 3, 5 ], [ 1, 4 ] ]
gap> d1:=KnotComplementWithBoundary(granny); 
Map of regular CW-complexes

gap> Invgranny:=FirstHomologyCoveringCokernels(d1,6); 
[ [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], 
  [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], 
  [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], [ 0, 0, 0, 2 ], 
  [ 0, 0, 0, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], 
  [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2 ], [ 0, 0, 2, 2, 2 ], 
  [ 0, 0, 2, 2, 2 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], 
  [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], 
  [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], 
  [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], [ 0, 0, 2, 3 ], 
  [ 0, 0, 2, 3 ], [ 0, 0, 3, 3 ], [ 0, 0, 3, 3 ], [ 0, 2 ], [ 0, 2 ], 
  [ 0, 2 ], [ 0, 2 ], [ 0, 2 ], [ 0, 2 ], [ 0, 2 ], [ 0, 2 ], 
  [ 0, 2, 2, 2, 3 ], [ 0, 2, 2, 2, 3 ], [ 0, 2, 2, 3 ], [ 0, 2, 2, 3 ], 
  [ 0, 2, 2, 3, 3 ], [ 0, 2, 2, 3, 3 ], [ 0, 3 ], [ 0, 3 ], [ 0, 3 ], 
  [ 0, 3 ], [ 0, 3 ], [ 0, 3 ], [ 0, 3, 3, 3 ], [ 0, 3, 3, 3 ], 
  [ 0, 3, 3, 3 ], [ 0, 3, 3, 3 ], [ 0, 3, 3, 3 ], [ 0, 3, 3, 3 ], [ 0, 8 ], 
  [ 0, 8 ], [ 0, 8 ], [ 0, 8 ], [ 2, 2 ], [ 2, 2 ], [ 2, 2 ], [ 2, 2 ], 
  [ 2, 2, 2, 2, 2, 2 ], [ 2, 2, 2, 2, 2, 2 ], [ 2, 2, 2, 2, 2, 2 ], 
  [ 2, 2, 2, 2, 2, 2 ], [ 2, 2, 2, 2, 2, 4 ], [ 2, 2, 2, 2, 2, 4 ], 
  [ 2, 2, 2, 2, 2, 4 ], [ 2, 2, 2, 3, 3 ], [ 2, 2, 2, 3, 3 ], [ 2, 2, 3, 3 ], 
  [ 2, 2, 3, 3 ], [ 2, 2, 3, 3 ], [ 2, 2, 3, 3 ], [ 2, 2, 3, 5 ], 
  [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], 
  [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], 
  [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], 
  [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], 
  [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 5 ], [ 2, 2, 3, 8 ], 
  [ 2, 2, 3, 8 ], [ 2, 2, 3, 8 ], [ 2, 2, 3, 8 ], [ 2, 2, 3, 8 ], 
  [ 2, 2, 3, 8 ], [ 2, 2, 3, 8 ], [ 2, 2, 3, 8 ], [ 3, 3, 3, 9 ], 
  [ 3, 3, 3, 9 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

