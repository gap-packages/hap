gap> START_TEST("HAP library");
gap> R:=ResolutionSL2Z(1,7);;
gap> List([1..6],i->Homology(TensorWithIntegers(R),i));
[ [ 12 ], [  ], [ 12 ], [  ], [ 12 ], [  ] ]
gap> R!.dimension(2);R!.dimension(6);
2
2
gap> R!.boundary(2,1);R!.boundary(6,1);
[ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 5 ], [ 1, 7 ], [ 1, 8 ] ]
[ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 5 ], [ 1, 7 ], [ 1, 8 ] ]
gap> R!.boundary(2,2);R!.boundary(6,2);
[ [ -2, 6 ], [ -2, 7 ], [ -1, 10 ], [ 1, 1 ], [ 1, 2 ] ]
[ [ -2, 6 ], [ -2, 7 ], [ -1, 10 ], [ 1, 1 ], [ 1, 2 ] ]
gap> STOP_TEST( "tst.tst", 1000 );


