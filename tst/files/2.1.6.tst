gap> START_TEST("HAP library");
gap> B:=[ [ [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ] ],
>     [ [ 2, 2, 1 ], [ 2, 3, 1 ], [ 2, 4, 2 ], [ 2, 4, 3 ] ],
>     [ [ 4, 4, 2, 1, 3 ] ], [ ] ];;
gap> Y:=RegularCWComplex(B);
Regular CW-complex of dimension 2

gap> OrientRegularCWComplex(Y);
gap> Y!.orientation;
[ [ [ 1 ], [ 1 ], [ 1 ], [ 1 ] ], 
  [ [ 1, -1 ], [ 1, -1 ], [ 1, -1 ], [ 1, -1 ] ], [ [ 1, 1, -1, -1 ] ], [  ] ]
gap> STOP_TEST( "tst.tst", 1000 );


