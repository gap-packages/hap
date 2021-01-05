gap> START_TEST("HAP library");
gap> fn:=function(i) return i^2; end;;
gap> L:=LazyList(fn);;
gap> L[4];
16
gap> Length(L);
infinity
gap> Position(L,36);
6
gap> STOP_TEST( "tst.tst", 1000 );


