#E
gap> START_TEST("HAP library");
gap> points:=Orbit(SymmetricGroup(4),[-2,-1,1,2],Permuted);;
gap> Y:=RegularCWPolytope(points);
Regular CW-complex of dimension 3

gap> NumbersOfCells:=List([0..3],Y!.nrCells);
[ 24, 36, 14, 1 ]
gap> STOP_TEST( "tst.tst", 1000 );


