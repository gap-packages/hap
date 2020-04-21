gap> START_TEST("HAP library");
gap> R:=ResolutionFiniteGroup(SymmetricGroup(5),5);;
gap> S:=SimplifiedComplex(R);;
gap> List([0..5],R!.dimension);
[ 1, 4, 10, 20, 35, 56 ]
gap> List([0..5],S!.dimension);
[ 1, 4, 7, 8, 9, 39 ]
gap> Size(R);
[ 8, 38, 100, 204, 340 ]
gap> Size(S);
[ 8, 32, 62, 83, 204 ]
gap> STOP_TEST( "tst.tst", 1000 );


