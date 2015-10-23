gap> START_TEST("HAP library");
gap> R:=ResolutionFiniteGroup(SymmetricGroup(5),7);;
gap> S:=SimplifiedComplex(R);;
gap> List([0..7],R!.dimension);
[ 1, 4, 10, 20, 35, 56, 76, 94 ]
gap> List([0..7],S!.dimension);
[ 1, 4, 7, 8, 9, 11, 13, 59 ]
gap> Size(R);
[ 8, 38, 100, 204, 340, 446, 648 ]
gap> Size(S);
[ 8, 32, 62, 79, 68, 75, 254 ]
gap> STOP_TEST( "tst.tst", 1000 );


