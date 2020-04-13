gap> START_TEST("HAP library");
gap> SetAssertionLevel(0);;
gap> SetCrystGroupDefaultAction( RightAction );
gap> G:=SpaceGroup(4,122);;
gap> R:=ResolutionCubicalCrystGroup(G,5);;
gap> C:=HomToIntegers(R);;
gap> Cohomology(C,0);
[ 0 ]
gap> Cohomology(C,1);
[  ]
gap> Cohomology(C,2);
[ 2, 4, 4, 0 ]
gap> Cohomology(C,3);
[ 0, 0 ]
gap> Cohomology(C,4);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


