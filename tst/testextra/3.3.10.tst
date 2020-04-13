gap> START_TEST("HAP library");
gap> SetAssertionLevel(0);;
gap> SetCrystGroupDefaultAction( RightAction );
gap> G:=SpaceGroup(2,10);;
gap> GroupHomology(G,3);
[ 2, 4, 4 ]
gap> STOP_TEST( "tst.tst", 1000 );


