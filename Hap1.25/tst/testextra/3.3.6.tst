gap> START_TEST("HAP library");
gap> SetAssertionLevel(0);;
gap> SetCrystGroupDefaultAction( RightAction );
gap> G:=SpaceGroup(4,2000);;
gap> R:=ResolutionCubicalCrystGroup(G,20);;
gap> Homology(TensorWithIntegers(R),19);
[ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4 ]
gap> STOP_TEST( "tst.tst", 1000 );


