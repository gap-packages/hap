#D
gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> R:=ResolutionPrimePowerGroup(G,3);;
gap> C:=RadicalSeries(R);;
gap> P:=PersistentBettiNumbers(C,2);
[ [ 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 23, 1, 1, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 32, 1, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 33, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 33, 1, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 33, 1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 33, 1, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 27, 1, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 1, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
gap> #BarCodeCompactDisplay(P);
gap> STOP_TEST( "tst.tst", 1000 );


