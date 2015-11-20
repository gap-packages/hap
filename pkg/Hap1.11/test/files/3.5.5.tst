gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> R:=ResolutionPrimePowerGroup(G,3);;
gap> C:=RadicalSeries(R);;
gap> P:=PersistentBettiNumbers(C,2);;
gap> BarCodeCompactDisplay(P);
gap> STOP_TEST( "tst.tst", 1000 );


