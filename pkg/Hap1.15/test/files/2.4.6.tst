gap> START_TEST("HAP library");
gap> Read(Concatenation(dir,"data246.txt"));
gap> F:=ThickeningFiltration(T,20);
gap> P1:=PersistentBettiNumbers(F,1);;
gap> BarCodeCompactDisplay(P1);
gap> P2:=PersistentBettiNumbers(F,1);;
gap> BarCodeCompactDisplay(P2);
gap> STOP_TEST( "tst.tst", 1000 );


