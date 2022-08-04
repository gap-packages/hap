#D
gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/data246.txt");;
gap> F:=ThickeningFiltration(T,5);;
gap> P1:=PersistentBettiNumbers(F,1);
[ [ 0, 0, 0, 0 ], [ 0, 364, 105, 13 ], [ 0, 0, 122, 13 ], [ 0, 0, 0, 13 ] ]
gap> #BarCodeCompactDisplay(P1);
gap> #P2:=PersistentBettiNumbers(F,1);;
gap> #BarCodeCompactDisplay(P2);
gap> STOP_TEST( "tst.tst", 1000 );


