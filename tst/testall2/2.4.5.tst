#D
gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/data245.txt");;
gap> G:=SymmetricMatrixToFilteredGraph(A,5,10);;
gap> K:=CliqueComplex(G,2);;
gap> P1:=PersistentBettiNumbers(K,1);
[ [ 0, 0, 0, 0, 0 ], [ 0, 8, 1, 0, 0 ], [ 0, 0, 9, 1, 1 ], [ 0, 0, 0, 16, 8 ],
  [ 0, 0, 0, 0, 22 ] ]
gap> #BarCodeCompactDisplay(P1);
gap> STOP_TEST( "tst.tst", 1000 );


