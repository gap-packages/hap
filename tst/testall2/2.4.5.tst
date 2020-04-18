#D
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/testall")[1],"data245.txt");;
gap> Read(dir);
gap> G:=SymmetricMatrixToFilteredGraph(A,5,10);;
gap> K:=CliqueComplex(G,2);;
gap> P1:=PersistentBettiNumbers(K,1);
[ [ 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0 ], [ 0, 0, 8, 1, 0 ], [ 0, 0, 0, 9, 1 ], 
  [ 0, 0, 0, 0, 16 ] ]
gap> #BarCodeCompactDisplay(P1);
gap> STOP_TEST( "tst.tst", 1000 );


