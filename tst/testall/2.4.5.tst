#D
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/testall")[1],"data245.txt");;
gap> Read(dir);
gap> G:=SymmetricMatrixToFilteredGraph(A,5,40);;
gap> K:=CliqueComplex(G,2);;
gap> P1:=PersistentBettiNumbers(K,1);
[ [ 0, 0, 0, 0, 0 ], [ 0, 131, 2, 0, 0 ], [ 0, 0, 11, 0, 0 ], 
  [ 0, 0, 0, 2, 1 ], [ 0, 0, 0, 0, 1 ] ]
gap> #BarCodeCompactDisplay(P1);
gap> STOP_TEST( "tst.tst", 1000 );


